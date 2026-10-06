#include "htrip/leaf_tile.hpp"
#include <queue>
#include <algorithm>
#ifndef __AVX2__
#error "AVX2 instruction set is required. Compile with -mavx2 or -march=native."
#endif
#include <immintrin.h>

namespace htrip {

std::vector<Point> getTilePorts(int R0, int C0, int s) {
    std::vector<Point> ports;
    ports.reserve(static_cast<size_t>(4 * s - 4));

    // 1. Top perimeter row: (R0, C0 + c) for c in [0, s)
    for (int c = 0; c < s; ++c) ports.push_back({R0, C0 + c});

    // 2. Bottom perimeter row: (R0 + s - 1, C0 + c) for c in [0, s)
    for (int c = 0; c < s; ++c) ports.push_back({R0 + s - 1, C0 + c});

    // 3. Left perimeter column (excluding corners to prevent duplicate ports):
    // (R0 + r, C0) for r in [1, s - 1)
    for (int r = 1; r < s - 1; ++r) ports.push_back({R0 + r, C0});

    // 4. Right perimeter column (excluding corners):
    // (R0 + r, C0 + s - 1) for r in [1, s - 1)
    for (int r = 1; r < s - 1; ++r) ports.push_back({R0 + r, C0 + s - 1});

    return ports;
}

LeafTile::LeafTile(int tx_, int ty_, int b_) : tx(tx_), ty(ty_), b(b_) {
    R0 = ty * b;
    C0 = tx * b;
    ports = getTilePorts(R0, C0, b);

    // Intra-tile APSP matrix initialization:
    // Dimensions: (b^2) x (b^2).
    // Initialized to tropical identity: D[i, i] = 0, D[i, j] = INF for i != j.
    const int nCells = b * b;
    D.assign(static_cast<size_t>(nCells * nCells), INF);
    for (int i = 0; i < nCells; ++i) {
        D[static_cast<size_t>(i * nCells + i)] = 0;
    }

    // Boundary summary matrix initialization:
    // Dimensions: |ports| x |ports|, where |ports| = 4b - 4.
    // Initialized to tropical identity: S[i, i] = 0, S[i, j] = INF for i != j.
    const int nPorts = static_cast<int>(ports.size());
    S.assign(static_cast<size_t>(nPorts * nPorts), INF);
    for (int i = 0; i < nPorts; ++i) {
        S[static_cast<size_t>(i * nPorts + i)] = 0;
    }
}

void LeafTile::recomputeFromGrid(const GridMap& grid) {
    const int nCells = b * b;
    const int nPorts = static_cast<int>(ports.size());

    // Reset intra-tile metric matrix D
    if (D.size() != static_cast<size_t>(nCells * nCells)) {
        D.assign(static_cast<size_t>(nCells * nCells), INF);
    } else {
        std::fill(D.begin(), D.end(), INF);
    }
    for (int i = 0; i < nCells; ++i) {
        D[static_cast<size_t>(i * nCells + i)] = 0;
    }

    // Reset boundary summary matrix S
    if (S.size() != static_cast<size_t>(nPorts * nPorts)) {
        S.assign(static_cast<size_t>(nPorts * nPorts), INF);
    } else {
        std::fill(S.begin(), S.end(), INF);
    }
    for (int i = 0; i < nPorts; ++i) {
        S[static_cast<size_t>(i * nPorts + i)] = 0;
    }

    // Identify passable cells and determine whether weighted search is required
    std::vector<bool> passable(static_cast<size_t>(nCells), false);
    bool any_passable = false;
    bool has_nonunit = false;
    for (int i = 0; i < nCells; ++i) {
        int r = R0 + i / b;
        int c = C0 + i % b;
        bool p = grid.isPassable(r, c);
        passable[static_cast<size_t>(i)] = p;
        if (p) {
            any_passable = true;
            if (grid.getCellWeight(r, c) != 1) has_nonunit = true;
        }
    }
    if (!any_passable) return;

    if (!has_nonunit) {
        // Fast unweighted BFS for uniform unit-weight tiles.
        // Performs nCells BFS passes (one per passable cell).
        // Uses contiguous vector queue buffer with head pointer to avoid dynamic allocations.
        std::vector<int> q;
        q.reserve(static_cast<size_t>(nCells));

        for (int i = 0; i < nCells; ++i) {
            if (!passable[static_cast<size_t>(i)]) continue;

            dist_t* row_i = &D[static_cast<size_t>(i * nCells)];
            q.clear();
            q.push_back(i);
            size_t head = 0;

            while (head < q.size()) {
                int u = q[head++];
                dist_t d_u = row_i[u];
                dist_t alt = static_cast<dist_t>(d_u + 1);

                int ur = u / b;
                int uc = u % b;

                // 4-connected cardinal neighbors within the local b x b subgrid:
                // Up
                if (ur > 0) {
                    int v = u - b;
                    if (passable[static_cast<size_t>(v)] && row_i[v] == INF) {
                        row_i[v] = alt;
                        q.push_back(v);
                    }
                }
                // Down
                if (ur < b - 1) {
                    int v = u + b;
                    if (passable[static_cast<size_t>(v)] && row_i[v] == INF) {
                        row_i[v] = alt;
                        q.push_back(v);
                    }
                }
                // Left
                if (uc > 0) {
                    int v = u - 1;
                    if (passable[static_cast<size_t>(v)] && row_i[v] == INF) {
                        row_i[v] = alt;
                        q.push_back(v);
                    }
                }
                // Right
                if (uc < b - 1) {
                    int v = u + 1;
                    if (passable[static_cast<size_t>(v)] && row_i[v] == INF) {
                        row_i[v] = alt;
                        q.push_back(v);
                    }
                }
            }
        }
    } else {
        // Dijkstra search for non-unit weighted tiles.
        // Uses preallocated priority queue buffer hoisted outside the loop to avoid per-node allocations.
        using PQItem = std::pair<dist_t, int>;
        std::vector<PQItem> pq_buf;
        pq_buf.reserve(static_cast<size_t>(nCells * 2));
        static const int dr[4] = {-1, 1, 0, 0};
        static const int dc[4] = {0, 0, -1, 1};

        for (int i = 0; i < nCells; ++i) {
            Point src = localPoint(i);
            if (!grid.isPassable(src)) continue;

            pq_buf.clear();
            pq_buf.push_back({0, i});

            while (!pq_buf.empty()) {
                std::pop_heap(pq_buf.begin(), pq_buf.end(), std::greater<PQItem>());
                auto [d_u, u] = pq_buf.back();
                pq_buf.pop_back();
                if (d_u > D[static_cast<size_t>(i * nCells + u)]) continue;

                Point pu = localPoint(u);
                for (int k = 0; k < 4; ++k) {
                    int vr = pu.r + dr[k];
                    int vc = pu.c + dc[k];
                    if (vr < R0 || vr >= R0 + b || vc < C0 || vc >= C0 + b) continue;
                    if (!grid.isPassable(vr, vc)) continue;

                    dist_t w_edge = grid.edgeWeight(pu, {vr, vc});
                    if (w_edge >= INF) continue;

                    int v = localIdx(vr, vc);
                    dist_t alt = add_dist(d_u, w_edge);
                    if (alt < D[static_cast<size_t>(i * nCells + v)]) {
                        D[static_cast<size_t>(i * nCells + v)] = alt;
                        pq_buf.push_back({alt, v});
                        std::push_heap(pq_buf.begin(), pq_buf.end(), std::greater<PQItem>());
                    }
                }
            }
        }
    }

    // Project computed cell-to-cell shortest path metric D onto perimeter ports:
    // S[i, j] = D[localIdx(ports[i]), localIdx(ports[j])]
    S.assign(static_cast<size_t>(nPorts * nPorts), INF);
    for (int i = 0; i < nPorts; ++i) {
        int u = localIdx(ports[static_cast<size_t>(i)].r, ports[static_cast<size_t>(i)].c);
        for (int j = 0; j < nPorts; ++j) {
            int v = localIdx(ports[static_cast<size_t>(j)].r, ports[static_cast<size_t>(j)].c);
            S[static_cast<size_t>(i * nPorts + j)] = D[static_cast<size_t>(u * nCells + v)];
        }
    }
}

void LeafTile::updateRank1(int u, int v, dist_t w) {
    if (u == v) return;
    const int n = b * b;
    if (w >= D[static_cast<size_t>(u * n + v)]) return;

    const dist_t* row_v = &D[static_cast<size_t>(v * n)];
    const dist_t* row_u = &D[static_cast<size_t>(u * n)];

    // Collect active cells that have finite distance to either u or v
    active_cells_scratch.clear();
    if (active_cells_scratch.capacity() < static_cast<size_t>(n)) {
        active_cells_scratch.reserve(static_cast<size_t>(n));
    }
    for (int j = 0; j < n; ++j) {
        if (row_v[j] < INF || row_u[j] < INF) {
            active_cells_scratch.push_back(j);
        }
    }
    if (active_cells_scratch.empty()) return;

    // Symmetric rank-1 tropical outer product:
    // D'[i, j] = min(D[i, j], D[i, u] + w + D[v, j], D[i, v] + w + D[u, j])
    for (int i : active_cells_scratch) {
        dist_t d_iu = D[static_cast<size_t>(i * n + u)];
        dist_t d_iv = D[static_cast<size_t>(i * n + v)];

        dist_t base1 = add_dist(d_iu, w);
        dist_t base2 = add_dist(d_iv, w);
        dist_t* row_i = &D[static_cast<size_t>(i * n)];

        for (int j : active_cells_scratch) {
            dist_t d_vj = row_v[j];
            dist_t d_uj = row_u[j];
            dist_t alt1 = add_dist(base1, d_vj);
            dist_t alt2 = add_dist(base2, d_uj);
            dist_t best = std::min(alt1, alt2);
            if (best < row_i[j]) {
                row_i[j] = best;
            }
        }
    }
}

std::vector<PortDecrease> LeafTile::refreshBoundarySummary() {
    const int nPorts = static_cast<int>(ports.size());
    const int nCells = b * b;
    std::vector<PortDecrease> decreased;
    decreased.reserve(static_cast<size_t>(nPorts * 2));

    // Inspect projected port distances:
    // For any port pair (i, j) where D[port_i, port_j] < S[i, j], update S and emit decrease.
    for (int i = 0; i < nPorts; ++i) {
        int u = localIdx(ports[static_cast<size_t>(i)].r, ports[static_cast<size_t>(i)].c);
        for (int j = 0; j < nPorts; ++j) {
            int v = localIdx(ports[static_cast<size_t>(j)].r, ports[static_cast<size_t>(j)].c);
            dist_t newDist = D[static_cast<size_t>(u * nCells + v)];
            if (newDist < S[static_cast<size_t>(i * nPorts + j)]) {
                S[static_cast<size_t>(i * nPorts + j)] = newDist;
                decreased.push_back({i, j, newDist});
            }
        }
    }
    return decreased;
}

std::vector<PortDecrease> LeafTile::onCellOpened(Point p, const GridMap& grid) {
    const int u = localIdx(p);
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    // Insert newly formed edges between cell p and each passable 4-neighbor
    for (int k = 0; k < 4; ++k) {
        int vr = p.r + dr[k];
        int vc = p.c + dc[k];
        if (vr < R0 || vr >= R0 + b || vc < C0 || vc >= C0 + b) continue;
        if (grid.isPassable(vr, vc)) {
            dist_t w = grid.edgeWeight(p, {vr, vc});
            if (w < INF) {
                int v = localIdx(vr, vc);
                updateRank1(u, v, w);
            }
        }
    }

    return refreshBoundarySummary();
}

std::vector<PortDecrease> LeafTile::onBatchCellsOpened(const std::vector<Point>& cells, const GridMap& grid) {
    last_pivot_count = 0;
    if (cells.empty()) return {};

    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};
    const int n = b * b;

    // Fast deduplicated pivot collection
    std::vector<uint8_t> in_K(static_cast<size_t>(n), 0);
    std::vector<int> K;
    K.reserve(cells.size() * 4);

    // 1. Insert direct newly opened edges and collect affected boundary vertices K
    for (Point p : cells) {
        int u = localIdx(p);
        for (int k = 0; k < 4; ++k) {
            int vr = p.r + dr[k];
            int vc = p.c + dc[k];
            if (vr < R0 || vr >= R0 + b || vc < C0 || vc >= C0 + b) continue;
            if (grid.isPassable(vr, vc)) {
                dist_t w = grid.edgeWeight(p, {vr, vc});
                if (w < INF) {
                    int v = localIdx(vr, vc);
                    if (w < D[static_cast<size_t>(u * n + v)]) {
                        D[static_cast<size_t>(u * n + v)] = w;
                        D[static_cast<size_t>(v * n + u)] = w;
                        if (!in_K[static_cast<size_t>(u)]) { in_K[static_cast<size_t>(u)] = 1; K.push_back(u); }
                        if (!in_K[static_cast<size_t>(v)]) { in_K[static_cast<size_t>(v)] = 1; K.push_back(v); }
                    }
                }
            }
        }
    }

    last_pivot_count = static_cast<int>(K.size());
    if (K.empty()) return {};

    // 2. Batched Conway star relaxation over pivot set K with AVX2 SIMD acceleration.
    //
    // Algorithmic Structure:
    // For each pivot vertex k in K:
    //   row_k is loaded. For each row i:
    //     d_ik = D[i, k]. If d_ik < INF, we relax row_i:
    //       row_i[j] = min(row_i[j], d_ik + row_k[j]) for all j in [0, n).
    //
    // AVX2 SIMD Vectorization:
    // - _mm256_set1_epi16 broadcasts d_ik across 16 parallel 16-bit lanes.
    // - _mm256_adds_epu16 performs saturated 16-bit tropical addition: d_ik + row_k[j..j+15].
    // - _mm256_min_epu16 performs saturated 16-bit tropical minimum: min(row_i, alternative).
    // - Zero branch instructions are issued in the inner SIMD loop.
    for (int k : K) {
        const dist_t* row_k = &D[static_cast<size_t>(k * n)];
        for (int i = 0; i < n; ++i) {
            dist_t d_ik = row_k[i];
            if (d_ik >= INF) continue;
            dist_t* row_i = &D[static_cast<size_t>(i * n)];
            __m256i v_dik = _mm256_set1_epi16(static_cast<short>(d_ik));
            int j = 0;
            for (; j + 16 <= n; j += 16) {
                __m256i v_row_k = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_k[j]));
                __m256i v_row_i = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_i[j]));
                __m256i v_alt = _mm256_adds_epu16(v_dik, v_row_k);
                __m256i v_best = _mm256_min_epu16(v_row_i, v_alt);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(&row_i[j]), v_best);
            }
            // Scalar tail loop for remaining non-multiple-of-16 elements
            for (; j < n; ++j) {
                dist_t alt = add_dist(d_ik, row_k[j]);
                if (alt < row_i[j]) {
                    row_i[j] = alt;
                }
            }
        }
    }

    return refreshBoundarySummary();
}

std::vector<PortDecrease> LeafTile::onCellBlocked(Point /*p*/, const GridMap& grid) {
    // Non-monotone edge deletion: recompute intra-tile metric from grid to restore bit-exactness
    recomputeFromGrid(grid);
    return refreshBoundarySummary();
}

} // namespace htrip

