#include "htrip/internal_node.hpp"
#include "htrip/leaf_tile.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#ifndef __AVX2__
#error "AVX2 instruction set is required. Compile with -mavx2 or -march=native."
#endif
#include <immintrin.h>

namespace htrip {

InternalNode::InternalNode(int lvl, int tx_, int ty_, int b_, int gw, int gh)
    : level(lvl), tx(tx_), ty(ty_), group_w(gw), group_h(gh) {
    // INVARIANT:
    // H-TRIP hierarchy strictly requires isotropic square aggregation branching factors (gw == gh).
    assert(gw == gh && "H-TRIP InternalNode requires square branching factors (gw == gh)");

    // Determine spatial side length of child regions: child_s = b * gw^(lvl - 1)
    int child_s = b_;
    for (int l = 1; l < lvl; ++l) {
        child_s *= gw;
    }
    // Composite super-tile side length: s = child_s * gw
    s = child_s * gw;
    R0 = ty * s;
    C0 = tx * s;

    // 1. Gather child perimeter ports and record cumulative prefix offsets
    const int num_children = gw * gh;
    childOffset.resize(static_cast<size_t>(num_children + 1), 0);

    for (int k = 0; k < num_children; ++k) {
        int cy = k / gw;
        int cx = k % gw;
        int child_r0 = R0 + cy * child_s;
        int child_c0 = C0 + cx * child_s;
        std::vector<Point> cp = getTilePorts(child_r0, child_c0, child_s);
        childOffset[static_cast<size_t>(k)] = static_cast<int>(V_U.size());
        V_U.insert(V_U.end(), cp.begin(), cp.end());
    }
    childOffset[static_cast<size_t>(num_children)] = static_cast<int>(V_U.size());

    // 2. Identify exterior boundary ports of this composite super-tile (4s - 4)
    extPorts = getTilePorts(R0, C0, s);
    extPortIdx.reserve(extPorts.size());
    for (const auto& ep : extPorts) {
        int found_idx = -1;
        for (size_t i = 0; i < V_U.size(); ++i) {
            if (V_U[i] == ep) {
                found_idx = static_cast<int>(i);
                break;
            }
        }
        extPortIdx.push_back(found_idx);
    }

    // 3. Initialize super-tile APSP metric D_U to tropical identity
    const int nV = static_cast<int>(V_U.size());
    D_U.assign(static_cast<size_t>(nV * nV), INF);
    for (int i = 0; i < nV; ++i) {
        D_U[static_cast<size_t>(i * nV + i)] = 0;
    }

    // 4. Initialize exterior boundary summary S_U to tropical identity
    const int nExt = static_cast<int>(extPorts.size());
    S_U.assign(static_cast<size_t>(nExt * nExt), INF);
    for (int i = 0; i < nExt; ++i) {
        S_U[static_cast<size_t>(i * nExt + i)] = 0;
    }

    // 5. Precompute static geometric touching port pairs between sibling children.
    // Two ports u in child_a and v in child_b form a seam candidate iff they belong to
    // different children (ca != cb) and are directly 4-connected adjacent (Manhattan distance == 1).
    auto get_child = [&](int u) noexcept {
        for (int k = 0; k < num_children; ++k) {
            if (u >= childOffset[static_cast<size_t>(k)] && u < childOffset[static_cast<size_t>(k + 1)]) return k;
        }
        return -1;
    };
    for (int i = 0; i < nV; ++i) {
        int ca = get_child(i);
        for (int j = i + 1; j < nV; ++j) {
            int cb = get_child(j);
            if (ca != cb && V_U[static_cast<size_t>(i)].manhattan(V_U[static_cast<size_t>(j)]) == 1) {
                seam_candidates.push_back({i, j});
            }
        }
    }

    refreshChildToExt();
}

void InternalNode::recomputeFromChildren(
    const std::vector<std::span<const dist_t>>& child_summaries,
    const GridMap& grid) {

    const int nV = static_cast<int>(V_U.size());
    if (D_U.size() != static_cast<size_t>(nV * nV)) {
        D_U.assign(static_cast<size_t>(nV * nV), INF);
    } else {
        std::fill(D_U.begin(), D_U.end(), INF);
    }
    for (int i = 0; i < nV; ++i) {
        D_U[static_cast<size_t>(i * nV + i)] = 0;
    }

    // 1. Insert intra-child boundary summary edges into block-diagonal submatrices of D_U
    const int num_children = group_w * group_h;
    bool has_any_child_edge = false;
    for (int k = 0; k < num_children; ++k) {
        int offset = childOffset[static_cast<size_t>(k)];
        int next_offset = childOffset[static_cast<size_t>(k + 1)];
        int n_child_ports = next_offset - offset;
        const auto& summary = child_summaries[static_cast<size_t>(k)];
        if (static_cast<int>(summary.size()) < n_child_ports * n_child_ports) continue;

        // Quick inspection: does this child have ANY non-trivial connected port pairs?
        bool child_has_edges = false;
        for (int a = 0; a < n_child_ports; ++a) {
            const dist_t* row_a = &summary[static_cast<size_t>(a * n_child_ports)];
            for (int b = 0; b < n_child_ports; ++b) {
                if (a != b && row_a[b] < INF) {
                    child_has_edges = true;
                    break;
                }
            }
            if (child_has_edges) break;
        }

        if (!child_has_edges) continue;

        has_any_child_edge = true;
        for (int a = 0; a < n_child_ports; ++a) {
            const dist_t* row_a = &summary[static_cast<size_t>(a * n_child_ports)];
            int u = offset + a;
            dist_t* row_u = &D_U[static_cast<size_t>(u * nV)];
            for (int b = 0; b < n_child_ports; ++b) {
                dist_t w = row_a[b];
                if (w < INF) {
                    row_u[offset + b] = std::min(row_u[offset + b], w);
                }
            }
        }
    }

    // 2. Insert inter-child seam edges (touching ports between adjacent children)
    active_seam_scratch.clear();
    if (in_active_seam_scratch.size() < static_cast<size_t>(nV)) {
        in_active_seam_scratch.assign(static_cast<size_t>(nV), 0);
    } else {
        std::fill_n(in_active_seam_scratch.data(), nV, 0);
    }

    for (const auto& sc : seam_candidates) {
        Point pi = V_U[static_cast<size_t>(sc.u)];
        Point pj = V_U[static_cast<size_t>(sc.v)];
        if (!grid.isPassable(pi) || !grid.isPassable(pj)) continue;

        dist_t ew = grid.edgeWeight(pi, pj);
        if (ew < INF) {
            D_U[static_cast<size_t>(sc.u * nV + sc.v)] = std::min(D_U[static_cast<size_t>(sc.u * nV + sc.v)], ew);
            D_U[static_cast<size_t>(sc.v * nV + sc.u)] = std::min(D_U[static_cast<size_t>(sc.v * nV + sc.u)], ew);
            if (!in_active_seam_scratch[static_cast<size_t>(sc.u)]) {
                in_active_seam_scratch[static_cast<size_t>(sc.u)] = 1;
                active_seam_scratch.push_back(sc.u);
            }
            if (!in_active_seam_scratch[static_cast<size_t>(sc.v)]) {
                in_active_seam_scratch[static_cast<size_t>(sc.v)] = 1;
                active_seam_scratch.push_back(sc.v);
            }
        }
    }

    // Short-circuit if node is completely disconnected (all-INF fog of war)
    if (!has_any_child_edge && active_seam_scratch.empty()) {
        const int nExt = static_cast<int>(extPorts.size());
        if (S_U.size() != static_cast<size_t>(nExt * nExt)) {
            S_U.assign(static_cast<size_t>(nExt * nExt), INF);
        } else {
            std::fill(S_U.begin(), S_U.end(), INF);
        }
        for (int i = 0; i < nExt; ++i) {
            S_U[static_cast<size_t>(i * nExt + i)] = 0;
        }
        child_to_ext_dirty = true;
        refreshChildToExt();
        return;
    }

    // 3. AVX2 Min-Plus Kleene Star Closure over active seams.
    //
    // Algorithmic Rationale:
    // Standard all-pairs Floyd-Warshall over V_U would cost O(|V_U|^3).
    // Because intra-child distances are already fully closed, only the boundary seam vertices
    // active_seam_scratch can introduce new paths connecting different child regions.
    // Restricting intermediate pivots to active_seam_scratch reduces relaxation cost to
    // O(|seams| * |V_U|^2 / 16) with zero accuracy loss.
    reach_k_scratch.reserve(static_cast<size_t>(nV));

    for (int k : active_seam_scratch) {
        const dist_t* row_k = &D_U[static_cast<size_t>(k * nV)];
        reach_k_scratch.clear();
        for (int i = 0; i < nV; ++i) {
            if (row_k[i] < INF) {
                reach_k_scratch.push_back(i);
            }
        }
        if (reach_k_scratch.empty()) continue;

        for (int i : reach_k_scratch) {
            dist_t d_ik = row_k[i];
            dist_t* row_i = &D_U[static_cast<size_t>(i * nV)];
            __m256i v_dik = _mm256_set1_epi16(static_cast<short>(d_ik));
            int j = 0;
            for (; j + 16 <= nV; j += 16) {
                __m256i v_row_k = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_k[j]));
                __m256i v_row_i = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_i[j]));
                __m256i v_alt = _mm256_adds_epu16(v_dik, v_row_k);
                __m256i v_best = _mm256_min_epu16(v_row_i, v_alt);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(&row_i[j]), v_best);
            }
            // Scalar tail loop
            for (; j < nV; ++j) {
                dist_t alt = add_dist(d_ik, row_k[j]);
                if (alt < row_i[j]) {
                    row_i[j] = alt;
                }
            }
        }
    }

    // 4. Refresh exterior boundary summary S_U and child_to_ext submatrices
    const int nExt = static_cast<int>(extPorts.size());
    if (nExt > 0 && !extPortIdx.empty()) {
        for (int i = 0; i < nExt; ++i) {
            int u = extPortIdx[static_cast<size_t>(i)];
            for (int j = 0; j < nExt; ++j) {
                int v = extPortIdx[static_cast<size_t>(j)];
                S_U[static_cast<size_t>(i * nExt + j)] = D_U[static_cast<size_t>(u * nV + v)];
            }
        }
        refreshChildToExt();
    }
}

void InternalNode::updateRank1(int u, int v, dist_t w) {
    if (u == v) return;
    const int n = static_cast<int>(V_U.size());
    if (w >= D_U[static_cast<size_t>(u * n + v)]) return;
    child_to_ext_dirty = true;

    const dist_t* row_v = &D_U[static_cast<size_t>(v * n)];
    const dist_t* row_u = &D_U[static_cast<size_t>(u * n)];

    active_i_scratch.clear();
    if (active_i_scratch.capacity() < static_cast<size_t>(n)) {
        active_i_scratch.reserve(static_cast<size_t>(n));
    }
    for (int i = 0; i < n; ++i) {
        if (D_U[static_cast<size_t>(i * n + u)] < INF || D_U[static_cast<size_t>(i * n + v)] < INF) {
            active_i_scratch.push_back(i);
        }
    }

    for (int i : active_i_scratch) {
        dist_t d_iu = D_U[static_cast<size_t>(i * n + u)];
        dist_t d_iv = D_U[static_cast<size_t>(i * n + v)];

        dist_t base1 = add_dist(d_iu, w);
        dist_t base2 = add_dist(d_iv, w);
        dist_t* row_i = &D_U[static_cast<size_t>(i * n)];

        for (int j : active_i_scratch) {
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

/**
 * @brief Computes canonical perimeter port index within a local tile of side s.
 */
static int getPortIndexInTile(int r, int c, int R0, int C0, int s) noexcept {
    if (r == R0 && c >= C0 && c < C0 + s) {
        return c - C0;
    }
    if (r == R0 + s - 1 && c >= C0 && c < C0 + s) {
        return s + (c - C0);
    }
    if (c == C0 && r > R0 && r < R0 + s - 1) {
        return 2 * s + (r - R0 - 1);
    }
    if (c == C0 + s - 1 && r > R0 && r < R0 + s - 1) {
        return 3 * s - 2 + (r - R0 - 1);
    }
    return -1;
}

int InternalNode::findVertexIndex(Point p) const noexcept {
    const int child_s = s / group_w;
    if (child_s <= 0) return -1;
    const int cy = (p.r - R0) / child_s;
    const int cx = (p.c - C0) / child_s;
    if (cx < 0 || cx >= group_w || cy < 0 || cy >= group_h) return -1;

    const int quad = cy * group_w + cx;
    const int cr0 = R0 + cy * child_s;
    const int cc0 = C0 + cx * child_s;
    const int p_idx = getPortIndexInTile(p.r, p.c, cr0, cc0, child_s);
    if (p_idx < 0) return -1;
    return childOffset[static_cast<size_t>(quad)] + p_idx;
}

bool InternalNode::insertSeamEdge(Point p1, Point p2, dist_t w) {
    if (p1 == p2) return false;
    int u = findVertexIndex(p1);
    int v = findVertexIndex(p2);
    if (u != -1 && v != -1) {
        updateRank1(u, v, w);
        return true;
    }
    return false;
}

std::vector<PortDecrease> InternalNode::applyConwayPerturbations(
    const std::vector<ChildDecrease>& child_decreases,
    const std::vector<SeamEdge>& seam_edges,
    bool is_root) {

    const int n = static_cast<int>(V_U.size());
    std::vector<uint8_t> in_K(static_cast<size_t>(n), 0);
    std::vector<int> K;
    K.reserve(child_decreases.size() * 2 + seam_edges.size() * 2);

    // 1. Incorporate child port perturbations via direct edge relaxation into D_U
    for (const auto& dec : child_decreases) {
        int offset = childOffset[static_cast<size_t>(dec.quad)];
        int u = offset + dec.p_a;
        int v = offset + dec.p_b;
        if (dec.w < D_U[static_cast<size_t>(u * n + v)]) {
            D_U[static_cast<size_t>(u * n + v)] = dec.w;
            D_U[static_cast<size_t>(v * n + u)] = dec.w;
            if (!in_K[static_cast<size_t>(u)]) { in_K[static_cast<size_t>(u)] = 1; K.push_back(u); }
            if (!in_K[static_cast<size_t>(v)]) { in_K[static_cast<size_t>(v)] = 1; K.push_back(v); }
        }
    }

    // 2. Incorporate newly opened inter-child seam edges enclosed at this level
    for (const auto& se : seam_edges) {
        int u = findVertexIndex(se.p1);
        int v = findVertexIndex(se.p2);
        if (u != -1 && v != -1 && se.w < D_U[static_cast<size_t>(u * n + v)]) {
            D_U[static_cast<size_t>(u * n + v)] = se.w;
            D_U[static_cast<size_t>(v * n + u)] = se.w;
            if (!in_K[static_cast<size_t>(u)]) { in_K[static_cast<size_t>(u)] = 1; K.push_back(u); }
            if (!in_K[static_cast<size_t>(v)]) { in_K[static_cast<size_t>(v)] = 1; K.push_back(v); }
        }
    }

    if (K.empty()) {
        return {};
    }

    // 3. Batched Conway star relaxation over dirty vertices K: D'_U = D_U (min) [ D_U(:, K) + D_U(K, :) ]
    for (int k : K) {
        const dist_t* row_k = &D_U[static_cast<size_t>(k * n)];
        active_i_scratch.clear();
        for (int j = 0; j < n; ++j) {
            if (row_k[j] < INF) active_i_scratch.push_back(j);
        }
        if (active_i_scratch.empty()) continue;

        for (int i : active_i_scratch) {
            dist_t d_ik = row_k[i];
            dist_t* row_i = &D_U[static_cast<size_t>(i * n)];
            __m256i v_dik = _mm256_set1_epi16(static_cast<short>(d_ik));
            int j = 0;
            for (; j + 16 <= n; j += 16) {
                __m256i v_row_k = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_k[j]));
                __m256i v_row_i = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_i[j]));
                __m256i v_alt = _mm256_adds_epu16(v_dik, v_row_k);
                __m256i v_best = _mm256_min_epu16(v_row_i, v_alt);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(&row_i[j]), v_best);
            }
            for (; j < n; ++j) {
                dist_t alt = add_dist(d_ik, row_k[j]);
                if (alt < row_i[j]) {
                    row_i[j] = alt;
                }
            }
        }
    }

    child_to_ext_dirty = true;
    ++version;

    // 4. Update exterior boundary summary S_U and return exterior decreases
    if (!extPorts.empty() && !extPortIdx.empty()) {
        refreshChildToExt();
        auto ext_decs = refreshBoundarySummary();
        if (!is_root) {
            return ext_decs;
        }
    }
    return {};
}

std::vector<PortDecrease> InternalNode::refreshBoundarySummary() {
    const int nExt = static_cast<int>(extPorts.size());
    const int nV = static_cast<int>(V_U.size());
    std::vector<PortDecrease> decreased;
    decreased.reserve(static_cast<size_t>(nExt * 2));

    for (int i = 0; i < nExt; ++i) {
        int u = extPortIdx[static_cast<size_t>(i)];
        for (int j = 0; j < nExt; ++j) {
            int v = extPortIdx[static_cast<size_t>(j)];
            dist_t newDist = D_U[static_cast<size_t>(u * nV + v)];
            if (newDist < S_U[static_cast<size_t>(i * nExt + j)]) {
                S_U[static_cast<size_t>(i * nExt + j)] = newDist;
                decreased.push_back({i, j, newDist});
            }
        }
    }
    return decreased;
}

void InternalNode::refreshChildToExt() {
    const int nExt = static_cast<int>(extPorts.size());
    const int nV = static_cast<int>(V_U.size());
    const int num_children = group_w * group_h;
    if (child_to_ext.size() != static_cast<size_t>(num_children)) {
        child_to_ext.resize(static_cast<size_t>(num_children));
    }
    for (int q = 0; q < num_children; ++q) {
        int s_offset = childOffset[static_cast<size_t>(q)];
        int nChildPorts = childOffset[static_cast<size_t>(q + 1)] - s_offset;
        auto& mat = child_to_ext[static_cast<size_t>(q)];
        mat.resize(static_cast<size_t>(nChildPorts * nExt));
        for (int p = 0; p < nChildPorts; ++p) {
            const dist_t* row = &D_U[static_cast<size_t>((s_offset + p) * nV)];
            for (int e = 0; e < nExt; ++e) {
                mat[static_cast<size_t>(p * nExt + e)] = row[extPortIdx[static_cast<size_t>(e)]];
            }
        }
    }
    child_to_ext_dirty = false;
}

} // namespace htrip

