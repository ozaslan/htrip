#include "htrip/dijkstra_oracle.hpp"
#include <queue>
#include <algorithm>

namespace htrip {

// ============================================================================
// Single-Pair Dijkstra Query on Weighted Grids
// ============================================================================

dist_t DijkstraOracle::queryDistance(const GridMap& grid, Point s, Point t) {
    // 1. Boundary & passability verification
    if (!grid.isPassable(s) || !grid.isPassable(t)) return INF;
    if (s == t) return 0;

    const int R = grid.rows;
    const int C = grid.cols;
    const int total_cells = R * C;
    std::vector<dist_t> dist(static_cast<size_t>(total_cells), INF);

    const int s_idx = grid.index(s);
    const int t_idx = grid.index(t);

    // Min-heap storing pairs of (distance, cell_index)
    using PQItem = std::pair<dist_t, int>;
    std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;

    dist[static_cast<size_t>(s_idx)] = 0;
    pq.push({0, s_idx});

    // 4-connected grid neighborhood offsets: {North, South, West, East}
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    while (!pq.empty()) {
        const auto [d_u, u] = pq.top();
        pq.pop();

        // Early termination: in Dijkstra on non-negative edge graphs, the first time
        // target t is extracted from the min-heap, its tentative distance is globally optimal.
        if (u == t_idx) return dist[static_cast<size_t>(t_idx)];

        // Prune stale queue entries that have already been relaxed via shorter paths
        if (d_u > dist[static_cast<size_t>(u)]) continue;

        const int ur = u / C;
        const int uc = u % C;
        const Point pu(ur, uc);

        for (int i = 0; i < 4; ++i) {
            const int vr = ur + dr[i];
            const int vc = uc + dc[i];
            if (vr < 0 || vr >= R || vc < 0 || vc >= C) continue;
            if (!grid.isPassable(vr, vc)) continue;

            // Retrieve non-uniform edge weight or costmap penalty
            const dist_t step_w = grid.edgeWeight(pu, {vr, vc});
            if (step_w >= INF) continue;

            const int v = vr * C + vc;
            const dist_t alt = add_dist(d_u, step_w);
            if (alt < dist[static_cast<size_t>(v)]) {
                dist[static_cast<size_t>(v)] = alt;
                pq.push({alt, v});
            }
        }
    }

    return dist[static_cast<size_t>(t_idx)];
}

// ============================================================================
// One-to-All Dijkstra Query on Weighted Grids
// ============================================================================

std::vector<dist_t> DijkstraOracle::queryOneToAll(
    const GridMap& grid, Point s, const std::vector<Point>& targets) {

    const size_t k = targets.size();
    std::vector<dist_t> result(k, INF);
    if (!grid.isPassable(s) || targets.empty()) return result;

    const int R = grid.rows;
    const int C = grid.cols;
    const int total_cells = R * C;
    std::vector<dist_t> dist(static_cast<size_t>(total_cells), INF);

    const int s_idx = grid.index(s);
    using PQItem = std::pair<dist_t, int>;
    std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;

    dist[static_cast<size_t>(s_idx)] = 0;
    pq.push({0, s_idx});

    size_t targets_remaining = 0;
    std::vector<int> target_indices(k, -1);
    for (size_t i = 0; i < k; ++i) {
        if (!grid.inBounds(targets[i]) || !grid.isPassable(targets[i])) {
            result[i] = INF;
            continue;
        }
        target_indices[i] = grid.index(targets[i]);
        if (targets[i] == s) {
            result[i] = 0;
        } else {
            targets_remaining++;
        }
    }

    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    while (!pq.empty() && targets_remaining > 0) {
        const auto [d_u, u] = pq.top();
        pq.pop();

        // Prune stale queue entries
        if (d_u > dist[static_cast<size_t>(u)]) continue;

        // Check if popped node corresponds to one or more targets
        for (size_t ti = 0; ti < k; ++ti) {
            if (target_indices[ti] == u && result[ti] == INF) {
                result[ti] = d_u;
                targets_remaining--;
            }
        }
        if (targets_remaining == 0) break; // All requested targets finalized!

        const int ur = u / C;
        const int uc = u % C;
        const Point pu(ur, uc);

        for (int i = 0; i < 4; ++i) {
            const int vr = ur + dr[i];
            const int vc = uc + dc[i];
            if (vr < 0 || vr >= R || vc < 0 || vc >= C) continue;
            if (!grid.isPassable(vr, vc)) continue;

            const dist_t step_w = grid.edgeWeight(pu, {vr, vc});
            if (step_w >= INF) continue;

            const int v = vr * C + vc;
            const dist_t alt = add_dist(d_u, step_w);
            if (alt < dist[static_cast<size_t>(v)]) {
                dist[static_cast<size_t>(v)] = alt;
                pq.push({alt, v});
            }
        }
    }

    return result;
}

// ============================================================================
// All-Pairs Dijkstra Query via k Sequential Sweeps
// ============================================================================

std::vector<dist_t> DijkstraOracle::queryAllPairs(
    const GridMap& grid, const std::vector<Point>& viewpoints) {

    const size_t k = viewpoints.size();
    std::vector<dist_t> matrix(k * k, INF);
    if (k == 0) return matrix;

    // Compute each row i via a dedicated one-to-all Dijkstra search from viewpoints[i]
    for (size_t i = 0; i < k; ++i) {
        const std::vector<dist_t> row = queryOneToAll(grid, viewpoints[i], viewpoints);
        for (size_t j = 0; j < k; ++j) {
            matrix[i * k + j] = row[j];
        }
    }

    return matrix;
}

} // namespace htrip
