#include "htrip/bfs_oracle.hpp"
#include <queue>

namespace htrip {

// ============================================================================
// Single-Pair BFS Query
// ============================================================================

dist_t BfsOracle::queryDistance(const GridMap& grid, Point s, Point t) {
    // 1. Boundary & passability verification
    if (!grid.isPassable(s) || !grid.isPassable(t)) return INF;
    if (s == t) return 0;

    const int R = grid.rows;
    const int C = grid.cols;
    const int total_cells = R * C;
    std::vector<dist_t> dist(static_cast<size_t>(total_cells), INF);

    const int s_idx = grid.index(s);
    const int t_idx = grid.index(t);

    // FIFO queue for unweighted BFS wavefront propagation
    std::queue<int> q;
    dist[static_cast<size_t>(s_idx)] = 0;
    q.push(s_idx);

    // 4-connected grid neighborhood offsets: {North, South, West, East}
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    while (!q.empty()) {
        const int u = q.front();
        q.pop();

        // Early termination: in an unweighted graph, the first time target t is popped,
        // its recorded distance is mathematically guaranteed to be the global minimum.
        if (u == t_idx) return dist[static_cast<size_t>(t_idx)];

        const dist_t d_u = dist[static_cast<size_t>(u)];
        const int ur = u / C;
        const int uc = u % C;

        for (int i = 0; i < 4; ++i) {
            const int vr = ur + dr[i];
            const int vc = uc + dc[i];
            if (vr < 0 || vr >= R || vc < 0 || vc >= C) continue;
            if (!grid.isPassable(vr, vc)) continue;

            const int v = vr * C + vc;
            if (dist[static_cast<size_t>(v)] == INF) {
                dist[static_cast<size_t>(v)] = add_dist(d_u, UNIT_STEP);
                q.push(v);
            }
        }
    }

    return dist[static_cast<size_t>(t_idx)];
}

// ============================================================================
// One-to-All BFS Query with Direct Target Map
// ============================================================================

namespace {

struct BfsThreadLocalScratch {
    std::vector<dist_t> dist;
    std::vector<int> target_map;
    std::vector<int> q;
    std::vector<int> visited;
    std::vector<int> target_indices;
};

thread_local BfsThreadLocalScratch tl_bfs_scratch;

} // namespace

std::vector<dist_t> BfsOracle::queryOneToAll(
    const GridMap& grid, Point s, const std::vector<Point>& targets) {

    const size_t k = targets.size();
    std::vector<dist_t> result(k, INF);
    if (!grid.isPassable(s) || targets.empty()) return result;

    const int R = grid.rows;
    const int C = grid.cols;
    const int total_cells = R * C;
    if (total_cells <= 0) return result;

    auto& scratch = tl_bfs_scratch;
    if (scratch.dist.size() != static_cast<size_t>(total_cells)) {
        scratch.dist.assign(static_cast<size_t>(total_cells), INF);
        scratch.target_map.assign(static_cast<size_t>(total_cells), -1);
    }
    scratch.visited.clear();
    scratch.q.clear();
    if (scratch.target_indices.size() < k) {
        scratch.target_indices.resize(k);
    }

    const int s_idx = grid.index(s);
    scratch.dist[static_cast<size_t>(s_idx)] = 0;
    scratch.visited.push_back(s_idx);
    scratch.q.push_back(s_idx);

    // Track targets to allow early termination if all are reached
    size_t targets_remaining = 0;
    for (size_t i = 0; i < k; ++i) {
        if (!grid.inBounds(targets[i]) || !grid.isPassable(targets[i])) {
            scratch.target_indices[i] = -1;
            result[i] = INF;
            continue;
        }
        int tidx = grid.index(targets[i]);
        scratch.target_indices[i] = tidx;
        if (targets[i] == s) {
            result[i] = 0;
        } else {
            targets_remaining++;
            if (tidx >= 0 && tidx < total_cells) {
                scratch.target_map[static_cast<size_t>(tidx)] = static_cast<int>(i);
            }
        }
    }

    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    size_t head = 0;
    while (head < scratch.q.size() && targets_remaining > 0) {
        const int u = scratch.q[head++];

        const dist_t d_u = scratch.dist[static_cast<size_t>(u)];
        const int ur = u / C;
        const int uc = u % C;

        for (int i = 0; i < 4; ++i) {
            const int vr = ur + dr[i];
            const int vc = uc + dc[i];
            if (vr < 0 || vr >= R || vc < 0 || vc >= C) continue;
            if (!grid.isPassable(vr, vc)) continue;

            const int v = vr * C + vc;
            if (scratch.dist[static_cast<size_t>(v)] == INF) {
                const dist_t alt = add_dist(d_u, UNIT_STEP);
                scratch.dist[static_cast<size_t>(v)] = alt;
                scratch.visited.push_back(v);
                scratch.q.push_back(v);

                // O(1) check if this discovered cell is one of our remaining targets
                const int tid = scratch.target_map[static_cast<size_t>(v)];
                if (tid >= 0 && result[static_cast<size_t>(tid)] == INF) {
                    result[static_cast<size_t>(tid)] = alt;
                    targets_remaining--;
                }
            }
        }
    }

    for (size_t t = 0; t < k; ++t) {
        int tidx = scratch.target_indices[t];
        if (tidx >= 0 && tidx < total_cells) {
            result[t] = scratch.dist[static_cast<size_t>(tidx)];
        }
    }

    // Clean up scratch to guarantee clean state for next traversal
    for (int v : scratch.visited) {
        scratch.dist[static_cast<size_t>(v)] = INF;
    }
    scratch.visited.clear();
    for (size_t i = 0; i < k; ++i) {
        int tidx = scratch.target_indices[i];
        if (tidx >= 0 && tidx < total_cells) {
            scratch.target_map[static_cast<size_t>(tidx)] = -1;
        }
    }

    return result;
}

// ============================================================================
// All-Pairs BFS Query via k Sequential Traversals
// ============================================================================

std::vector<dist_t> BfsOracle::queryAllPairs(
    const GridMap& grid, const std::vector<Point>& viewpoints) {

    const size_t k = viewpoints.size();
    std::vector<dist_t> matrix(k * k, INF);

    // Compute each row i via a dedicated one-to-all BFS search from viewpoints[i]
    for (size_t i = 0; i < k; ++i) {
        const auto row = queryOneToAll(grid, viewpoints[i], viewpoints);
        for (size_t j = 0; j < k; ++j) {
            matrix[i * k + j] = row[j];
        }
    }

    return matrix;
}

} // namespace htrip
