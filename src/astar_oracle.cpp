#include "htrip/astar_oracle.hpp"
#include <queue>
#include <algorithm>
#include <tuple>

namespace htrip {

namespace {

/**
 * @struct PQElement
 * @brief Search node entry in the A* priority queue.
 */
struct PQElement {
    dist_t f; ///< Total estimated cost: f = g + h
    dist_t g; ///< Cost accumulated from the start node
    int u;    ///< Flattened grid cell index: u = r * cols + c

    /**
     * @brief Min-heap comparator with goal-directed tie-breaking.
     * @details If two nodes have equal f-scores, preference is given to the node
     * with the larger g-score (closer to the goal target t). This significantly
     * prunes search expansions in open corridors and planar grid regions.
     */
    bool operator>(const PQElement& o) const noexcept {
        if (f != o.f) return f > o.f;
        return g < o.g; // Tie-breaker: prefer larger g (closer to goal)
    }
};

/**
 * @brief Core A* graph search routine with epoch-stamped state reuse.
 * @param grid The 2D grid map containing obstacles and edge weights.
 * @param s Start coordinate.
 * @param t Target coordinate.
 * @param g_score Vector storing tentative shortest path costs from s.
 * @param epoch Vector storing the epoch index when each cell was last visited.
 * @param cur_epoch Active query epoch index.
 * @return Exact shortest path distance from s to t, or INF if unreachable.
 */
inline dist_t astarInternal(const GridMap& grid, Point s, Point t,
                            std::vector<dist_t>& g_score,
                            std::vector<uint32_t>& epoch,
                            uint32_t cur_epoch) {
    if (!grid.isPassable(s) || !grid.isPassable(t)) return INF;
    if (s == t) return 0;

    const int R = grid.rows;
    const int C = grid.cols;
    const int s_idx = grid.index(s);
    const int t_idx = grid.index(t);

    std::priority_queue<PQElement, std::vector<PQElement>, std::greater<PQElement>> pq;

    // Initialize start node
    epoch[static_cast<size_t>(s_idx)] = cur_epoch;
    g_score[static_cast<size_t>(s_idx)] = 0;
    const dist_t h0 = static_cast<dist_t>(s.manhattan(t) * UNIT_STEP);
    pq.push({h0, 0, s_idx});

    // 4-connected grid neighborhood offsets: {North, South, West, East}
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    while (!pq.empty()) {
        const auto [f_u, g_u, u] = pq.top();
        pq.pop();

        // Goal reached: because Manhattan distance is consistent, first expansion of t is optimal
        if (u == t_idx) return g_u;

        // Prune stale heap entries
        if (epoch[static_cast<size_t>(u)] == cur_epoch && g_u > g_score[static_cast<size_t>(u)]) {
            continue;
        }

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
            const dist_t alt_g = add_dist(g_u, step_w);

            // Relax neighbor if unvisited in the current epoch or if a strictly shorter path is found
            if (epoch[static_cast<size_t>(v)] != cur_epoch || alt_g < g_score[static_cast<size_t>(v)]) {
                epoch[static_cast<size_t>(v)] = cur_epoch;
                g_score[static_cast<size_t>(v)] = alt_g;
                const Point pv(vr, vc);
                const dist_t h = static_cast<dist_t>(pv.manhattan(t) * UNIT_STEP);
                const dist_t f = add_dist(alt_g, h);
                pq.push({f, alt_g, v});
            }
        }
    }

    return INF;
}

} // anonymous namespace

// ============================================================================
// Single-Pair A* Query
// ============================================================================

dist_t AStarOracle::queryDistance(const GridMap& grid, Point s, Point t) {
    const int total_cells = grid.rows * grid.cols;
    std::vector<dist_t> g_score(static_cast<size_t>(total_cells), INF);
    std::vector<uint32_t> epoch(static_cast<size_t>(total_cells), 0);
    return astarInternal(grid, s, t, g_score, epoch, 1);
}

// ============================================================================
// All-Pairs A* Matrix Query
// ============================================================================

std::vector<dist_t> AStarOracle::queryAllPairs(
    const GridMap& grid, const std::vector<Point>& viewpoints) {

    const size_t k = viewpoints.size();
    std::vector<dist_t> matrix(k * k, INF);
    for (size_t i = 0; i < k; ++i) {
        matrix[i * k + i] = 0;
    }

    // Allocate scratch vectors once and reuse across all k(k-1)/2 pairwise queries
    const int total_cells = grid.rows * grid.cols;
    std::vector<dist_t> g_score(static_cast<size_t>(total_cells), INF);
    std::vector<uint32_t> epoch(static_cast<size_t>(total_cells), 0);
    uint32_t cur_epoch = 0;

    // Symmetrical evaluation: compute upper triangle and mirror to lower triangle
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            cur_epoch++;
            if (cur_epoch == 0) { // 32-bit integer overflow guard: reset epoch table
                std::fill(epoch.begin(), epoch.end(), 0);
                cur_epoch = 1;
            }
            const dist_t d = astarInternal(grid, viewpoints[i], viewpoints[j], g_score, epoch, cur_epoch);
            matrix[i * k + j] = d;
            matrix[j * k + i] = d;
        }
    }

    return matrix;
}

} // namespace htrip
