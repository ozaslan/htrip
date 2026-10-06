#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>

namespace htrip {

/**
 * @class AStarOracle
 * @brief Point-to-point heuristic search baseline using the Manhattan distance heuristic.
 *
 * @details
 * AStarOracle computes exact shortest paths between pairs of points on 4-connected grid graphs
 * using the A* heuristic search algorithm.
 *
 * Algorithmic Features:
 * 1. Admissible & Consistent Heuristic: Employs the Manhattan distance metric
 *    h(u, t) = |u.r - t.r| + |u.c - t.c| (scaled by UNIT_STEP). Because edge weights
 *    are >= UNIT_STEP on 4-connected grids, Manhattan distance is provably admissible
 *    (never overestimates) and monotone consistent (satisfies triangle inequality h(u) <= c(u, v) + h(v)).
 *    Consequently, the first time a goal node is expanded, its path is guaranteed optimal.
 * 2. Goal-Directed Tie-Breaking: When multiple priority queue elements have identical
 *    f-scores (f = g + h), tie-breaking prioritizes larger g-scores (nodes closer to the goal),
 *    dramatically reducing queue expansions along flat search valleys.
 * 3. Epoch-Stamped Allocation Reuse: In all-pairs matrix computations, instead of allocating
 *    or memsetting O(|V|) state vectors for each of the k(k-1)/2 searches, the algorithm maintains
 *    an integer `epoch` vector and increments `cur_epoch`, achieving O(1) state reset per query.
 *
 * Complexity Characteristics:
 * - Single-pair query: O(|E| log |V|) worst-case, with practical runtime proportional to the
 *   area enclosed by the heuristic search ellipse.
 * - All-pairs k x k query: Performs k(k-1)/2 independent point-to-point A* searches and mirrors
 *   results across the diagonal, yielding O(k^2 * |E| log |V|) total runtime.
 */
class AStarOracle {
public:
    /**
     * @brief Computes the exact point-to-point shortest path distance between s and t via A*.
     * @param grid The 2D grid map containing obstacle configurations and edge weights.
     * @param s The source coordinate.
     * @param t The target coordinate.
     * @return The exact shortest path distance, or INF if unreachable.
     */
    [[nodiscard]] static dist_t queryDistance(const GridMap& grid, Point s, Point t);

    /**
     * @brief Computes the dense all-pairs k x k distance matrix via k(k-1)/2 pairwise A* queries.
     * @details Exploits undirected grid symmetry (dist(s, t) == dist(t, s)) to evaluate only
     * the upper triangle (j > i) and mirror values across the diagonal.
     * @param grid The 2D grid map containing obstacle configurations.
     * @param viewpoints Array of k viewpoint coordinates.
     * @return Row-major k x k distance matrix where matrix[i * k + j] = dist(vp[i], vp[j]).
     */
    [[nodiscard]] static std::vector<dist_t> queryAllPairs(
        const GridMap& grid, const std::vector<Point>& viewpoints);
};

} // namespace htrip
