#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>

namespace htrip {

/**
 * @class DijkstraOracle
 * @brief Ground-truth shortest-path oracle for weighted 4-connected grid graphs.
 *
 * @details
 * DijkstraOracle computes exact shortest paths on grid graphs featuring non-uniform edge weights,
 * discrete clearance tiers, or dynamic costmap layers via Dijkstra's priority-queue algorithm.
 *
 * Algorithmic Roles:
 * 1. Weighted Verification Oracle: Serves as the canonical ground-truth baseline in
 *    weighted exactness campaigns (`test_weighted.cpp`, `campaign_exactness.cpp`) to prove
 *    100% bit-level agreement with H-TRIP under arbitrary non-uniform edge costs.
 * 2. Small-k Fallback on Weighted Grids: When `config.hybrid_k_threshold > 0` and
 *    `grid.hasWeights() == true`, H-TRIP delegates small-k all-pairs queries to DijkstraOracle.
 * 3. Leaf Tile Construction: Used within `LeafTile::recomputeFromGrid` as an initial build
 *    subroutine when the underlying grid contains non-unit edge costs.
 *
 * Complexity Characteristics:
 * - Single-pair query: O(|E| log |V|) worst case using a binary min-heap (`std::priority_queue`).
 *   Terminates early as soon as target t is extracted from the priority queue.
 * - One-to-all query: Single Dijkstra wavefront visiting up to |V| cells, terminating early
 *   once all k target destinations have been popped and finalized.
 * - All-pairs k x k query: Performs k sequential one-to-all Dijkstra traversals,
 *   yielding O(k * |E| log |V|) total runtime.
 */
class DijkstraOracle {
public:
    /**
     * @brief Computes the exact single-pair shortest path distance between s and t.
     * @param grid The 2D grid map providing edge costs via `grid.edgeWeight()`.
     * @param s The source coordinate.
     * @param t The target coordinate.
     * @return The exact shortest path distance, or INF if unreachable.
     */
    [[nodiscard]] static dist_t queryDistance(const GridMap& grid, Point s, Point t);

    /**
     * @brief Computes single-source shortest path distances from s to multiple targets.
     * @param grid The 2D grid map providing edge costs via `grid.edgeWeight()`.
     * @param s The source coordinate.
     * @param targets Array of destination coordinates.
     * @return Vector of distances corresponding element-wise to @p targets.
     */
    [[nodiscard]] static std::vector<dist_t> queryOneToAll(
        const GridMap& grid, Point s, const std::vector<Point>& targets);

    /**
     * @brief Computes the dense all-pairs k x k distance matrix via k sequential Dijkstra sweeps.
     * @param grid The 2D grid map providing edge costs via `grid.edgeWeight()`.
     * @param viewpoints Array of k viewpoint coordinates.
     * @return Row-major k x k distance matrix where matrix[i * k + j] = dist(vp[i], vp[j]).
     */
    [[nodiscard]] static std::vector<dist_t> queryAllPairs(
        const GridMap& grid, const std::vector<Point>& viewpoints);
};

} // namespace htrip
