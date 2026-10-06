#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>

namespace htrip {

/**
 * @class BfsOracle
 * @brief Ground-truth shortest-path oracle for unweighted 4-connected grid graphs.
 *
 * @details
 * BfsOracle provides exact shortest path distances using standard Breadth-First Search
 * wavefront propagation over the 4-connected grid topology.
 *
 * Algorithmic Roles:
 * 1. Verification Oracle: Serves as the canonical ground-truth baseline in exactness
 *    verification test suites (`test_exactness.cpp`) to prove 100% bit-level agreement
 *    with H-TRIP across thousands of randomized obstacle perturbations and queries.
 * 2. Small-k Fallback: In hybrid execution mode (`config.hybrid_k_threshold > 0`),
 *    when the number of query viewpoints k is below the threshold, H-TRIP delegates
 *    all-pairs matrix queries to sequential BFS searches to eliminate hierarchy overhead.
 *
 * Complexity Characteristics:
 * - Single-pair query: O(|V|) worst case, where |V| = rows * cols. Terminates early
 *   as soon as target t is extracted from the FIFO queue.
 * - One-to-all query: Single BFS traversal visiting up to |V| cells, tracking targets
 *   via an O(1) table lookup with early exit when all k targets are reached.
 * - All-pairs k x k query: Performs k sequential one-to-all BFS traversals, resulting
 *   in O(k * |V|) total runtime.
 */
class BfsOracle {
public:
    /**
     * @brief Computes the exact single-pair shortest path distance between s and t.
     * @param grid The 2D grid map containing obstacle configurations.
     * @param s The source coordinate.
     * @param t The target coordinate.
     * @return The exact shortest path distance in unit steps, or INF if unreachable.
     */
    [[nodiscard]] static dist_t queryDistance(const GridMap& grid, Point s, Point t);

    /**
     * @brief Computes single-source shortest path distances from s to multiple targets.
     * @details Executes a single expanding BFS wavefront from source s. Utilizes an O(1)
     * direct target lookup table to detect target arrivals and terminates early as soon
     * as all valid targets have been discovered.
     * @param grid The 2D grid map containing obstacle configurations.
     * @param s The source coordinate.
     * @param targets Array of destination coordinates.
     * @return Vector of distances corresponding element-wise to @p targets.
     */
    [[nodiscard]] static std::vector<dist_t> queryOneToAll(
        const GridMap& grid, Point s, const std::vector<Point>& targets);

    /**
     * @brief Computes the dense all-pairs k x k distance matrix via k sequential BFS sweeps.
     * @param grid The 2D grid map containing obstacle configurations.
     * @param viewpoints Array of k viewpoint coordinates.
     * @return Row-major k x k distance matrix where matrix[i * k + j] = dist(vp[i], vp[j]).
     */
    [[nodiscard]] static std::vector<dist_t> queryAllPairs(
        const GridMap& grid, const std::vector<Point>& viewpoints);
};

} // namespace htrip
