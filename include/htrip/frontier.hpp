#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>

namespace htrip {

/**
 * @struct FrontierCluster
 * @brief Contiguous boundary component separating known free space from unknown territory.
 *
 * @details
 * Represents an aggregated set of explorable frontier cells along with a safe,
 * clearance-aware representative viewpoint located in traversable free space.
 */
struct FrontierCluster {
    int id = 0;                         ///< Unique cluster identifier
    Point representative;               ///< Clearance-safe viewpoint coordinate p(F_i) in traversable free space
    std::vector<Point> cells;           ///< Constituent frontier grid cells forming the cluster

    /**
     * @brief Returns the number of frontier cells in this cluster.
     */
    [[nodiscard]] int size() const noexcept {
        return static_cast<int>(cells.size());
    }
};

/**
 * @class FrontierManager
 * @brief Detects, filters, and clusters explorable frontier boundaries in a grid map.
 *
 * @details
 * The FrontierManager is responsible for identifying the interface between known free space
 * and unexplored unknown territory during autonomous exploration.
 *
 * Key Algorithmic Stages:
 * 1. Frontier Classification: Classifies free cells as frontiers if they touch unknown space
 *    and satisfy robust 5x5 neighborhood depth criteria to eliminate isolated wall shadows.
 * 2. Clearance & Safety Filtering: Discards spurious frontier candidates that hug obstacle
 *    walls without opening into a substantial unknown void.
 * 3. 8-Connected Clustering: Groups contiguous frontier cells into spatial components.
 * 4. Component Partitioning: Subdivides oversized frontier lines (size > max_cluster_size)
 *    into compact local viewpoint clusters.
 * 5. Clearance-Aware Medoid Computation: Selects a representative viewpoint for each cluster
 *    using an L1 medoid objective penalized by obstacle proximity, ensuring the chosen viewpoint
 *    has safe clearance (clearance >= 2) for robot navigation.
 */
class FrontierManager {
public:
    int min_cluster_size = 2;   ///< Minimum number of cells required to form a valid cluster
    int max_cluster_size = 16;  ///< Maximum cluster size before partitioning into sub-clusters

    /**
     * @brief Constructs a FrontierManager with specified size thresholds.
     * @param min_size Minimum cell count per cluster (default: 2).
     * @param max_size Maximum cell count before sub-partitioning (default: 16).
     */
    FrontierManager(int min_size = 2, int max_size = 16)
        : min_cluster_size(min_size), max_cluster_size(max_size) {}

    /**
     * @brief Scans the entire grid map, detects valid frontier cells, and clusters them.
     * @param grid The active 2D grid map with observed cell states and clearances.
     * @return Vector of valid FrontierCluster objects with assigned representative viewpoints.
     */
    std::vector<FrontierCluster> detectAndCluster(const GridMap& grid) const;

    /**
     * @brief Tests whether a given grid cell satisfies all criteria for a valid frontier cell.
     * @details Checks 4-connected adjacency to unknown space, 5x5 unknown volume, and wall clearance.
     * @param grid The active grid map.
     * @param p Cell coordinate to test.
     * @return True if cell is a valid explorable frontier, false otherwise.
     */
    [[nodiscard]] static bool isFrontierCell(const GridMap& grid, Point p) noexcept;

    /**
     * @brief Computes a clearance-safe representative viewpoint (medoid) for a cluster.
     * @param grid The active grid map.
     * @param cluster_cells Array of coordinates belonging to the cluster.
     * @return Representative viewpoint coordinate located in safe free space.
     */
    [[nodiscard]] static Point computeMedoid(const GridMap& grid, const std::vector<Point>& cluster_cells) noexcept;
};

} // namespace htrip
