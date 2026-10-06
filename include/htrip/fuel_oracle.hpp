#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>
#include <chrono>

namespace htrip {

/**
 * @struct FuelStalenessReport
 * @brief Detailed diagnostic metrics evaluating distance cache staleness against ground truth.
 *
 * @details
 * In heuristic exploration frameworks like FUEL (Zhou et al., IEEE RA-L / ICRA 2021),
 * frontier viewpoint distance matrices are updated only within a local sensor bounding box.
 * When global shortcuts (such as opened doorways or cleared corridors) emerge outside the
 * active sensor bounding box, the cached distance between surviving distant viewpoints
 * is never updated and remains stale (an overestimation error).
 *
 * This struct quantifies the severity, prevalence, and computational cost of this staleness
 * when compared against an exact ground-truth shortest-path oracle (e.g., BfsOracle).
 */
struct FuelStalenessReport {
    size_t total_pairs = 0;              ///< Total number of distinct viewpoint pairs evaluated: k * (k - 1)
    size_t stale_pairs = 0;              ///< Number of pairs where cached_dist > ground_truth_dist
    double stale_pair_fraction = 0.0;    ///< Ratio of stale pairs: stale_pairs / total_pairs
    dist_t max_absolute_error = 0;       ///< Maximum absolute error: max(cached_dist - ground_truth_dist)
    double max_relative_error = 0.0;     ///< Maximum relative error: max((cached - gt) / gt)
    double mean_relative_error = 0.0;    ///< Mean relative error averaged across all pairs
    size_t bfs_searches_executed = 0;    ///< Number of single-source BFS wavefront searches run during update
    double update_time_us = 0.0;         ///< Wall-clock duration of the update in microseconds
};

/**
 * @class FuelCacheOracle
 * @brief Faithful model of the Incremental Frontier Information Structure (FIS) from FUEL.
 *
 * @details
 * This oracle faithfully reproduces the distance matrix maintenance algorithm published in
 * FUEL: Fast, UAV Explorations using Incremental Frontier Information Structure
 * (Zhou et al., IEEE Robotics and Automation Letters / ICRA 2021, `FrontierFinder::updateFrontierCostMatrix`).
 *
 * The FUEL Caching Protocol:
 * 1. Initial State: Computes all-pairs pairwise distances between k viewpoints via k BFS searches.
 * 2. Local Update: When an obstacle cell is opened or explored within a local sensor field-of-view (FOV)
 *    bounding box of half-width @p sensor_radius_, only viewpoints situated inside the bounding box
 *    (or newly instantiated viewpoints) are marked as modified.
 * 3. Selective Search: Single-source Dijkstra/BFS searches are executed ONLY from modified viewpoints
 *    to all other viewpoints.
 * 4. Stale Preservation: Viewpoint pairs where BOTH viewpoints reside outside the bounding box
 *    are completely bypassed; their previous distance entries are retained verbatim.
 *
 * Why Staleness Occurs:
 * If a newly opened doorway connects two previously distant wings of an environment, but both
 * viewpoints lie outside the immediate sensor bounding box, the new shortcut is completely missed
 * by the cache. The cached distance remains obsolete, leading to suboptimal TSP tour sequencing.
 *
 * H-TRIP Comparison:
 * Unlike FUEL's bounding-box heuristic, H-TRIP dynamically propagates Kleene-star min-plus updates
 * through hierarchical boundary summaries in sub-millisecond time, achieving zero staleness (100% exactness)
 * without performing any graph search.
 */
class FuelCacheOracle {
public:
    /**
     * @brief Constructs the FUEL cache oracle with a specified sensor FOV radius.
     * @param sensor_fov_radius Half-width of the square sensor bounding box (default: 12 cells).
     */
    explicit FuelCacheOracle(int sensor_fov_radius = 12)
        : sensor_radius_(sensor_fov_radius) {}

    /**
     * @brief Initializes the all-pairs distance cache from scratch via k BFS searches.
     * @param grid The 2D grid map.
     * @param viewpoints Array of k viewpoint coordinates.
     */
    void build(const GridMap& grid, const std::vector<Point>& viewpoints);

    /**
     * @brief Triggers a local FUEL cache update when a grid cell is opened.
     * @details Sets the sensor bounding box centered at @p cell_opened with radius @p sensor_radius_.
     * @param cell_opened Coordinate of the newly opened grid cell.
     * @param grid The updated grid map.
     * @param viewpoints Array of active viewpoint coordinates.
     */
    void onCellOpened(Point cell_opened, const GridMap& grid, const std::vector<Point>& viewpoints);

    /**
     * @brief Generalized sensor bounding box update centered at an arbitrary coordinate.
     * @param sensor_pos Center of the sensor field of view.
     * @param box_radius Half-width of the axis-aligned bounding box.
     * @param grid The updated grid map.
     * @param viewpoints Array of active viewpoint coordinates.
     */
    void onSensorUpdate(Point sensor_pos, int box_radius, const GridMap& grid, const std::vector<Point>& viewpoints);

    /**
     * @brief Retrieves the cached distance between viewpoint i and viewpoint j.
     * @param i Index of the source viewpoint.
     * @param j Index of the target viewpoint.
     * @return Cached distance in unit steps, or INF if out of bounds or unreachable.
     */
    [[nodiscard]] dist_t queryDistance(size_t i, size_t j) const noexcept;

    /**
     * @brief Returns a reference to the entire row-major k x k cost matrix.
     */
    [[nodiscard]] const std::vector<dist_t>& getCostMatrix() const noexcept {
        return cost_matrix_;
    }

    /**
     * @brief Returns the number of active viewpoints currently tracked.
     */
    [[nodiscard]] size_t size() const noexcept {
        return viewpoints_.size();
    }

    /**
     * @brief Evaluates the cached matrix against ground truth and generates a staleness report.
     * @param ground_truth_matrix Dense k x k ground truth matrix (e.g., from BfsOracle).
     * @return FuelStalenessReport containing error fractions, maxima, and execution times.
     */
    [[nodiscard]] FuelStalenessReport evaluateAgainstGroundTruth(
        const std::vector<dist_t>& ground_truth_matrix) const;

    /**
     * @brief Returns the configured sensor bounding-box radius.
     */
    [[nodiscard]] int sensorRadius() const noexcept { return sensor_radius_; }

    /**
     * @brief Sets the sensor bounding-box radius.
     */
    void setSensorRadius(int r) noexcept { sensor_radius_ = r; }

    /**
     * @brief Returns the duration of the most recent cache update in microseconds.
     */
    [[nodiscard]] double lastUpdateTimeUs() const noexcept { return last_update_us_; }

    /**
     * @brief Returns the number of BFS searches executed in the most recent update.
     */
    [[nodiscard]] size_t lastBfsCount() const noexcept { return last_bfs_count_; }

private:
    int sensor_radius_ = 12;
    std::vector<Point> viewpoints_;
    std::vector<dist_t> cost_matrix_; // k x k matrix
    double last_update_us_ = 0.0;
    size_t last_bfs_count_ = 0;
};

} // namespace htrip
