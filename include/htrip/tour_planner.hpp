#pragma once

#include "types.hpp"
#include <vector>
#include <cstddef>

namespace htrip {

/**
 * @struct RobotTour
 * @brief Represents an ordered sequential exploration tour for an individual robot.
 *
 * @details
 * Tracks the current target waypoint and remaining path of frontier medoids assigned to
 * a specific robot within a multi-robot exploration fleet.
 */
struct RobotTour {
    int robot_id = -1;              ///< Unique robot index within the fleet (0..N-1)
    std::vector<Point> waypoints;   ///< Ordered sequence of frontier medoids to visit
    size_t current_wp_idx = 0;      ///< Index of the active target waypoint along the tour
    uint32_t total_tour_cost = 0;   ///< Total topological travel distance of the tour in unit steps

    /**
     * @brief Checks if all waypoints along the tour have been visited or cleared.
     */
    [[nodiscard]] bool empty() const noexcept {
        return current_wp_idx >= waypoints.size();
    }

    /**
     * @brief Returns the number of unvisited waypoints remaining on this tour.
     */
    [[nodiscard]] size_t remainingWaypoints() const noexcept {
        if (current_wp_idx >= waypoints.size()) return 0;
        return waypoints.size() - current_wp_idx;
    }

    /**
     * @brief Retrieves the active goal coordinate for the robot.
     * @return Current waypoint coordinate, or {-1, -1} if tour is complete.
     */
    [[nodiscard]] Point currentGoal() const noexcept {
        if (current_wp_idx < waypoints.size()) {
            return waypoints[current_wp_idx];
        }
        return {-1, -1};
    }

    /**
     * @brief Advances to the next sequential waypoint along the tour.
     */
    void advance() noexcept {
        if (current_wp_idx < waypoints.size()) {
            current_wp_idx++;
        }
    }

    /**
     * @brief Clears all waypoints and resets tour state to empty.
     */
    void clear() {
        waypoints.clear();
        current_wp_idx = 0;
        total_tour_cost = 0;
    }
};

/**
 * @class MTSPPlanner
 * @brief High-performance Multi-Robot Traveling Salesperson (m-TSP) Planner.
 *
 * @details
 * Computes coordinated open TSP exploration routes for a fleet of N robots targeting K frontier
 * viewpoints. The planner operates directly on the exact geodesic distance matrix computed
 * by H-TRIP in microseconds, completely decoupled from map search.
 *
 * Planning Workflow:
 * 1. Capacity-Constrained Partitioning: Assigns each frontier viewpoint to the nearest reachable
 *    robot, constrained by an adaptive capacity limit to prevent single-robot workload dominance.
 * 2. Idle Robot Coverage Fill: Ensures 100% fleet utilization by allocating the least-loaded
 *    frontier to any robot that would otherwise remain idle.
 * 3. Open TSP Optimization: For each robot's assigned frontier subset, constructs an initial route
 *    via Greedy Nearest Neighbor from the robot's current pose, then iteratively refines the tour
 *    using 2-opt segment reversals to eliminate crossings.
 */
class MTSPPlanner {
public:
    /**
     * @brief Partitions K frontiers among N robots and plans an open TSP tour for each robot.
     *
     * @param robot_positions Starting coordinates of the N robots (matrix indices 0..N-1).
     * @param frontier_points Coordinates of the K frontier candidates (matrix indices N..N+K-1).
     * @param dist_matrix (N+K) x (N+K) row-major geodesic distance matrix from H-TRIP.
     * @param out_tours Output vector of size N populated with the planned tour for each robot.
     * @param max_frontiers_per_robot Maximum unique frontiers assigned to any single robot (default: 32).
     */
    static void planFleetTours(
        const std::vector<Point>& robot_positions,
        const std::vector<Point>& frontier_points,
        const std::vector<dist_t>& dist_matrix,
        std::vector<RobotTour>& out_tours,
        size_t max_frontiers_per_robot = 32);

    /**
     * @brief Computes an open TSP tour starting from a robot pose and visiting assigned frontiers.
     *
     * @param robot_idx Matrix index of the robot (0..N-1).
     * @param assigned_frontier_indices List of frontier indices (0..K-1) assigned to this robot.
     * @param num_robots Total number of robots N.
     * @param total_points Total dimension of the matrix M = N + K.
     * @param dist_matrix (N+K) x (N+K) distance matrix.
     * @return Ordered permutation of assigned_frontier_indices minimizing total open tour distance.
     */
    static std::vector<int> solveOpenTSP(
        int robot_idx,
        const std::vector<int>& assigned_frontier_indices,
        size_t num_robots,
        size_t total_points,
        const std::vector<dist_t>& dist_matrix);
};

} // namespace htrip
