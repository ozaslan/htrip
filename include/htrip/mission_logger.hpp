#pragma once

#include "types.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <utility>
#include <cstdint>
#include <mutex>

namespace htrip {

class MultiRobotSimulator;

/**
 * @brief High-resolution scientific telemetry logger for multi-robot exploration missions.
 *
 * Captures comprehensive per-step and aggregated experimental data for academic publications:
 * 1. mission_telemetry.csv: Time-series of exploration coverage, active frontiers, query/update
 *    latencies (H-TRIP vs. BFS), speedup factors, fleet tour planning time, and control frequencies.
 * 2. robot_trajectories.csv: Continuous robot SE(2) trajectories (x, y, theta, linear/angular velocity),
 *    goal commitments, waypoint progression, stuck durations, and cumulative odometry distance.
 * 3. frontier_events.csv: Frontier cluster centroids, spatial sizes, and fleet assignment matching.
 * 4. mission_summary.json: Aggregated metrics, latency distributions (mean, median, p95), total
 *    distance traveled per robot, and map/hierarchy configuration parameters.
 */
class MissionLogger {
public:
    MissionLogger() = default;
    ~MissionLogger();

    MissionLogger(const MissionLogger&) = delete;
    MissionLogger& operator=(const MissionLogger&) = delete;

    /**
     * @brief Start a telemetry logging session in the specified output directory.
     *
     * Creates the directory if it does not exist and initializes CSV file headers.
     *
     * @param output_dir Target folder for CSV and JSON files.
     * @param session_prefix Optional filename prefix or identifier.
     * @return true if files were successfully opened, false otherwise.
     */
    bool startSession(const std::string& output_dir,
                      const std::string& session_prefix = "mission_trace",
                      bool timestamp_subdir = true);

    /**
     * @brief Record telemetry snapshot for the current simulation step.
     *
     * @param sim Reference to the active MultiRobotSimulator.
     */
    void logStep(const MultiRobotSimulator& sim);

    /**
     * @brief Append one paired H-TRIP vs BFS all-pairs timing (same map snapshot and point set).
     *
     * Safe to call from the BFS worker thread. Each row is flushed immediately so a crash
     * still leaves the completed comparisons on disk.
     */
    void logBfsCompare(uint64_t job_id, uint64_t tick, size_t num_points,
                       double htrip_query_us, double bfs_query_us, double astar_query_us = 0.0);

    /**
     * @brief Finalize and flush all telemetry streams, computing mission summary statistics.
     *
     * @param sim Reference to the active MultiRobotSimulator.
     */
    void endSession(const MultiRobotSimulator& sim);

    [[nodiscard]] bool isLogging() const noexcept { return is_active; }
    [[nodiscard]] const std::string& getSessionDir() const noexcept { return session_dir; }
    [[nodiscard]] uint64_t getLoggedStepsCount() const noexcept { return logged_steps; }

private:
    bool is_active = false;
    std::string session_dir;
    uint64_t logged_steps = 0;

    std::ofstream telemetry_csv;
    std::ofstream trajectories_csv;
    std::ofstream frontiers_csv;
    std::ofstream bfs_compare_csv;
    std::mutex bfs_compare_mutex;

    // Running performance statistics for academic summary
    std::vector<double> htrip_query_latencies_us;
    std::vector<double> bfs_query_latencies_us;
    std::vector<double> htrip_update_latencies_us;
    std::vector<double> replica_lag_ticks_samples;
    std::vector<double> mtsp_plan_latencies_us;
    std::vector<double> speedup_samples;
    std::vector<double> robot_total_distance_m;
    std::vector<std::pair<double, double>> last_robot_poses;

    std::chrono::steady_clock::time_point session_start_time{};
};

} // namespace htrip
