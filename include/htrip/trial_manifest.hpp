#pragma once

#include "robot_simulator.hpp"
#include "arena_recorder.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace htrip {

/**
 * @brief Complete, repeatable description of one H-TRIP multi-robot trial.
 *
 * Written to experiment.json so a run can be replayed from the CLI or dataset
 * browser without reconstructing knobs from memory.
 */
struct TrialManifest {
    std::string command_line;
    std::string map_path;
    std::string map_sha256;
    std::string map_id;
    int scale = 1;
    double cell_size_cm = 5.0;
    int tile_size_b = 16;
    int group_size_g = 2;
    int num_robots = 4;
    int max_frontiers = 1024;
    uint32_t seed = 0;
    int max_steps = 4000;
    double stop_coverage = 0.90;
    double dt = 0.05;
    int sim_step_multiplier = 1;
    double lidar_range_m = 8.0;
    double lidar_span_deg = 360.0;
    int lidar_num_rays = 150;
    double v_max = 0.60;
    double omega_max = 1.80;
    double radius_m = 0.18;
    int replan_period_ticks = 20;
    int persistence_window_ticks = 50;
    bool enable_async = true;
    bool enable_auto_respawn = true;
    int auto_respawn_ticks = 100;
    bool enable_live_frontier_projection = true;
    bool enable_fog_of_war = true;
    bool enable_bfs = false;
    bool enable_astar = false;
    int bfs_compare_interval = 0;
    bool auto_clock_seed_on_reset = false;
    std::string video_quality = "medium";
    int video_fps = 20;
    int video_resolution_px = 768;
    bool video_delete_ppm = true;
    int stills = 5;
    std::string log_dir;

    std::string git_commit;
    bool git_dirty = false;
    std::string hostname;
    std::string cpu_model;
    bool avx2 = false;

    std::string status = "running";
    std::string stop_reason;
    uint64_t total_ticks = 0;
    double final_coverage = 0.0;
    uint64_t bfs_compare_count = 0;
    double bfs_query_mean_us = 0.0;
    double paired_htrip_mean_us = 0.0;
    double mean_speedup_vs_bfs = 0.0;
    int grid_rows = 0;
    int grid_cols = 0;
    double world_w_m = 0.0;
    double world_h_m = 0.0;
    int tree_depth = 0;
    int total_leaves = 0;
    double wall_time_s = 0.0;
    std::string video_path;
    std::vector<std::string> still_paths;
};

/**
 * @brief Escapes control characters and quotation marks for JSON serialization.
 */
[[nodiscard]] std::string jsonEscape(const std::string& s);

/**
 * @brief Captures the current Git commit hash (HEAD) via popen.
 */
[[nodiscard]] std::string captureGitCommit();

/**
 * @brief Checks if the working git tree has uncommitted modifications.
 */
[[nodiscard]] bool captureGitDirty();

/**
 * @brief Retrieves the network hostname of the executing host.
 */
[[nodiscard]] std::string captureHostname();

/**
 * @brief Extracts the CPU model string from /proc/cpuinfo.
 */
[[nodiscard]] std::string captureCpuModel();

/**
 * @brief Computes the SHA-256 cryptographic checksum of a file on disk.
 */
[[nodiscard]] std::string sha256File(const std::string& path);

/**
 * @brief Reconstructs the complete invocation command line from argc and argv.
 */
[[nodiscard]] std::string joinCommandLine(int argc, char** argv);

/**
 * @brief Parses video quality string ("low", "medium", "high") into VideoQuality enum.
 */
[[nodiscard]] VideoQuality parseVideoQuality(const std::string& name);

/**
 * @brief Populates Git, host, and CPU hardware provenance fields into TrialManifest.
 */
void fillRuntimeEnvironment(TrialManifest& m);

/**
 * @brief Extracts map, fleet, hierarchy, and runtime metrics from active MultiRobotSimulator.
 */
void fillFromSimulator(TrialManifest& m, const MultiRobotSimulator& sim);

/**
 * @brief Serializes the complete TrialManifest into an experiment.json file.
 */
bool writeTrialManifest(const std::string& path, const TrialManifest& m);

}  // namespace htrip
