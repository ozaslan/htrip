#include "htrip/mission_logger.hpp"
#include "htrip/robot_simulator.hpp"
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ctime>

namespace htrip {

MissionLogger::~MissionLogger() {
    if (is_active) {
        telemetry_csv.flush();
        trajectories_csv.flush();
        frontiers_csv.flush();
        bfs_compare_csv.flush();
    }
}

/**
 * @brief Computes the sample median of a numeric vector by partial sorting.
 */
static double calculateMedian(std::vector<double>& vals) {
    if (vals.empty()) return 0.0;
    size_t n = vals.size();
    std::sort(vals.begin(), vals.end());
    if (n % 2 == 1) return vals[n / 2];
    return (vals[n / 2 - 1] + vals[n / 2]) * 0.5;
}

/**
 * @brief Computes the specified percentile (e.g., 95th percentile) via linear interpolation.
 */
static double calculatePercentile(std::vector<double>& vals, double p) {
    if (vals.empty()) return 0.0;
    std::sort(vals.begin(), vals.end());
    double idx = (p / 100.0) * static_cast<double>(vals.size() - 1);
    size_t i0 = static_cast<size_t>(std::floor(idx));
    size_t i1 = std::min(i0 + 1, vals.size() - 1);
    double frac = idx - static_cast<double>(i0);
    return vals[i0] + frac * (vals[i1] - vals[i0]);
}

/**
 * @brief Initializes a telemetry session, creates the directory hierarchy, and writes CSV headers.
 */
bool MissionLogger::startSession(const std::string& output_dir,
                                 const std::string& session_prefix,
                                 bool timestamp_subdir) {
    std::filesystem::path dir = output_dir;
    if (timestamp_subdir) {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::tm bt{};
#if defined(_WIN32)
        localtime_s(&bt, &in_time_t);
#else
        localtime_r(&in_time_t, &bt);
#endif

        char time_str[64];
        std::strftime(time_str, sizeof(time_str), "%Y%m%d_%H%M%S", &bt);
        dir = std::filesystem::path(output_dir) / (session_prefix + "_" + time_str);
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) return false;

    session_dir = dir.string();

    telemetry_csv.open(dir / "mission_telemetry.csv");
    if (!telemetry_csv.is_open()) return false;

    trajectories_csv.open(dir / "robot_trajectories.csv");
    if (!trajectories_csv.is_open()) return false;

    frontiers_csv.open(dir / "frontier_events.csv");
    if (!frontiers_csv.is_open()) return false;

    bfs_compare_csv.open(dir / "bfs_compare.csv");
    if (!bfs_compare_csv.is_open()) return false;

    // CSV Headers
    telemetry_csv << "tick,sim_time_s,wall_time_s,num_robots,explored_cells,total_free_cells,coverage_pct,"
                  << "newly_opened_cells,active_frontiers,active_viewpoints,htrip_update_us,htrip_query_us,"
                  << "bfs_query_us,speedup,mtsp_plan_us,sim_step_freq_hz,map_update_freq_hz,query_freq_hz,"
                  << "total_respawns,replica_lag_ticks,pending_openings,updating\n";

    trajectories_csv << "tick,sim_time_s,robot_id,x_m,y_m,theta_rad,vx_mps,w_radps,has_goal,"
                     << "goal_r,goal_c,goal_x_m,goal_y_m,tour_waypoints_left,stuck_ticks,"
                     << "persistence_ticks,distance_traveled_m\n";

    frontiers_csv << "tick,sim_time_s,frontier_id,rep_r,rep_c,rep_x_m,rep_y_m,cell_count,assigned_robot_id\n";

    bfs_compare_csv << "job_id,tick,num_points,htrip_query_us,bfs_query_us,speedup,astar_query_us,astar_speedup\n";

    logged_steps = 0;
    htrip_query_latencies_us.clear();
    bfs_query_latencies_us.clear();
    htrip_update_latencies_us.clear();
    replica_lag_ticks_samples.clear();
    mtsp_plan_latencies_us.clear();
    speedup_samples.clear();
    robot_total_distance_m.clear();
    last_robot_poses.clear();

    session_start_time = std::chrono::steady_clock::now();
    is_active = true;
    return true;
}

void MissionLogger::logStep(const MultiRobotSimulator& sim) {
    if (!is_active) return;

    double sim_time_s = static_cast<double>(sim.total_ticks) * 0.05;
    auto now = std::chrono::steady_clock::now();
    double wall_time_s = std::chrono::duration<double>(now - session_start_time).count();

    size_t N = sim.robots.size();
    if (robot_total_distance_m.size() != N) {
        robot_total_distance_m.assign(N, 0.0);
        last_robot_poses.resize(N);
        for (size_t i = 0; i < N; ++i) {
            last_robot_poses[i] = {sim.robots[i].x, sim.robots[i].y};
        }
    } else {
        for (size_t i = 0; i < N; ++i) {
            double dx = sim.robots[i].x - last_robot_poses[i].first;
            double dy = sim.robots[i].y - last_robot_poses[i].second;
            double dist = std::hypot(dx, dy);
            // Ignore instantaneous relocation (respawn) jumps > 1.5m
            if (dist < 1.5) {
                robot_total_distance_m[i] += dist;
            }
            last_robot_poses[i] = {sim.robots[i].x, sim.robots[i].y};
        }
    }

    double speedup_val = sim.last_paired_speedup.load();
    if (sim.last_htrip_query_us > 0.0) {
        htrip_query_latencies_us.push_back(sim.last_htrip_query_us);
    }
    if (sim.last_htrip_update_us > 0.0) {
        htrip_update_latencies_us.push_back(sim.last_htrip_update_us);
    }
    replica_lag_ticks_samples.push_back(static_cast<double>(sim.replica_lag_ticks.load()));
    if (sim.last_mtsp_plan_us > 0.0) {
        mtsp_plan_latencies_us.push_back(sim.last_mtsp_plan_us);
    }

    // 1. Write telemetry row
    telemetry_csv << sim.total_ticks << ","
                  << std::fixed << std::setprecision(3) << sim_time_s << ","
                  << wall_time_s << ","
                  << N << ","
                  << sim.explored_free_cells << ","
                  << sim.total_free_cells_gt << ","
                  << std::setprecision(2) << sim.explorationFraction() * 100.0 << ","
                  << sim.last_revealed_cells_count << ","
                  << sim.frontiers.size() << ","
                  << sim.viewpoints.size() << ","
                  << std::setprecision(1) << sim.last_htrip_update_us << ","
                  << sim.last_htrip_query_us << ","
                  << sim.last_bfs_query_us << ","
                  << std::setprecision(2) << speedup_val << ","
                  << std::setprecision(1) << sim.last_mtsp_plan_us << ","
                  << sim.sim_step_freq_hz << ","
                  << sim.map_update_freq_hz << ","
                  << sim.htrip_query_freq_hz << ","
                  << sim.total_robot_respawns << ","
                  << sim.replica_lag_ticks.load() << ","
                  << sim.pending_openings.load() << ","
                  << (sim.is_updating_htrip.load() ? 1 : 0) << "\n";

    // 2. Write robot trajectories rows
    for (size_t i = 0; i < N; ++i) {
        const auto& bot = sim.robots[i];
        double gx_m = 0.0, gy_m = 0.0;
        if (bot.has_goal && bot.assigned_frontier.r >= 0) {
            auto [gx, gy] = sim.cellToWorld(bot.assigned_frontier);
            gx_m = gx;
            gy_m = gy;
        }

        trajectories_csv << sim.total_ticks << ","
                         << std::fixed << std::setprecision(3) << sim_time_s << ","
                         << i << ","
                         << std::setprecision(4) << bot.x << ","
                         << bot.y << ","
                         << bot.theta << ","
                         << std::setprecision(3) << bot.v << ","
                         << bot.omega << ","
                         << (bot.has_goal ? 1 : 0) << ","
                         << bot.assigned_frontier.r << ","
                         << bot.assigned_frontier.c << ","
                         << std::setprecision(3) << gx_m << ","
                         << gy_m << ","
                         << bot.tour.waypoints.size() << ","
                         << bot.stuck_ticks << ","
                         << bot.persistence_ticks << ","
                         << std::setprecision(3) << robot_total_distance_m[i] << "\n";
    }

    // 3. Write active frontier clusters (sample every 4 ticks to keep log size compact)
    if (sim.total_ticks % 4 == 0) {
        for (size_t f = 0; f < sim.frontiers.size(); ++f) {
            const auto& fr = sim.frontiers[f];
            auto [fx, fy] = sim.cellToWorld(fr.representative);
            int assigned_bot = -1;
            for (size_t i = 0; i < N; ++i) {
                if (sim.robots[i].assigned_frontier == fr.representative) {
                    assigned_bot = static_cast<int>(i);
                    break;
                }
            }

            frontiers_csv << sim.total_ticks << ","
                          << std::fixed << std::setprecision(3) << sim_time_s << ","
                          << f << ","
                          << fr.representative.r << ","
                          << fr.representative.c << ","
                          << std::setprecision(3) << fx << ","
                          << fy << ","
                          << fr.cells.size() << ","
                          << assigned_bot << "\n";
        }
    }

    logged_steps++;
}

void MissionLogger::logBfsCompare(uint64_t job_id, uint64_t tick, size_t num_points,
                                 double htrip_query_us, double bfs_query_us, double astar_query_us) {
    if (!is_active) return;
    const double speedup = (htrip_query_us > 0.0) ? (bfs_query_us / htrip_query_us) : 0.0;
    const double astar_speedup = (htrip_query_us > 0.0) ? (astar_query_us / htrip_query_us) : 0.0;
    std::lock_guard<std::mutex> lock(bfs_compare_mutex);
    if (!bfs_compare_csv.is_open()) return;
    bfs_compare_csv << job_id << ","
                    << tick << ","
                    << num_points << ","
                    << std::fixed << std::setprecision(3)
                    << htrip_query_us << ","
                    << bfs_query_us << ","
                    << std::setprecision(6) << speedup << ","
                    << std::fixed << std::setprecision(3)
                    << astar_query_us << ","
                    << std::setprecision(6) << astar_speedup << "\n";
    bfs_compare_csv.flush();
}

void MissionLogger::endSession(const MultiRobotSimulator& sim) {
    if (!is_active) return;

    telemetry_csv.flush();
    trajectories_csv.flush();
    frontiers_csv.flush();
    {
        std::lock_guard<std::mutex> lock(bfs_compare_mutex);
        bfs_compare_csv.flush();
        bfs_compare_csv.close();
    }

    telemetry_csv.close();
    trajectories_csv.close();
    frontiers_csv.close();

    // Generate comprehensive academic summary JSON
    double sim_time_s = static_cast<double>(sim.total_ticks) * 0.05;
    auto now = std::chrono::steady_clock::now();
    double wall_time_s = std::chrono::duration<double>(now - session_start_time).count();

    double mean_q_us = 0.0;
    double median_q_us = 0.0;
    double p95_q_us = 0.0;
    if (!htrip_query_latencies_us.empty()) {
        double sum = std::accumulate(htrip_query_latencies_us.begin(), htrip_query_latencies_us.end(), 0.0);
        mean_q_us = sum / static_cast<double>(htrip_query_latencies_us.size());
        auto copy = htrip_query_latencies_us;
        median_q_us = calculateMedian(copy);
        p95_q_us = calculatePercentile(copy, 95.0);
    }

    std::vector<double> paired_htrip_us;
    std::vector<double> paired_bfs_us;
    std::vector<double> paired_speedup;
    paired_htrip_us.reserve(sim.bfs_compare_samples.size());
    paired_bfs_us.reserve(sim.bfs_compare_samples.size());
    paired_speedup.reserve(sim.bfs_compare_samples.size());
    for (const auto& sample : sim.bfs_compare_samples) {
        paired_htrip_us.push_back(sample.htrip_query_us);
        paired_bfs_us.push_back(sample.bfs_query_us);
        paired_speedup.push_back(sample.speedup);
    }

    double mean_bfs_us = 0.0;
    double median_bfs_us = 0.0;
    if (!paired_bfs_us.empty()) {
        double sum = std::accumulate(paired_bfs_us.begin(), paired_bfs_us.end(), 0.0);
        mean_bfs_us = sum / static_cast<double>(paired_bfs_us.size());
        auto copy = paired_bfs_us;
        median_bfs_us = calculateMedian(copy);
    }

    double mean_paired_htrip_us = 0.0;
    if (!paired_htrip_us.empty()) {
        double sum = std::accumulate(paired_htrip_us.begin(), paired_htrip_us.end(), 0.0);
        mean_paired_htrip_us = sum / static_cast<double>(paired_htrip_us.size());
    }

    double mean_speedup = 0.0;
    double max_speedup = 0.0;
    if (!paired_speedup.empty()) {
        double sum = std::accumulate(paired_speedup.begin(), paired_speedup.end(), 0.0);
        mean_speedup = sum / static_cast<double>(paired_speedup.size());
        max_speedup = *std::max_element(paired_speedup.begin(), paired_speedup.end());
    }

    double mean_upd_us = 0.0;
    double median_upd_us = 0.0;
    double p95_upd_us = 0.0;
    if (!htrip_update_latencies_us.empty()) {
        double sum = std::accumulate(htrip_update_latencies_us.begin(), htrip_update_latencies_us.end(), 0.0);
        mean_upd_us = sum / static_cast<double>(htrip_update_latencies_us.size());
        auto copy = htrip_update_latencies_us;
        median_upd_us = calculateMedian(copy);
        p95_upd_us = calculatePercentile(copy, 95.0);
    }

    double lag_p95 = 0.0;
    double lag_max = 0.0;
    double lag_busy_frac = 0.0;
    if (!replica_lag_ticks_samples.empty()) {
        auto copy = replica_lag_ticks_samples;
        lag_p95 = calculatePercentile(copy, 95.0);
        lag_max = *std::max_element(replica_lag_ticks_samples.begin(), replica_lag_ticks_samples.end());
        size_t busy = 0;
        for (double v : replica_lag_ticks_samples) {
            if (v > 0.0) ++busy;
        }
        lag_busy_frac = static_cast<double>(busy) / static_cast<double>(replica_lag_ticks_samples.size());
    }

    double total_fleet_dist = std::accumulate(robot_total_distance_m.begin(), robot_total_distance_m.end(), 0.0);

    std::ofstream summary_json(std::filesystem::path(session_dir) / "mission_summary.json");
    if (summary_json.is_open()) {
        summary_json << "{\n"
                     << "  \"mission_metadata\": {\n"
                     << "    \"session_dir\": \"" << session_dir << "\",\n"
                     << "    \"total_ticks\": " << sim.total_ticks << ",\n"
                     << "    \"total_sim_time_s\": " << sim_time_s << ",\n"
                     << "    \"total_wall_time_s\": " << wall_time_s << ",\n"
                     << "    \"logged_steps\": " << logged_steps << "\n"
                     << "  },\n"
                     << "  \"environment\": {\n"
                     << "    \"grid_rows\": " << sim.observed.rows << ",\n"
                     << "    \"grid_cols\": " << sim.observed.cols << ",\n"
                     << "    \"scale_factor\": " << sim.scale_factor << ",\n"
                     << "    \"cell_size_m\": " << sim.cell_size_m << ",\n"
                     << "    \"total_free_cells\": " << sim.total_free_cells_gt << ",\n"
                     << "    \"explored_free_cells\": " << sim.explored_free_cells << ",\n"
                     << "    \"final_coverage_pct\": " << sim.explorationFraction() * 100.0 << "\n"
                     << "  },\n"
                     << "  \"htrip_hierarchy\": {\n"
                     << "    \"tile_size_b\": " << sim.oracle.config.tile_size_b << ",\n"
                     << "    \"group_w\": " << sim.oracle.config.group_w << ",\n"
                     << "    \"tree_depth\": " << sim.oracle.depth() + 1 << ",\n"
                     << "    \"total_leaves\": " << sim.oracle.leaves.size() << "\n"
                     << "  },\n"
                     << "  \"repeatability\": {\n"
                     << "    \"robot_seed\": " << sim.getSeed() << ",\n"
                     << "    \"max_frontiers_eval\": " << sim.max_frontiers_eval << ",\n"
                     << "    \"replan_period_ticks\": " << sim.replan_period_ticks << ",\n"
                     << "    \"persistence_window_ticks\": " << sim.persistence_window_ticks << ",\n"
                     << "    \"enable_async_updates\": " << (sim.enable_async_updates ? "true" : "false") << ",\n"
                     << "    \"enable_auto_respawn\": " << (sim.enable_auto_respawn ? "true" : "false") << ",\n"
                     << "    \"auto_respawn_ticks\": " << sim.auto_respawn_ticks << ",\n"
                     << "    \"enable_live_frontier_projection\": " << (sim.enable_live_frontier_projection ? "true" : "false") << ",\n"
                     << "    \"enable_bfs_baseline\": " << (sim.enable_bfs_baseline ? "true" : "false") << ",\n"
                     << "    \"bfs_compare_interval\": " << sim.bfs_compare_interval << ",\n"
                     << "    \"lidar_range_m\": " << sim.robot_config.lidar.range_m << ",\n"
                     << "    \"lidar_span_deg\": " << sim.robot_config.lidar.span_deg << ",\n"
                     << "    \"lidar_num_rays\": " << sim.robot_config.lidar.num_rays << ",\n"
                     << "    \"v_max\": " << sim.robot_config.v_max << ",\n"
                     << "    \"omega_max\": " << sim.robot_config.omega_max << ",\n"
                     << "    \"radius_m\": " << sim.robot_config.radius_m << "\n"
                     << "  },\n"
                     << "  \"fleet\": {\n"
                     << "    \"robot_count\": " << sim.robots.size() << ",\n"
                     << "    \"total_respawns\": " << sim.total_robot_respawns << ",\n"
                     << "    \"total_fleet_distance_m\": " << total_fleet_dist << ",\n"
                     << "    \"per_robot_distance_m\": [";
        for (size_t i = 0; i < robot_total_distance_m.size(); ++i) {
            summary_json << robot_total_distance_m[i] << (i + 1 < robot_total_distance_m.size() ? ", " : "");
        }
        summary_json << "]\n"
                     << "  },\n"
                     << "  \"performance_metrics\": {\n"
                     << "    \"htrip_query_mean_us\": " << mean_q_us << ",\n"
                     << "    \"htrip_query_median_us\": " << median_q_us << ",\n"
                     << "    \"htrip_query_p95_us\": " << p95_q_us << ",\n"
                     << "    \"bfs_compare_count\": " << sim.bfs_compare_samples.size() << ",\n"
                     << "    \"paired_htrip_query_mean_us\": " << mean_paired_htrip_us << ",\n"
                     << "    \"bfs_query_mean_us\": " << mean_bfs_us << ",\n"
                     << "    \"bfs_query_median_us\": " << median_bfs_us << ",\n"
                     << "    \"mean_speedup_vs_bfs\": " << mean_speedup << ",\n"
                     << "    \"max_speedup_vs_bfs\": " << max_speedup << ",\n"
                     << "    \"htrip_update_mean_us\": " << mean_upd_us << ",\n"
                     << "    \"htrip_update_median_us\": " << median_upd_us << ",\n"
                     << "    \"htrip_update_p95_us\": " << p95_upd_us << ",\n"
                     << "    \"replica_lag_p95_ticks\": " << lag_p95 << ",\n"
                     << "    \"replica_lag_max_ticks\": " << lag_max << ",\n"
                     << "    \"replica_lag_busy_frac\": " << lag_busy_frac << "\n"
                     << "  }\n"
                     << "}\n";
    }

    is_active = false;
}

} // namespace htrip
