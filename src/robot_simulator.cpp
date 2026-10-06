#include "htrip/robot_simulator.hpp"
#include "htrip/bfs_oracle.hpp"
#include "htrip/astar_oracle.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>
#include <random>

namespace htrip {

namespace {

inline double normalizeAngle(double a) noexcept {
    while (a > M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

inline bool hasLineOfSight(const GridMap& grid, Point a, Point b) noexcept {
    if (a == b) return true;
    int r0 = a.r, c0 = a.c;
    int r1 = b.r, c1 = b.c;
    int dr = std::abs(r1 - r0), dc = std::abs(c1 - c0);
    int sr = (r0 < r1) ? 1 : -1, sc = (c0 < c1) ? 1 : -1;
    int err = dc - dr;
    while (true) {
        if (!grid.inBounds(r0, c0) || !grid.isPassable(r0, c0)) return false;
        if (r0 == r1 && c0 == c1) break;
        int e2 = 2 * err;
        if (e2 > -dr) { err -= dr; c0 += sc; }
        if (e2 < dc)  { err += dc; r0 += sr; }
    }
    return true;
}

const uint32_t ROBOT_PALETTE[] = {
    0xFFEE8833, // Vibrant Cyan/Blue
    0xFF33CC55, // Emerald Green
    0xFF3388EE, // Coral Orange
    0xFFDD44BB, // Orchid Purple
    0xFF22DDBB, // Mint Green
    0xFF4455EE, // Crimson Red
    0xFFEECC22, // Sky Blue
    0xFF8844DD  // Rose Magenta
};

} // anonymous namespace

void MultiRobotSimulator::init(const GridMap& raw_gt,
                               double raw_cell_size_m,
                               int scale,
                               int num_robots,
                               const RobotConfig& cfg,
                               int b,
                               int g,
                               uint32_t seed,
                               bool randomize_poses,
                               bool enable_fog) {
    // 0. Gracefully stop and join any active background worker threads
    stopBackgroundThreads();
    enable_fog_of_war = enable_fog;
    scale_factor = std::clamp(scale, 1, 4);
    cell_size_m = raw_cell_size_m / static_cast<double>(scale_factor);
    robot_config = cfg;
    total_ticks = 0;

    // Seed initialization: user-provided seed or clock time for reproducibility
    if (seed != 0) {
        robot_seed = seed;
    } else if (robot_seed == 0) {
        robot_seed = static_cast<uint32_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }

    int R = raw_gt.rows * scale_factor;
    int C = raw_gt.cols * scale_factor;
    world_w_m = static_cast<double>(C) * cell_size_m;
    world_h_m = static_cast<double>(R) * cell_size_m;

    // 1. Build the upscaled ground truth and observed grids:
    // - ground_truth: Complete obstacle geometry loaded from MovingAI benchmark.
    // - observed: Perception layer revealed in real time as LiDAR rays sweep (initially Unknown).
    // - indexed_map: Strictly tracks the exact grid state reflected in the H-TRIP hierarchy.
    ground_truth = GridMap(R, C, CellState::Free);
    if (!enable_fog_of_war) {
        observed = GridMap(R, C, CellState::Free);
        indexed_map = GridMap(R, C, CellState::Free);
    } else {
        observed = GridMap(R, C, CellState::Unknown);
        indexed_map = GridMap(R, C, CellState::Unknown);
    }

    total_free_cells_gt = 0;
    for (int r0 = 0; r0 < raw_gt.rows; ++r0) {
        for (int c0 = 0; c0 < raw_gt.cols; ++c0) {
            CellState cs = raw_gt.getObserved(r0, c0);
            if (cs == CellState::Unknown && !raw_gt.ground_truth.empty()) {
                cs = raw_gt.ground_truth[static_cast<size_t>(raw_gt.index(r0, c0))];
            }
            for (int dr = 0; dr < scale_factor; ++dr) {
                for (int dc = 0; dc < scale_factor; ++dc) {
                    int r = r0 * scale_factor + dr;
                    int c = c0 * scale_factor + dc;
                    ground_truth.setObserved(r, c, cs);
                    if (!enable_fog_of_war) {
                        observed.setObserved(r, c, cs);
                        indexed_map.setObserved(r, c, cs);
                    }
                    if (cs == CellState::Free) {
                        total_free_cells_gt++;
                    }
                }
            }
        }
    }

    explored_free_cells = (!enable_fog_of_war) ? total_free_cells_gt : 0;
    total_ticks = 0;
    last_replan_tick = 0;
    frontiers.clear();
    viewpoints.clear();
    all_pairs_matrix.clear();
    newly_opened_cells.clear();
    updated_leaf_indices.clear();
    last_htrip_update_us = 0.0;
    last_htrip_query_us = 0.0;
    last_mtsp_plan_us = 0.0;
    last_bfs_query_us = 0.0;
    last_astar_query_us = 0.0;
    last_paired_htrip_query_us = 0.0;
    last_paired_speedup = 0.0;
    last_paired_astar_speedup = 0.0;
    bfs_compare_count = 0;
    bfs_compare_samples.clear();
    bfs_inflight_job = {};
    bfs_pending_job = {};
    bfs_await_htrip_job = {};
    last_query_pts.clear();

    // 2. Spawn N robots in free, obstacle-cleared locations (reveals surroundings into observed)
    robots.clear();
    setRobotCount(num_robots, randomize_poses, robot_seed);

    // 3. Initialize H-TRIP Oracle directly on observed map
    HierarchyConfig hcfg;
    hcfg.tile_size_b = b;
    hcfg.group_w = g;
    hcfg.group_h = g;
    oracle = MinPlusOracle(hcfg);
    oracle.build(observed);
    indexed_map = observed;
    grid_tile_size_b = oracle.config.tile_size_b;
    grid_num_leaf_tiles_x = oracle.num_leaf_tiles_x;
    grid_num_leaf_tiles_y = oracle.num_leaf_tiles_y;
    grid_num_leaves = oracle.leaves.size();
    double_buffer_dirty = false;
    {
        std::lock_guard<std::mutex> q_lock(query_oracle_mutex);
        query_oracle = oracle;
        query_map = indexed_map;
        active_query_oracle_ptr.store(&query_oracle);
        active_query_map_ptr.store(&query_map);
        active_update_oracle_ptr.store(&oracle);
        active_update_map_ptr.store(&indexed_map);
    }
    {
        std::lock_guard<std::mutex> lock(commit_mutex);
        completed_indexed_cells.clear();
    }

    newly_opened_cells.clear();
    updated_leaf_indices.clear();
    {
        std::lock_guard<std::mutex> lock(decay_mutex);
        active_tile_decay.assign(oracle.leaves.size(), 0);
        active_index_decay.assign(oracle.leaves.size(), 0);
    }
    total_map_updates = 0;
    last_update_tick = 0;
    replica_lag_ticks = 0;
    pending_openings = 0;
    map_update_freq_hz = 0.0;
    sim_step_freq_hz = 0.0;
    last_step_time = {};
    staged_all_pairs_matrix.clear();
    staged_fleet_tours.clear();
    staged_viewpoints.clear();
    new_tours_ready.store(false);
    htrip_query_freq_hz = 0.0;
    total_htrip_queries = 0;
    last_query_tick = 0;

    if (!enable_fog_of_war) {
        initPatrolTargets(static_cast<size_t>(max_frontiers_eval > 0 ? max_frontiers_eval : 25));
    }
}

MultiRobotSimulator::~MultiRobotSimulator() {
    if (mission_logger.isLogging()) {
        stopMissionLogging();
    }
    stopBackgroundThreads();
}

void MultiRobotSimulator::stopBackgroundThreads() noexcept {
    while (is_updating_htrip.load() || is_running_bfs.load() || is_querying_pairwise.load()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    try {
        if (update_worker_thread.joinable()) {
            update_worker_thread.join();
        }
        if (query_worker_thread.joinable()) {
            query_worker_thread.join();
        }
        if (bfs_worker_thread.joinable()) {
            bfs_worker_thread.join();
        }
    } catch (...) {
        // Guard against std::system_error in noexcept
    }
    std::lock_guard<std::mutex> lock(commit_mutex);
    for (const auto& pt : completed_indexed_cells) {
        indexed_map.setObserved(pt, CellState::Free);
    }
    completed_indexed_cells.clear();
}

void MultiRobotSimulator::reset(bool randomize_poses, uint32_t new_seed) {
    stopBackgroundThreads();
    std::lock_guard<std::mutex> lock(oracle_mutex);
    if (new_seed != 0) {
        robot_seed = new_seed;
    }
    async_pending_cells.clear();
    // Reset observed and indexed grid
    if (!enable_fog_of_war) {
        observed = ground_truth;
        indexed_map = ground_truth;
        explored_free_cells = total_free_cells_gt;
    } else {
        for (int r = 0; r < observed.rows; ++r) {
            for (int c = 0; c < observed.cols; ++c) {
                observed.setObserved(r, c, CellState::Unknown);
                indexed_map.setObserved(r, c, CellState::Unknown);
            }
        }
        explored_free_cells = 0;
    }
    total_ticks = 0;
    last_replan_tick = 0;
    newly_opened_cells.clear();
    updated_leaf_indices.clear();
    frontiers.clear();
    viewpoints.clear();
    all_pairs_matrix.clear();
    last_query_pts.clear();

    // Respawn robots and reveal initial surroundings
    int n = static_cast<int>(robots.size());
    robots.clear();
    setRobotCount(n, randomize_poses, robot_seed);

    // Rebuild Oracle directly on observed map
    oracle.build(observed);
    indexed_map = observed;
    grid_tile_size_b = oracle.config.tile_size_b;
    grid_num_leaf_tiles_x = oracle.num_leaf_tiles_x;
    grid_num_leaf_tiles_y = oracle.num_leaf_tiles_y;
    grid_num_leaves = oracle.leaves.size();
    double_buffer_dirty = false;
    {
        std::lock_guard<std::mutex> q_lock(query_oracle_mutex);
        query_oracle = oracle;
        query_map = indexed_map;
        active_query_oracle_ptr.store(&query_oracle);
        active_query_map_ptr.store(&query_map);
        active_update_oracle_ptr.store(&oracle);
        active_update_map_ptr.store(&indexed_map);
    }
    {
        std::lock_guard<std::mutex> lock(commit_mutex);
        completed_indexed_cells.clear();
    }

    newly_opened_cells.clear();
    updated_leaf_indices.clear();
    {
        std::lock_guard<std::mutex> lock(decay_mutex);
        active_tile_decay.assign(oracle.leaves.size(), 0);
        active_index_decay.assign(oracle.leaves.size(), 0);
    }
    total_map_updates = 0;
    last_update_tick = 0;
    replica_lag_ticks = 0;
    pending_openings = 0;
    map_update_freq_hz = 0.0;
    sim_step_freq_hz = 0.0;
    last_step_time = {};
    last_htrip_update_us = 0.0;
    last_htrip_query_us = 0.0;
    last_mtsp_plan_us = 0.0;
    last_bfs_query_us = 0.0;
    last_astar_query_us = 0.0;
    last_paired_htrip_query_us = 0.0;
    last_paired_speedup = 0.0;
    last_paired_astar_speedup = 0.0;
    bfs_compare_count = 0;
    bfs_compare_samples.clear();
    bfs_inflight_job = {};
    bfs_pending_job = {};
    bfs_await_htrip_job = {};
    staged_all_pairs_matrix.clear();
    staged_fleet_tours.clear();
    staged_viewpoints.clear();
    new_tours_ready.store(false);
    htrip_query_freq_hz = 0.0;
    total_htrip_queries = 0;
    last_query_tick = 0;

    if (!enable_fog_of_war) {
        initPatrolTargets(static_cast<size_t>(max_frontiers_eval > 0 ? max_frontiers_eval : 25));
    }
}

void MultiRobotSimulator::initPatrolTargets(size_t count) {
    patrol_targets.clear();
    targets_visited_count = 0;
    if (count == 0) return;

    std::vector<Point> candidates;
    for (int r = 0; r < ground_truth.rows; ++r) {
        for (int c = 0; c < ground_truth.cols; ++c) {
            if (ground_truth.isPassable(r, c) && ground_truth.getClearance(Point{r, c}) >= 1) {
                candidates.push_back({r, c});
            }
        }
    }
    if (candidates.empty()) return;

    std::mt19937_64 rng(robot_seed ^ 0x9e3779b97f4a7c15ULL);

    int grid_div = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(count))));
    grid_div = std::max(1, grid_div);
    int sector_h = (ground_truth.rows + grid_div - 1) / grid_div;
    int sector_w = (ground_truth.cols + grid_div - 1) / grid_div;

    std::vector<std::vector<Point>> sectors(static_cast<size_t>(grid_div * grid_div));
    for (const auto& pt : candidates) {
        int sr = std::clamp(pt.r / sector_h, 0, grid_div - 1);
        int sc = std::clamp(pt.c / sector_w, 0, grid_div - 1);
        sectors[static_cast<size_t>(sr * grid_div + sc)].push_back(pt);
    }

    for (auto& sec : sectors) {
        if (!sec.empty()) {
            std::shuffle(sec.begin(), sec.end(), rng);
        }
    }

    size_t round = 0;
    while (patrol_targets.size() < count) {
        bool added_any = false;
        for (auto& sec : sectors) {
            if (round < sec.size()) {
                patrol_targets.push_back(sec[round]);
                added_any = true;
                if (patrol_targets.size() == count) break;
            }
        }
        if (!added_any) break;
        round++;
    }

    frontiers.clear();
    for (size_t i = 0; i < patrol_targets.size(); ++i) {
        FrontierCluster cl;
        cl.id = static_cast<int>(i + 1);
        cl.representative = patrol_targets[i];
        cl.cells = {patrol_targets[i]};
        frontiers.push_back(std::move(cl));
    }
    viewpoints = patrol_targets;
}

void MultiRobotSimulator::replenishPatrolTargets(size_t target_count) {
    if (patrol_targets.size() >= target_count) return;

    std::vector<Point> candidates;
    for (int r = 0; r < ground_truth.rows; ++r) {
        for (int c = 0; c < ground_truth.cols; ++c) {
            if (ground_truth.isPassable(r, c) && ground_truth.getClearance(Point{r, c}) >= 1) {
                candidates.push_back({r, c});
            }
        }
    }
    if (candidates.empty()) return;

    std::mt19937_64 rng(robot_seed + static_cast<uint32_t>(total_ticks * 31 + targets_visited_count * 17));
    std::shuffle(candidates.begin(), candidates.end(), rng);

    for (const auto& cand : candidates) {
        if (patrol_targets.size() >= target_count) break;
        bool too_close = false;
        for (const auto& existing : patrol_targets) {
            if (std::hypot(cand.r - existing.r, cand.c - existing.c) < 16.0) {
                too_close = true;
                break;
            }
        }
        if (!too_close) {
            patrol_targets.push_back(cand);
        }
    }
}

void MultiRobotSimulator::setRobotCount(int num_robots, bool randomize_poses, uint32_t seed) {
    num_robots = std::max(1, num_robots);
    if (seed != 0) {
        robot_seed = seed;
    }
    int current_n = static_cast<int>(robots.size());
    if (current_n == num_robots && !randomize_poses) return;

    if (num_robots < current_n && !randomize_poses) {
        robots.resize(static_cast<size_t>(num_robots));
        return;
    }

    robots.resize(static_cast<size_t>(num_robots));

    // Find candidate free cells for robot spawning
    double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
    std::vector<Point> spawn_candidates;
    int step = std::max(1, ground_truth.rows / 32);
    for (int r = 1; r < ground_truth.rows - 1; r += step) {
        for (int c = 1; c < ground_truth.cols - 1; c += step) {
            if (ground_truth.isPassable(r, c)) {
                auto [wx, wy] = cellToWorld({r, c});
                if (checkRobotCollision(wx, wy, eff_r)) {
                    spawn_candidates.push_back({r, c});
                }
            }
        }
    }
    if (spawn_candidates.empty()) {
        for (int r = 1; r < ground_truth.rows - 1; ++r) {
            for (int c = 1; c < ground_truth.cols - 1; ++c) {
                if (ground_truth.isPassable(r, c)) {
                    spawn_candidates.push_back({r, c});
                }
            }
        }
    }

    if (randomize_poses && !spawn_candidates.empty()) {
        uint32_t effective_seed = (seed != 0) ? seed : robot_seed;
        if (effective_seed == 0) {
            effective_seed = static_cast<uint32_t>(
                std::chrono::high_resolution_clock::now().time_since_epoch().count());
            robot_seed = effective_seed;
        }
        std::mt19937 rng(effective_seed);
        std::shuffle(spawn_candidates.begin(), spawn_candidates.end(), rng);
        std::uniform_real_distribution<double> angle_dist(-M_PI, M_PI);

        // Filter for well-spaced candidates so robots don't spawn overlapping
        std::vector<Point> chosen_spawns;
        for (const auto& cand : spawn_candidates) {
            if (chosen_spawns.size() == static_cast<size_t>(num_robots)) break;
            auto [cx, cy] = cellToWorld(cand);
            bool too_close = false;
            for (const auto& prev : chosen_spawns) {
                auto [px, py] = cellToWorld(prev);
                double dist = std::hypot(cx - px, cy - py);
                if (dist < 2.0 * robot_config.radius_m) {
                    too_close = true;
                    break;
                }
            }
            if (!too_close) {
                chosen_spawns.push_back(cand);
            }
        }
        size_t cand_idx = 0;
        while (chosen_spawns.size() < static_cast<size_t>(num_robots) && !spawn_candidates.empty()) {
            chosen_spawns.push_back(spawn_candidates[cand_idx % spawn_candidates.size()]);
            cand_idx++;
        }

        for (int i = 0; i < num_robots; ++i) {
            RobotState& bot = robots[static_cast<size_t>(i)];
            bot.id = i;
            bot.color = ROBOT_PALETTE[static_cast<size_t>(i % 8)];
            bot.v = 0.0;
            bot.omega = 0.0;
            bot.stuck_ticks = 0;
            bot.escape_ticks = 0;
            bot.escape_turn_dir = (i % 2 == 0) ? 1.0 : -1.0;
            bot.trail.clear();
            bot.assigned_frontier = {-1, -1};
            bot.assigned_frontier_id = -1;
            bot.has_goal = false;
            bot.next_waypoint = {-1, -1};
            bot.backtrack_ticks = 0;
            bot.backtrack_target_x = 0.0;
            bot.backtrack_target_y = 0.0;
            bot.waypoint_stuck_count = 0;
            bot.tour.clear();
            bot.planned_path.clear();
            bot.path_index = 0;

            Point sp = chosen_spawns[static_cast<size_t>(i % chosen_spawns.size())];
            auto [wx, wy] = cellToWorld(sp);
            bot.x = wx;
            bot.y = wy;
            bot.stuck_ref_x = wx;
            bot.stuck_ref_y = wy;
            bot.stuck_translation_ticks = 0;
            bot.total_respawns = 0;
            bot.theta = angle_dist(rng);
            bot.wander_angle_offset = bot.theta;

            castLidarRays(bot);
        }
        return;
    }

    int start_idx = (current_n == 0) ? 0 : current_n;
    for (int i = start_idx; i < num_robots; ++i) {
        RobotState& bot = robots[static_cast<size_t>(i)];
        bot.id = i;
        bot.color = ROBOT_PALETTE[static_cast<size_t>(i % 8)];
        bot.v = 0.0;
        bot.omega = 0.0;
        bot.stuck_ticks = 0;
        bot.escape_ticks = 0;
        bot.escape_turn_dir = (i % 2 == 0) ? 1.0 : -1.0;
        bot.wander_angle_offset = (2.0 * M_PI * i) / static_cast<double>(num_robots);
        bot.trail.clear();
        bot.assigned_frontier = {-1, -1};
        bot.assigned_frontier_id = -1;
        bot.has_goal = false;
        bot.next_waypoint = {-1, -1};
        bot.tour.clear();
        bot.planned_path.clear();
        bot.path_index = 0;

        if (!spawn_candidates.empty()) {
            size_t stride = std::max(size_t{1}, spawn_candidates.size() / static_cast<size_t>(num_robots));
            size_t idx = (static_cast<size_t>(i) * stride) % spawn_candidates.size();
            Point sp = spawn_candidates[idx];
            auto [wx, wy] = cellToWorld(sp);
            bot.x = wx;
            bot.y = wy;
        } else {
            bot.x = world_w_m / 2.0;
            bot.y = world_h_m / 2.0;
        }
        bot.stuck_ref_x = bot.x;
        bot.stuck_ref_y = bot.y;
        bot.stuck_translation_ticks = 0;
        bot.total_respawns = 0;
        bot.theta = bot.wander_angle_offset;

        // Immediately cast initial scan to reveal surroundings
        castLidarRays(bot);
    }
}

void MultiRobotSimulator::updateLidarConfig(const LidarConfig& lcfg) {
    robot_config.lidar = lcfg;
    for (auto& bot : robots) {
        castLidarRays(bot);
    }
}

bool MultiRobotSimulator::respawnRobot(int robot_id) {
    if (robot_id < 0 || robot_id >= static_cast<int>(robots.size())) return false;
    RobotState& bot = robots[static_cast<size_t>(robot_id)];

    double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);

    // 1. Collect candidate passable cells in ground truth
    std::vector<Point> candidates;
    int step_size = std::max(1, ground_truth.rows / 32);
    for (int r = 1; r < ground_truth.rows - 1; r += step_size) {
        for (int c = 1; c < ground_truth.cols - 1; c += step_size) {
            if (ground_truth.isPassable(r, c)) {
                auto [wx, wy] = cellToWorld({r, c});
                if (checkRobotCollision(wx, wy, eff_r)) {
                    candidates.push_back({r, c});
                }
            }
        }
    }
    if (candidates.empty()) {
        for (int r = 1; r < ground_truth.rows - 1; ++r) {
            for (int c = 1; c < ground_truth.cols - 1; ++c) {
                if (ground_truth.isPassable(r, c)) {
                    auto [wx, wy] = cellToWorld({r, c});
                    if (checkRobotCollision(wx, wy, eff_r)) {
                        candidates.push_back({r, c});
                    }
                }
            }
        }
    }
    if (candidates.empty()) return false;

    // 2. Score candidates: pick a safe, open location away from the stuck position and other robots
    std::mt19937 rng(static_cast<uint32_t>(robot_seed + total_ticks * 31 + robot_id * 1009));
    std::shuffle(candidates.begin(), candidates.end(), rng);

    Point best_cand = candidates[0];
    double best_score = -1e9;

    for (const auto& cand : candidates) {
        auto [cx, cy] = cellToWorld(cand);
        double dist_from_stuck = std::hypot(cx - bot.x, cy - bot.y);
        if (dist_from_stuck < std::max(1.5, cell_size_m * 10.0) && candidates.size() > 5) {
            continue; // avoid respawning at the same stuck trap
        }

        double min_bot_dist = 1e9;
        for (const auto& other : robots) {
            if (other.id == robot_id) continue;
            double d = std::hypot(cx - other.x, cy - other.y);
            if (d < min_bot_dist) min_bot_dist = d;
        }

        if (min_bot_dist < 2.0 * robot_config.radius_m) {
            continue; // avoid overlapping with other robots
        }

        double min_frontier_d = 1e9;
        for (const auto& vp : viewpoints) {
            auto [vx, vy] = cellToWorld(vp);
            double d = std::hypot(cx - vx, cy - vy);
            if (d < min_frontier_d) min_frontier_d = d;
        }

        double score = std::min(min_bot_dist, 5.0) * 2.0 + std::min(dist_from_stuck, 10.0);
        if (min_frontier_d < 1e8) {
            score -= std::min(min_frontier_d, 20.0) * 0.4;
        }

        if (score > best_score) {
            best_score = score;
            best_cand = cand;
        }
    }

    auto [spawn_x, spawn_y] = cellToWorld(best_cand);
    bot.x = spawn_x;
    bot.y = spawn_y;
    bot.stuck_ref_x = spawn_x;
    bot.stuck_ref_y = spawn_y;
    bot.stuck_translation_ticks = 0;
    bot.stuck_ticks = 0;
    bot.v = 0.0;
    bot.omega = 0.0;

    std::uniform_real_distribution<double> angle_dist(-M_PI, M_PI);
    bot.theta = angle_dist(rng);
    bot.wander_angle_offset = bot.theta;

    bot.recovery_phase = RobotState::RecoveryPhase::None;
    bot.recovery_ticks = 0;
    bot.waypoint_stuck_count = 0;
    bot.has_goal = false;
    bot.assigned_frontier = {-1, -1};
    bot.assigned_frontier_id = -1;
    bot.planned_path.clear();
    bot.path_index = 0;
    bot.tour.clear();
    bot.trail.clear();
    bot.trail.push_back({static_cast<float>(bot.x), static_cast<float>(bot.y)});

    bot.total_respawns++;
    total_robot_respawns++;

    // Immediately reveal surroundings at new spawn point
    castLidarRays(bot);
    return true;
}

int MultiRobotSimulator::checkAndRespawnStuckRobots() {
    if (!enable_auto_respawn || auto_respawn_ticks <= 0) return 0;
    int count = 0;
    for (auto& bot : robots) {
        if (bot.stuck_translation_ticks >= auto_respawn_ticks) {
            if (respawnRobot(bot.id)) {
                count++;
            }
        }
    }
    return count;
}

bool MultiRobotSimulator::checkRobotCollision(double x, double y, double r) const noexcept {
    double eff_r = std::min(r, cell_size_m * 0.38);
    int min_c = static_cast<int>(std::floor((x - eff_r) / cell_size_m));
    int max_c = static_cast<int>(std::floor((x + eff_r) / cell_size_m));
    int min_r = static_cast<int>(std::floor((y - eff_r) / cell_size_m));
    int max_r = static_cast<int>(std::floor((y + eff_r) / cell_size_m));

    for (int cr = min_r; cr <= max_r; ++cr) {
        for (int cc = min_c; cc <= max_c; ++cc) {
            if (!ground_truth.inBounds(cr, cc)) return false; // Map boundary collision
            if (ground_truth.getObserved(cr, cc) == CellState::Obstacle) {
                // Circle-box overlap
                double cell_min_x = cc * cell_size_m;
                double cell_max_x = (cc + 1) * cell_size_m;
                double cell_min_y = cr * cell_size_m;
                double cell_max_y = (cr + 1) * cell_size_m;

                double closest_x = std::clamp(x, cell_min_x, cell_max_x);
                double closest_y = std::clamp(y, cell_min_y, cell_max_y);
                double dx = x - closest_x;
                double dy = y - closest_y;
                if ((dx * dx + dy * dy) < (eff_r * eff_r)) {
                    return false; // Collision detected
                }
            }
        }
    }
    return true; // Collision free
}

bool MultiRobotSimulator::isSegmentClear(double x1, double y1, double x2, double y2, double r) const noexcept {
    double dist = std::hypot(x2 - x1, y2 - y1);
    double eff_r = std::min(r, cell_size_m * 0.38);

    // Sample along the segment with fine spacing to guarantee no corner clipping
    double step_m = std::min(cell_size_m * 0.35, 0.035);
    int num_steps = static_cast<int>(std::ceil(dist / step_m));
    num_steps = std::max(1, num_steps);
    double dx = (x2 - x1) / num_steps;
    double dy = (y2 - y1) / num_steps;

    for (int i = 0; i <= num_steps; ++i) {
        double px = x1 + i * dx;
        double py = y1 + i * dy;

        int min_c = static_cast<int>(std::floor((px - eff_r) / cell_size_m));
        int max_c = static_cast<int>(std::floor((px + eff_r) / cell_size_m));
        int min_r = static_cast<int>(std::floor((py - eff_r) / cell_size_m));
        int max_r = static_cast<int>(std::floor((py + eff_r) / cell_size_m));

        for (int cr = min_r; cr <= max_r; ++cr) {
            for (int cc = min_c; cc <= max_c; ++cc) {
                if (!observed.inBounds(cr, cc)) return false;
                // In observed map, both Obstacle and Unknown (fog-of-war) block line of sight
                if (observed.getObserved(cr, cc) != CellState::Free) {
                    double cell_min_x = cc * cell_size_m;
                    double cell_max_x = (cc + 1) * cell_size_m;
                    double cell_min_y = cr * cell_size_m;
                    double cell_max_y = (cr + 1) * cell_size_m;

                    double closest_x = std::clamp(px, cell_min_x, cell_max_x);
                    double closest_y = std::clamp(py, cell_min_y, cell_max_y);
                    double ddx = px - closest_x;
                    double ddy = py - closest_y;
                    if ((ddx * ddx + ddy * ddy) < (eff_r * eff_r)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

/**
 * @brief Simulates planar LiDAR range sensing via discrete polar raymarching.
 *
 * For each beam uniformly spaced across the sensor FOV:
 * 1. Computes world beam angle ray_angle = normalizeAngle(robot.theta + rel_ang).
 * 2. Marches outward from robot pose along (cos(ray_angle), sin(ray_angle)) at step size 0.5 * cell_size_m.
 * 3. Checks the ground_truth map at the continuous endpoint:
 *    - If an obstacle cell is hit: terminates ray, records obstacle hit, and marks observed cell as Obstacle.
 *    - If free: marks cell as Free (if previously Unknown), appends to newly_opened_cells, and increments explored count.
 * 4. Stores polar hit payload in robot.scan for visualization and local obstacle repulsion.
 */
void MultiRobotSimulator::castLidarRays(RobotState& robot) {
    const auto& lcfg = robot_config.lidar;
    robot.scan.clear();
    robot.scan.reserve(static_cast<size_t>(lcfg.num_rays));

    double span_rad = lcfg.spanRad();
    double half_span = span_rad * 0.5;
    double start_angle = (lcfg.span_deg >= 360.0) ? -M_PI : -half_span;
    double delta_angle = (lcfg.span_deg >= 360.0)
        ? (2.0 * M_PI / static_cast<double>(lcfg.num_rays))
        : (span_rad / static_cast<double>(std::max(1, lcfg.num_rays - 1)));

    // Sub-cell step size guarantees that no thin diagonal obstacle walls are missed
    double step_size_m = cell_size_m * 0.5;
    int max_steps = static_cast<int>(std::ceil(lcfg.range_m / step_size_m));

    for (int i = 0; i < lcfg.num_rays; ++i) {
        double rel_ang = start_angle + static_cast<double>(i) * delta_angle;
        double ray_angle = normalizeAngle(robot.theta + rel_ang);
        double cos_a = std::cos(ray_angle);
        double sin_a = std::sin(ray_angle);

        double hit_dist = lcfg.range_m;
        bool hit_obstacle = false;
        double hx = robot.x + cos_a * lcfg.range_m;
        double hy = robot.y + sin_a * lcfg.range_m;

        for (int s = 1; s <= max_steps; ++s) {
            double cur_dist = s * step_size_m;
            if (cur_dist > lcfg.range_m) break;

            double cx = robot.x + cos_a * cur_dist;
            double cy = robot.y + sin_a * cur_dist;
            Point cp = worldToCell(cx, cy);

            CellState gt_state = ground_truth.getObserved(cp);
            if (gt_state == CellState::Obstacle) {
                hit_dist = cur_dist;
                hit_obstacle = true;
                hx = cx;
                hy = cy;

                if (observed.getObserved(cp) != CellState::Obstacle) {
                    observed.setObserved(cp, CellState::Obstacle);
                    indexed_map.setObserved(cp, CellState::Obstacle);
                }
                break;
            }

            // Cell is Free along the clear ray
            if (observed.getObserved(cp) == CellState::Unknown) {
                observed.setObserved(cp, CellState::Free);
                newly_opened_cells.push_back(cp);
                explored_free_cells++;
            }
        }

        robot.scan.push_back({ray_angle, hit_dist, hx, hy, hit_obstacle});
    }
}

void MultiRobotSimulator::updateDirtyTileIndices(const std::vector<Point>& cells) {
    int b = grid_tile_size_b;
    if (b <= 0 || grid_num_leaves == 0) return;
    if (dirty_leaf_stamp.size() != grid_num_leaves) {
        dirty_leaf_stamp.assign(grid_num_leaves, 0);
        dirty_stamp_epoch = 1;
    } else {
        dirty_stamp_epoch++;
        if (dirty_stamp_epoch == 0) {
            dirty_leaf_stamp.assign(grid_num_leaves, 0);
            dirty_stamp_epoch = 1;
        }
    }
    updated_leaf_indices.clear();
    std::lock_guard<std::mutex> decay_lock(decay_mutex);
    if (active_tile_decay.size() != grid_num_leaves) {
        active_tile_decay.assign(grid_num_leaves, 0);
    }
    if (active_index_decay.size() != grid_num_leaves) {
        active_index_decay.assign(grid_num_leaves, 0);
    }
    for (Point p : cells) {
        int tx = p.c / b;
        int ty = p.r / b;
        if (tx >= 0 && tx < grid_num_leaf_tiles_x && ty >= 0 && ty < grid_num_leaf_tiles_y) {
            int l_idx = ty * grid_num_leaf_tiles_x + tx;
            if (dirty_leaf_stamp[static_cast<size_t>(l_idx)] != dirty_stamp_epoch) {
                dirty_leaf_stamp[static_cast<size_t>(l_idx)] = dirty_stamp_epoch;
                updated_leaf_indices.push_back(l_idx);
                active_tile_decay[static_cast<size_t>(l_idx)] = 15; // Glow for 15 ticks (~0.75s)
            }
        }
    }
}

/**
 * @brief Flushes newly revealed free cells into the H-TRIP hierarchical index.
 *
 * Implements a non-blocking double-buffering protocol:
 * - Synchronous mode (enable_async_updates == false):
 *   Locks oracle_mutex, executes oracle.onBatchCellsOpened directly on the calling thread,
 *   marks double_buffer_dirty, and updates decay counters.
 * - Asynchronous mode (enable_async_updates == true):
 *   If the background update worker is busy, accumulates cells in async_pending_cells and returns immediately.
 *   When the worker is idle, gathers all accumulated cells, sets is_updating_htrip = true,
 *   and launches update_worker_thread to run dynamic rank-1 Conway updates in the background.
 *   Upon completion, the worker copies the updated state into query_oracle under query_oracle_mutex,
 *   guaranteeing that ongoing queries on query_oracle are never stalled.
 */
void MultiRobotSimulator::flushPendingUpdates() {
    if (newly_opened_cells.empty() && async_pending_cells.empty()) {
        return;
    }

    if (!enable_async_updates) {
        if (!async_pending_cells.empty()) {
            newly_opened_cells.insert(newly_opened_cells.end(), async_pending_cells.begin(), async_pending_cells.end());
            async_pending_cells.clear();
        }
        updateDirtyTileIndices(newly_opened_cells);
        last_revealed_cells_count = static_cast<int>(newly_opened_cells.size());

        total_map_updates++;
        uint64_t ticks_since = total_ticks - last_update_tick;
        last_update_tick = total_ticks;
        if (ticks_since > 0) {
            double inst_hz = 1.0 / (static_cast<double>(ticks_since) * 0.05);
            map_update_freq_hz = (map_update_freq_hz <= 0.0) ? inst_hz : (0.85 * map_update_freq_hz + 0.15 * inst_hz);
        }

        auto t0 = std::chrono::high_resolution_clock::now();
        {
            std::lock_guard<std::mutex> lock(oracle_mutex);
            oracle.onBatchCellsOpened(newly_opened_cells, observed);
            for (const auto& pt : newly_opened_cells) {
                indexed_map.setObserved(pt, CellState::Free);
            }
        }
        {
            std::lock_guard<std::mutex> q_lock(query_oracle_mutex);
            query_oracle.onBatchCellsOpened(newly_opened_cells, observed);
            for (const auto& pt : newly_opened_cells) {
                query_map.setObserved(pt, CellState::Free);
            }
        }
        double_buffer_dirty = false;
        auto t1 = std::chrono::high_resolution_clock::now();
        last_htrip_update_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

        // Safe decay glow update under decay_mutex
        {
            std::lock_guard<std::mutex> lock(decay_mutex);
            for (int l_idx : updated_leaf_indices) {
                if (l_idx >= 0 && static_cast<size_t>(l_idx) < active_tile_decay.size()) {
                    active_tile_decay[static_cast<size_t>(l_idx)] = 0;
                }
                if (l_idx >= 0 && static_cast<size_t>(l_idx) < active_index_decay.size()) {
                    active_index_decay[static_cast<size_t>(l_idx)] = 20; // Index updated: Electric emerald glow
                }
            }
        }
        newly_opened_cells.clear();
        return;
    }

    // In asynchronous mode:
    // If background thread is currently busy updating, buffer newly opened cells and return immediately!
    if (is_updating_htrip.load()) {
        if (!newly_opened_cells.empty()) {
            updateDirtyTileIndices(newly_opened_cells);
            async_pending_cells.insert(async_pending_cells.end(), newly_opened_cells.begin(), newly_opened_cells.end());
            newly_opened_cells.clear();
        }
        return;
    }

    // Worker is idle: collect all pending cells into a batch
    std::vector<Point> batch;
    batch.reserve(newly_opened_cells.size() + async_pending_cells.size());
    if (!async_pending_cells.empty()) {
        batch.insert(batch.end(), async_pending_cells.begin(), async_pending_cells.end());
        async_pending_cells.clear();
    }
    if (!newly_opened_cells.empty()) {
        batch.insert(batch.end(), newly_opened_cells.begin(), newly_opened_cells.end());
        newly_opened_cells.clear();
    }

    if (batch.empty()) return;

    updateDirtyTileIndices(batch);
    last_revealed_cells_count = static_cast<int>(batch.size());

    total_map_updates++;
    uint64_t ticks_since = total_ticks - last_update_tick;
    last_update_tick = total_ticks;
    if (ticks_since > 0) {
        double inst_hz = 1.0 / (static_cast<double>(ticks_since) * 0.05);
        map_update_freq_hz = (map_update_freq_hz <= 0.0) ? inst_hz : (0.85 * map_update_freq_hz + 0.15 * inst_hz);
    }

    if (update_worker_thread.joinable()) {
        update_worker_thread.join();
    }

    is_updating_htrip.store(true);
    update_worker_thread = std::thread([this, batch_cells = std::move(batch), dirty_indices = this->updated_leaf_indices, obs_snapshot = this->observed]() {
        auto t0 = std::chrono::high_resolution_clock::now();

        // 1. Update the background update buffer (queries continue concurrently on front buffer!)
        MinPlusOracle* u_oracle = this->active_update_oracle_ptr.load(std::memory_order_relaxed);
        GridMap* u_map = this->active_update_map_ptr.load(std::memory_order_relaxed);
        {
            std::lock_guard<std::mutex> lock(this->oracle_mutex);
            u_oracle->onBatchCellsOpened(batch_cells, obs_snapshot);
            for (const auto& pt : batch_cells) {
                u_map->setObserved(pt, CellState::Free);
            }
        }
        {
            std::lock_guard<std::mutex> lock(this->commit_mutex);
            this->completed_indexed_cells.insert(
                this->completed_indexed_cells.end(),
                batch_cells.begin(),
                batch_cells.end()
            );
        }

        // 2. TRUE BUFFER SWAP: Atomic pointer exchange!
        // Instantaneous pointer swap: query worker immediately sees the updated buffer!
        MinPlusOracle* old_q_oracle = nullptr;
        GridMap* old_q_map = nullptr;
        {
            std::lock_guard<std::mutex> q_lock(this->query_oracle_mutex);
            old_q_oracle = this->active_query_oracle_ptr.load(std::memory_order_relaxed);
            old_q_map = this->active_query_map_ptr.load(std::memory_order_relaxed);

            this->active_query_oracle_ptr.store(u_oracle, std::memory_order_release);
            this->active_query_map_ptr.store(u_map, std::memory_order_release);

            this->active_update_oracle_ptr.store(old_q_oracle, std::memory_order_release);
            this->active_update_map_ptr.store(old_q_map, std::memory_order_release);
        }

        // 3. Catch up the previous query buffer (now the back buffer) completely OFF-LOCK!
        // The query worker has already transitioned to u_oracle; this catch-up runs fully asynchronously!
        if (old_q_oracle != nullptr && old_q_map != nullptr) {
            old_q_oracle->onBatchCellsOpened(batch_cells, obs_snapshot);
            for (const auto& pt : batch_cells) {
                old_q_map->setObserved(pt, CellState::Free);
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        this->last_htrip_update_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

        // Safe decay glow update under decay_mutex (no race with main thread)
        {
            std::lock_guard<std::mutex> lock(this->decay_mutex);
            for (int l_idx : dirty_indices) {
                if (l_idx >= 0 && static_cast<size_t>(l_idx) < this->active_tile_decay.size()) {
                    this->active_tile_decay[static_cast<size_t>(l_idx)] = 0;
                }
                if (l_idx >= 0 && static_cast<size_t>(l_idx) < this->active_index_decay.size()) {
                    this->active_index_decay[static_cast<size_t>(l_idx)] = 20; // Index updated: Electric emerald glow
                }
            }
        }

        this->is_updating_htrip.store(false);
    });
}

std::vector<Point> MultiRobotSimulator::planPathAStar(Point start, Point goal) const {
    std::vector<Point> path;
    if (!observed.inBounds(start) || !observed.inBounds(goal)) return path;
    if (start == goal) {
        path.push_back(start);
        return path;
    }

    // If goal itself is not passable (e.g. unknown cell or overlapping obstacle edge), snap to nearest passable cell
    if (!observed.isPassable(goal)) {
        Point best_alt = goal;
        double best_dist = 1e9;
        for (int dr = -2; dr <= 2; ++dr) {
            for (int dc = -2; dc <= 2; ++dc) {
                int gr = goal.r + dr;
                int gc = goal.c + dc;
                if (!observed.inBounds(gr, gc) || !observed.isPassable(gr, gc)) continue;
                double d = std::hypot(dr, dc);
                if (d < best_dist) {
                    best_dist = d;
                    best_alt = {gr, gc};
                }
            }
        }
        if (observed.isPassable(best_alt)) {
            goal = best_alt;
        }
    }

    // If start itself is not passable (e.g. unknown cell or slightly overlapping obstacle edge), snap to nearest passable cell
    if (!observed.isPassable(start)) {
        Point best_start = start;
        double best_dist = 1e9;
        for (int dr = -2; dr <= 2; ++dr) {
            for (int dc = -2; dc <= 2; ++dc) {
                int sr = start.r + dr;
                int sc = start.c + dc;
                if (observed.inBounds(sr, sc) && observed.isPassable(sr, sc)) {
                    double d = std::hypot(dr, dc);
                    if (d < best_dist) {
                        best_dist = d;
                        best_start = {sr, sc};
                    }
                }
            }
        }
        if (observed.isPassable(best_start)) {
            start = best_start;
        }
    }

    if (!observed.isPassable(start) || !observed.isPassable(goal)) return path;

    int R = observed.rows;
    int C = observed.cols;
    int total_cells = R * C;

    static thread_local std::vector<float> g_score;
    static thread_local std::vector<int> parent;
    static thread_local std::vector<int> touched;

    if (g_score.size() < static_cast<size_t>(total_cells)) {
        g_score.assign(static_cast<size_t>(total_cells), 1e9f);
        parent.assign(static_cast<size_t>(total_cells), -1);
    }
    touched.clear();

    struct Node {
        int idx;
        float f;
        bool operator>(const Node& other) const noexcept { return f > other.f; }
    };
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> open_set;

    int start_idx = observed.index(start);
    int goal_idx = observed.index(goal);

    auto heuristic = [&](int r, int c) noexcept -> float {
        float dr = std::abs(r - goal.r);
        float dc = std::abs(c - goal.c);
        return (dr + dc) + (1.41421356f - 2.0f) * std::min(dr, dc);
    };

    g_score[static_cast<size_t>(start_idx)] = 0.0f;
    parent[static_cast<size_t>(start_idx)] = start_idx;
    touched.push_back(start_idx);
    open_set.push({start_idx, heuristic(start.r, start.c)});

    static const int dr[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
    static const int dc[8] = {0, 0, -1, 1, -1, 1, -1, 1};
    static const float move_cost[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.4142f, 1.4142f, 1.4142f, 1.4142f};

    bool reached = false;

    while (!open_set.empty()) {
        auto [u, f_u] = open_set.top();
        open_set.pop();

        if (u == goal_idx) {
            reached = true;
            break;
        }

        float g_u = g_score[static_cast<size_t>(u)];
        int ur = u / C;
        int uc = u % C;

        // Closed set / duplicate-pop check: skip stale queue entries
        if (f_u > g_u + heuristic(ur, uc) + 1e-4f) continue;

        for (int i = 0; i < 8; ++i) {
            int vr = ur + dr[i];
            int vc = uc + dc[i];
            if (!observed.inBounds(vr, vc)) continue;
            if (!observed.isPassable(vr, vc)) continue;

            // Diagonal corner-cutting prevention: both adjacent orthogonal cells must be passable
            if (i >= 4) {
                int adj1_r = ur + dr[i];
                int adj1_c = uc;
                int adj2_r = ur;
                int adj2_c = uc + dc[i];
                if (!observed.isPassable(adj1_r, adj1_c) || !observed.isPassable(adj2_r, adj2_c)) {
                    continue;
                }
            }

            uint8_t clr = observed.getClearance(vr, vc);
            // Clearance penalty: strongly favor wide-open corridors and give corner vertices wide berth
            float pen = (clr <= 1) ? 35.0f : ((clr == 2) ? 8.0f : 0.0f);
            float tentative_g = g_u + move_cost[i] + pen;

            int v = vr * C + vc;
            if (tentative_g < g_score[static_cast<size_t>(v)]) {
                if (g_score[static_cast<size_t>(v)] >= 1e8f) {
                    touched.push_back(v);
                }
                g_score[static_cast<size_t>(v)] = tentative_g;
                parent[static_cast<size_t>(v)] = u;
                open_set.push({v, tentative_g + heuristic(vr, vc)});
            }
        }
    }

    if (reached && goal_idx != start_idx) {
        int curr = goal_idx;
        while (curr != start_idx && curr >= 0) {
            path.push_back({curr / C, curr % C});
            int p = parent[static_cast<size_t>(curr)];
            if (p == curr) break;
            curr = p;
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
    }

    // Reset touched entries for next search
    for (int idx : touched) {
        g_score[static_cast<size_t>(idx)] = 1e9f;
        parent[static_cast<size_t>(idx)] = -1;
    }

    return path;
}

Point MultiRobotSimulator::findNextStepToGoal(Point start, Point goal) const {
    auto path = planPathAStar(start, goal);
    if (path.size() >= 2) return path[1];
    if (!path.empty()) return path[0];
    return start;
}

bool MultiRobotSimulator::hasLineOfSight(Point p1, Point p2) const noexcept {
    int dr = std::abs(p2.r - p1.r);
    int dc = std::abs(p2.c - p1.c);
    int sr = (p1.r < p2.r) ? 1 : -1;
    int sc = (p1.c < p2.c) ? 1 : -1;
    int err = dr - dc;

    int r = p1.r;
    int c = p1.c;
    while (true) {
        if (!ground_truth.inBounds(r, c)) return false;
        if (ground_truth.getObserved(r, c) == CellState::Obstacle) return false;
        if (r == p2.r && c == p2.c) break;

        int e2 = 2 * err;
        int prev_r = r;
        int prev_c = c;
        if (e2 > -dc) {
            err -= dc;
            r += sr;
        }
        if (e2 < dr) {
            err += dr;
            c += sc;
        }

        // Prevent diagonal corner-cutting against adjacent obstacles
        if (r != prev_r && c != prev_c) {
            if (ground_truth.getObserved(prev_r, c) == CellState::Obstacle ||
                ground_truth.getObserved(r, prev_c) == CellState::Obstacle) {
                return false;
            }
        }
    }
    return true;
}

bool MultiRobotSimulator::isFrontierAreaExplored(Point p, int radius_cells) const noexcept {
    if (p.r < 0 || p.c < 0 || !observed.inBounds(p)) return true;
    if (!enable_fog_of_war) {
        for (const auto& bot : robots) {
            Point bot_cell = worldToCell(bot.x, bot.y);
            double d = std::hypot(bot_cell.r - p.r, bot_cell.c - p.c);
            if (d <= static_cast<double>(radius_cells + 2)) {
                return true;
            }
        }
        return false;
    }
    for (int dr = -radius_cells; dr <= radius_cells; ++dr) {
        for (int dc = -radius_cells; dc <= radius_cells; ++dc) {
            int r = p.r + dr;
            int c = p.c + dc;
            if (observed.inBounds(r, c)) {
                if (observed.getObserved(r, c) == CellState::Unknown) {
                    return false;
                }
            }
        }
    }
    return true;
}

Point MultiRobotSimulator::findBestReachableFrontier(const RobotState& bot, std::vector<Point>& out_path) const {
    out_path.clear();
    if (viewpoints.empty()) return {-1, -1};

    Point cur_pt = worldToCell(bot.x, bot.y);

    // Collect other robots' current goals to prioritize unshared frontiers
    std::vector<Point> other_goals;
    other_goals.reserve(robots.size());
    for (const auto& other : robots) {
        if (other.id != bot.id && other.has_goal && other.assigned_frontier.r >= 0) {
            other_goals.push_back(other.assigned_frontier);
        }
    }

    struct Candidate {
        Point pt;
        double dist_sq;
        bool is_shared;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(viewpoints.size());

    for (const auto& vp : viewpoints) {
        if (!observed.inBounds(vp)) continue;
        double dr = vp.r - cur_pt.r;
        double dc = vp.c - cur_pt.c;
        double d2 = dr * dr + dc * dc;

        bool shared = false;
        for (const auto& og : other_goals) {
            if (std::hypot(vp.r - og.r, vp.c - og.c) <= 6.0) {
                shared = true;
                break;
            }
        }
        candidates.push_back({vp, d2, shared});
    }

    // Sort: unshared first, then closest distance
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        if (a.is_shared != b.is_shared) return !a.is_shared;
        return a.dist_sq < b.dist_sq;
    });

    // Check reachability via planPathAStar on top candidates until a valid path is found
    for (const auto& cand : candidates) {
        auto path = planPathAStar(cur_pt, cand.pt);
        if (!path.empty()) {
            out_path = std::move(path);
            return cand.pt;
        }
    }

    return {-1, -1};
}

void MultiRobotSimulator::applyFleetTours(std::vector<RobotTour>& fleet_tours, const std::vector<Point>& current_vps) {
    size_t N = robots.size();
    for (size_t i = 0; i < N; ++i) {
        if (i < fleet_tours.size() && !fleet_tours[i].empty()) {
            Point current_goal = robots[i].assigned_frontier;
            bool area_explored = isFrontierAreaExplored(current_goal, 4);

            // Check if current goal is still a valid active or shifted frontier
            bool has_active_nearby_vp = false;
            if (current_goal.r >= 0 && current_goal.c >= 0) {
                for (const auto& vp : current_vps) {
                    if (std::hypot(vp.r - current_goal.r, vp.c - current_goal.c) <= 6.0) {
                        has_active_nearby_vp = true;
                        break;
                    }
                }
            }

            // Persistence commitment policy:
            // Robot must persist on its assigned frontier if:
            // - Robot has an active goal
            // - Goal area is not yet fully explored
            // - Either the persistence window is active (> 0) OR robot is making progress (stuck_ticks < 8)
            //   and the frontier is still active nearby
            bool should_persist = robots[i].has_goal &&
                                  current_goal.r >= 0 && current_goal.c >= 0 &&
                                  observed.inBounds(current_goal) &&
                                  !area_explored &&
                                  (has_active_nearby_vp || robots[i].persistence_ticks > 0) &&
                                  robots[i].stuck_ticks < 8;

            Point prev_assigned = robots[i].assigned_frontier;
            robots[i].tour = std::move(fleet_tours[i]);
            robots[i].has_goal = true;

            if (should_persist) {
                // Retain commitment: pin current_goal to the front of this robot's tour
                if (robots[i].tour.currentGoal() != current_goal) {
                    robots[i].tour.waypoints.insert(robots[i].tour.waypoints.begin(), current_goal);
                }
                robots[i].assigned_frontier = current_goal;
            } else {
                robots[i].assigned_frontier = robots[i].tour.currentGoal();
                robots[i].persistence_ticks = persistence_window_ticks;
            }

            if (robots[i].assigned_frontier != prev_assigned) {
                robots[i].planned_path.clear();
                robots[i].path_index = 0;
                robots[i].persistence_ticks = persistence_window_ticks;
            }

            robots[i].assigned_frontier_id = -1;
            for (size_t f = 0; f < frontiers.size(); ++f) {
                if (frontiers[f].representative == robots[i].assigned_frontier) {
                    robots[i].assigned_frontier_id = static_cast<int>(f);
                    break;
                }
            }
        } else {
            // No tour assigned in m-TSP: immediately find a reachable frontier so robot never sits idle
            std::vector<Point> fallback_path;
            Point best_vp = findBestReachableFrontier(robots[i], fallback_path);
            if (best_vp.r >= 0 && !fallback_path.empty()) {
                robots[i].tour.waypoints = {best_vp};
                robots[i].tour.current_wp_idx = 0;
                robots[i].assigned_frontier = best_vp;
                robots[i].has_goal = true;
                robots[i].planned_path = std::move(fallback_path);
                robots[i].path_index = 0;
                robots[i].next_waypoint = robots[i].planned_path[0];
                robots[i].persistence_ticks = persistence_window_ticks;
                robots[i].waypoint_stuck_count = 0;
            } else {
                robots[i].tour.clear();
                robots[i].has_goal = false;
                robots[i].assigned_frontier = {-1, -1};
                robots[i].assigned_frontier_id = -1;
                robots[i].planned_path.clear();
                robots[i].path_index = 0;
                robots[i].persistence_ticks = 0;
            }
        }
    }
}

static Point snapPointToIndexedMap(const GridMap& map, Point p, int max_radius = 16) noexcept {
    if (map.inBounds(p) && map.isPassable(p)) return p;
    double min_d = 1e9;
    Point best_p = p;
    for (int dr = -max_radius; dr <= max_radius; ++dr) {
        for (int dc = -max_radius; dc <= max_radius; ++dc) {
            Point cand{p.r + dr, p.c + dc};
            if (map.inBounds(cand) && map.isPassable(cand)) {
                double d = static_cast<double>(dr * dr + dc * dc);
                if (d < min_d) {
                    min_d = d;
                    best_p = cand;
                }
            }
        }
    }
    if (min_d < 1e8) return best_p;
    // Fallback: wider scan up to 32 cells if map was sparsely indexed
    for (int dr = -32; dr <= 32; ++dr) {
        for (int dc = -32; dc <= 32; ++dc) {
            Point cand{p.r + dr, p.c + dc};
            if (map.inBounds(cand) && map.isPassable(cand)) {
                double d = static_cast<double>(dr * dr + dc * dc);
                if (d < min_d) {
                    min_d = d;
                    best_p = cand;
                }
            }
        }
    }
    return best_p;
}

/**
 * @brief Asynchronously dispatches H-TRIP all-pairs distance evaluation and m-TSP fleet planning.
 *
 * Execution flow:
 * 1. Snaps current robot poses and candidate frontier viewpoints to query_map under query_oracle_mutex.
 * 2. If benchmarking is active, creates a paired comparison job for BFS/A*.
 * 3. Launches query_worker_thread:
 *    - Evaluates query_oracle.queryAllPairs using caller-owned async_query_scratch (zero heap allocations).
 *    - Records query time and notifies baseline comparison job.
 *    - Executes MTSPPlanner::planFleetTours with 2-opt refinement.
 *    - Stages computed matrix, fleet tours, and viewpoints under query_mutex, setting new_tours_ready = true.
 */
void MultiRobotSimulator::startAsyncPairwiseQuery() {
    if (robots.empty() || viewpoints.empty()) return;
    if (is_querying_pairwise.load()) return;
    is_querying_pairwise.store(true);

    if (query_worker_thread.joinable()) {
        query_worker_thread.join();
    }

    size_t N = robots.size();
    const size_t k_cap = static_cast<size_t>(std::max(1, max_frontiers_eval));
    size_t K = std::min(viewpoints.size(), k_cap);

    std::vector<Point> combined_pts;
    combined_pts.reserve(N + K);
    std::vector<Point> robot_pts;
    robot_pts.reserve(N);
    std::vector<Point> snap_vps;
    snap_vps.reserve(K);

    // Snap robot and viewpoint coordinates to query_map under query_oracle_mutex
    {
        std::lock_guard<std::mutex> lock(query_oracle_mutex);
        const auto* q_map = active_query_map_ptr.load(std::memory_order_acquire);
        for (const auto& bot : robots) {
            Point r_pt = worldToCell(bot.x, bot.y);
            r_pt = snapPointToIndexedMap(*q_map, r_pt, 16);
            robot_pts.push_back(r_pt);
            combined_pts.push_back(r_pt);
        }

        for (size_t i = 0; i < K; ++i) {
            Point q_vp = snapPointToIndexedMap(*q_map, viewpoints[i], 16);
            combined_pts.push_back(q_vp);
            snap_vps.push_back(viewpoints[i]);
        }
    }

    last_query_pts = combined_pts;
    uint64_t bfs_job_id = 0;
    const bool want_bfs = (enable_bfs_baseline || enable_astar_baseline) && (bfs_compare_interval > 0) &&
                          ((total_htrip_queries + 1) % static_cast<uint64_t>(bfs_compare_interval) == 0) &&
                          combined_pts.size() >= 2;

    if (want_bfs && baseline_mode == BaselineExecutionMode::Asynchronous) {
        GridMap bfs_snap;
        {
            std::lock_guard<std::mutex> lock(query_oracle_mutex);
            bfs_snap = *active_query_map_ptr.load(std::memory_order_acquire);
        }
        bfs_job_id = createBfsCompareJob(std::move(bfs_snap), last_query_pts);
    }

    query_worker_thread = std::thread([this, c_pts = std::move(combined_pts), r_pts = std::move(robot_pts), vps = std::move(snap_vps), bfs_job_id, want_bfs]() {
        // H-TRIP all-pairs distance query on the active double-buffered snapshot.
        // Lock query_oracle_mutex ONLY during the query.
        // It NEVER locks oracle_mutex, so ongoing map updates NEVER block this query.
        auto t0 = std::chrono::high_resolution_clock::now();
        {
            std::lock_guard<std::mutex> lock(this->query_oracle_mutex);
            const auto* q_oracle = this->active_query_oracle_ptr.load(std::memory_order_acquire);
            const auto* q_map = this->active_query_map_ptr.load(std::memory_order_acquire);
            q_oracle->queryAllPairs(c_pts, *q_map, this->async_query_matrix, this->async_query_scratch);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double q_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

        if (want_bfs && this->baseline_mode == BaselineExecutionMode::Synchronous) {
            GridMap bfs_snap;
            {
                std::lock_guard<std::mutex> lock(this->query_oracle_mutex);
                bfs_snap = *this->active_query_map_ptr.load(std::memory_order_acquire);
            }
            this->runSyncBaselines(bfs_snap, c_pts, q_us, this->total_ticks);
        } else if (bfs_job_id != 0) {
            this->noteHtripQueryTime(bfs_job_id, q_us);
        }

        // Plan fleet tours using MTSPPlanner without holding any mutex
        std::vector<RobotTour> tours;
        auto t_mtsp0 = std::chrono::high_resolution_clock::now();
        MTSPPlanner::planFleetTours(r_pts, vps, this->async_query_matrix, tours, 32);
        auto t_mtsp1 = std::chrono::high_resolution_clock::now();
        double p_us = std::chrono::duration<double, std::micro>(t_mtsp1 - t_mtsp0).count();

        // Stage computed tours and matrix for simulator thread
        {
            std::lock_guard<std::mutex> lock(this->query_mutex);
            this->staged_all_pairs_matrix = this->async_query_matrix;
            this->staged_fleet_tours = std::move(tours);
            this->staged_viewpoints = std::move(vps);
            this->last_htrip_query_us = q_us;
            this->last_mtsp_plan_us = p_us;
            this->new_tours_ready.store(true);
        }
        this->is_querying_pairwise.store(false);
    });
}

void MultiRobotSimulator::assignFrontiersToFleet() {
    if (robots.empty()) return;
    if (baseline_mode == BaselineExecutionMode::Asynchronous) {
        pumpBfsPending();
    }

    // 0. Adopt fresh fleet tours if background pairwise distance query thread has completed
    if (new_tours_ready.load()) {
        std::vector<RobotTour> tours_to_apply;
        std::vector<Point> vps_to_apply;
        {
            std::lock_guard<std::mutex> lock(query_mutex);
            tours_to_apply = std::move(staged_fleet_tours);
            vps_to_apply = std::move(staged_viewpoints);
            all_pairs_matrix = std::move(staged_all_pairs_matrix);
            new_tours_ready.store(false);
        }
        total_htrip_queries++;
        auto q_done_now = std::chrono::steady_clock::now();
        if (last_query_completed_time.time_since_epoch().count() > 0) {
            double dt_actual_s = std::chrono::duration<double>(q_done_now - last_query_completed_time).count();
            if (dt_actual_s > 0.001 && dt_actual_s < 5.0) {
                double inst_hz = 1.0 / dt_actual_s;
                htrip_query_freq_hz = (htrip_query_freq_hz <= 0.0)
                    ? inst_hz
                    : (0.85 * htrip_query_freq_hz + 0.15 * inst_hz);
            }
        } else {
            htrip_query_freq_hz = target_query_hz;
        }
        last_query_completed_time = q_done_now;
        last_query_tick = total_ticks;
        applyFleetTours(tours_to_apply, vps_to_apply);
        last_replan_tick = total_ticks;
    }

    // 1. Advance tour if robot physically reached current goal or cluster dissolved
    for (auto& bot : robots) {
        if (!bot.has_goal) continue;

        auto [gx, gy] = cellToWorld(bot.assigned_frontier);
        double dist_to_goal_m = std::hypot(bot.x - gx, bot.y - gy);

        // A. Centroid proximity tracking:
        // As LiDAR reveals cells along frontier boundary, the cluster centroid shifts by 1-2 cells.
        // Check if there is an active viewpoint close to bot.assigned_frontier.
        Point best_matching_vp = {-1, -1};
        double min_vp_dist = 1e9;
        for (const auto& vp : viewpoints) {
            double d = std::hypot(vp.r - bot.assigned_frontier.r, vp.c - bot.assigned_frontier.c);
            if (d < min_vp_dist) {
                min_vp_dist = d;
                best_matching_vp = vp;
            }
        }

        // If an active viewpoint is within 6 cells, the cluster is still active (shifted slightly)
        if (min_vp_dist <= 6.0 && best_matching_vp.r >= 0) {
            bot.assigned_frontier = best_matching_vp;
        }

        // B. Reached or fully explored check:
        // Physically reached proximity of goal (within 0.35m or 3.5 cells)
        bool physically_reached = (dist_to_goal_m < std::max(0.35, cell_size_m * 3.5));

        // Frontier area fully explored: no Unknown cells remain in radius
        bool area_explored = isFrontierAreaExplored(bot.assigned_frontier, 4);

        // Only consider the cluster dissolved if no viewpoint exists nearby AND the area has been explored
        bool cluster_dissolved = (min_vp_dist > 6.0) && area_explored;

        if (physically_reached || cluster_dissolved) {
            bot.tour.advance();
            bot.planned_path.clear();
            bot.path_index = 0;
            bot.next_waypoint = {-1, -1};
            bot.waypoint_stuck_count = 0;
            if (!bot.tour.empty()) {
                bot.assigned_frontier = bot.tour.currentGoal();
                bot.has_goal = true;
                bot.persistence_ticks = persistence_window_ticks;
            } else {
                bot.has_goal = false;
                bot.assigned_frontier = {-1, -1};
                bot.assigned_frontier_id = -1;
                bot.persistence_ticks = 0;
            }
        }
    }

    // 1b. Ensure ZERO idle robots ("bosta kalan robot olmamali")
    // If any robot has no active goal or its tour is empty, immediately assign it
    // the best reachable frontier from viewpoints so it NEVER sits idle or wanders aimlessly.
    for (auto& bot : robots) {
        if (!bot.has_goal || bot.tour.empty() || bot.assigned_frontier.r < 0) {
            std::vector<Point> fallback_path;
            Point best_vp = findBestReachableFrontier(bot, fallback_path);
            if (best_vp.r >= 0 && !fallback_path.empty()) {
                bot.tour.waypoints = {best_vp};
                bot.tour.current_wp_idx = 0;
                bot.assigned_frontier = best_vp;
                bot.has_goal = true;
                bot.planned_path = std::move(fallback_path);
                bot.path_index = 0;
                bot.next_waypoint = bot.planned_path[0];
                bot.persistence_ticks = persistence_window_ticks;
                bot.waypoint_stuck_count = 0;
            }
        }
    }

    // 2. Decide if full fleet m-TSP tour re-planning is needed
    bool any_robot_starved = false;
    for (const auto& bot : robots) {
        if (!bot.has_goal || bot.tour.empty()) {
            any_robot_starved = true;
            break;
        }
    }

    uint64_t ticks_since_replan = total_ticks - last_replan_tick;
    bool periodic_replan = (replan_period_ticks <= 1) || (ticks_since_replan >= static_cast<uint64_t>(replan_period_ticks));
    bool need_replan = periodic_replan || (any_robot_starved && !viewpoints.empty());

    if (enable_async_updates) {
        // Continuous pairwise distance calculation on double-buffered previous index
        // Target cadence: 20 Hz (every 50 ms)
        auto now = std::chrono::steady_clock::now();
        double min_interval_ms = (target_query_hz > 0.0) ? (1000.0 / target_query_hz) : 50.0;
        double elapsed_ms = std::chrono::duration<double, std::milli>(now - last_pairwise_dispatch_time).count();

        if (elapsed_ms >= min_interval_ms || last_pairwise_dispatch_time.time_since_epoch().count() == 0) {
            if (!is_querying_pairwise.load() && !viewpoints.empty() && !robots.empty()) {
                last_pairwise_dispatch_time = now;
                startAsyncPairwiseQuery();
            }
        }
    } else {
        // Synchronous pairwise distance matrix query using double-buffered snapshot.
        // In synchronous baseline mode, H-TRIP and BFS execute on EVERY single step!
        const bool should_query = (need_replan && !viewpoints.empty()) || (need_replan && robots.size() >= 2);
        if (should_query) {
            size_t N = robots.size();
            const size_t k_cap = static_cast<size_t>(std::max(1, max_frontiers_eval));
            size_t K = std::min(viewpoints.size(), k_cap);

            std::vector<Point> combined_pts;
            combined_pts.reserve(N + K);
            std::vector<Point> robot_pts;
            robot_pts.reserve(N);

            const auto* q_map_ptr = active_query_map_ptr.load(std::memory_order_acquire);
            for (const auto& bot : robots) {
                Point r_pt = worldToCell(bot.x, bot.y);
                r_pt = snapPointToIndexedMap(*q_map_ptr, r_pt, 16);
                robot_pts.push_back(r_pt);
                combined_pts.push_back(r_pt);
            }

            for (size_t i = 0; i < K; ++i) {
                Point q_vp = snapPointToIndexedMap(*q_map_ptr, viewpoints[i], 16);
                combined_pts.push_back(q_vp);
            }

            last_query_pts = combined_pts;
            uint64_t bfs_job_id = 0;
            const bool want_bfs = (enable_bfs_baseline || enable_astar_baseline) &&
                                  ((bfs_compare_interval == 1) || (bfs_compare_interval > 0 &&
                                  ((total_htrip_queries + 1) % static_cast<uint64_t>(bfs_compare_interval) == 0))) &&
                                  combined_pts.size() >= 2;

            if (want_bfs && baseline_mode == BaselineExecutionMode::Asynchronous) {
                GridMap bfs_snap;
                {
                    std::lock_guard<std::mutex> lock(query_oracle_mutex);
                    bfs_snap = *active_query_map_ptr.load(std::memory_order_acquire);
                }
                bfs_job_id = createBfsCompareJob(std::move(bfs_snap), last_query_pts);
            }

            if (combined_pts.size() >= 2) {
                auto t0 = std::chrono::high_resolution_clock::now();
                {
                    std::lock_guard<std::mutex> lock(query_oracle_mutex);
                    const auto* q_oracle = active_query_oracle_ptr.load(std::memory_order_acquire);
                    const auto* q_map = active_query_map_ptr.load(std::memory_order_acquire);
                    q_oracle->queryAllPairs(combined_pts, *q_map, all_pairs_matrix, sync_query_scratch);
                }
                auto t1 = std::chrono::high_resolution_clock::now();
                double q_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
                last_htrip_query_us = q_us;
                total_htrip_queries++;

                if (want_bfs && baseline_mode == BaselineExecutionMode::Synchronous) {
                    GridMap bfs_snap;
                    {
                        std::lock_guard<std::mutex> lock(query_oracle_mutex);
                        bfs_snap = *active_query_map_ptr.load(std::memory_order_acquire);
                    }
                    runSyncBaselines(bfs_snap, combined_pts, q_us, total_ticks);
                } else if (bfs_job_id != 0) {
                    noteHtripQueryTime(bfs_job_id, q_us);
                }
            }

            if (need_replan && !viewpoints.empty() && combined_pts.size() >= 2) {
                std::vector<RobotTour> fleet_tours;
                auto t_mtsp0 = std::chrono::high_resolution_clock::now();
                MTSPPlanner::planFleetTours(robot_pts, viewpoints, all_pairs_matrix, fleet_tours, 32);
                auto t_mtsp1 = std::chrono::high_resolution_clock::now();
                last_mtsp_plan_us = std::chrono::duration<double, std::micro>(t_mtsp1 - t_mtsp0).count();

                applyFleetTours(fleet_tours, viewpoints);
                last_replan_tick = total_ticks;
            }
        }
    }

    // 3. Update planned A* path and next lookahead waypoint towards current assigned goal
    for (auto& bot : robots) {
        if (bot.has_goal && observed.inBounds(bot.assigned_frontier)) {
            Point cur_pt = worldToCell(bot.x, bot.y);

            bool need_replan = bot.planned_path.empty() ||
                               (total_ticks % 10 == 0) ||
                               (bot.path_index >= bot.planned_path.size());

            // Check if upcoming path steps hit newly discovered obstacle
            if (!need_replan && !bot.planned_path.empty()) {
                size_t check_end = std::min(bot.planned_path.size(), bot.path_index + 12);
                for (size_t i = bot.path_index; i < check_end; ++i) {
                    if (!observed.isPassable(bot.planned_path[i])) {
                        need_replan = true;
                        break;
                    }
                }
            }

            if (need_replan) {
                bot.planned_path = planPathAStar(cur_pt, bot.assigned_frontier);
                bot.path_index = 0;
            }

            if (bot.planned_path.empty()) {
                // Goal unreachable in observed grid: advance tour or find immediate reachable frontier!
                bot.tour.advance();
                if (!bot.tour.empty()) {
                    bot.assigned_frontier = bot.tour.currentGoal();
                    bot.planned_path = planPathAStar(cur_pt, bot.assigned_frontier);
                    bot.path_index = 0;
                }
                if (bot.planned_path.empty()) {
                    std::vector<Point> fb_path;
                    Point best_vp = findBestReachableFrontier(bot, fb_path);
                    if (best_vp.r >= 0 && !fb_path.empty()) {
                        bot.tour.waypoints = {best_vp};
                        bot.tour.current_wp_idx = 0;
                        bot.assigned_frontier = best_vp;
                        bot.has_goal = true;
                        bot.planned_path = std::move(fb_path);
                        bot.path_index = 0;
                        bot.next_waypoint = bot.planned_path[0];
                        bot.persistence_ticks = persistence_window_ticks;
                        bot.waypoint_stuck_count = 0;
                    } else {
                        bot.has_goal = false;
                        bot.assigned_frontier = {-1, -1};
                        bot.assigned_frontier_id = -1;
                        bot.next_waypoint = {-1, -1};
                    }
                }
            }

            if (!bot.planned_path.empty()) {
                double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
                double reach_dist = std::max(0.18, cell_size_m * 2.5);

                // Advance path_index past waypoints already reached
                while (bot.path_index < bot.planned_path.size()) {
                    auto [wpx, wpy] = cellToWorld(bot.planned_path[bot.path_index]);
                    double d = std::hypot(bot.x - wpx, bot.y - wpy);
                    if (d < reach_dist && bot.path_index + 1 < bot.planned_path.size()) {
                        bot.path_index++;
                    } else {
                        break;
                    }
                }

                // Check if final frontier physically reached or fully explored
                auto [gx, gy] = cellToWorld(bot.assigned_frontier);
                double dist_to_goal = std::hypot(bot.x - gx, bot.y - gy);
                bool goal_reached = (dist_to_goal < std::max(0.35, cell_size_m * 3.5));
                bool area_explored = isFrontierAreaExplored(bot.assigned_frontier, 4);

                if (goal_reached || area_explored) {
                    bot.tour.advance();
                    bot.planned_path.clear();
                    bot.path_index = 0;
                    bot.next_waypoint = {-1, -1};
                    bot.waypoint_stuck_count = 0;
                    if (!bot.tour.empty()) {
                        bot.assigned_frontier = bot.tour.currentGoal();
                        bot.persistence_ticks = persistence_window_ticks;
                    } else {
                        // Immediately find another reachable frontier if available
                        std::vector<Point> next_path;
                        Point best_vp = findBestReachableFrontier(bot, next_path);
                        if (best_vp.r >= 0 && !next_path.empty()) {
                            bot.tour.waypoints = {best_vp};
                            bot.tour.current_wp_idx = 0;
                            bot.assigned_frontier = best_vp;
                            bot.has_goal = true;
                            bot.planned_path = std::move(next_path);
                            bot.path_index = 0;
                            bot.next_waypoint = bot.planned_path[0];
                            bot.persistence_ticks = persistence_window_ticks;
                        } else {
                            bot.has_goal = false;
                            bot.assigned_frontier = {-1, -1};
                            bot.assigned_frontier_id = -1;
                            bot.persistence_ticks = 0;
                        }
                    }
                } else if (bot.path_index >= bot.planned_path.size()) {
                    // Path finished but goal not reached yet: replan path from current position
                    bot.planned_path = planPathAStar(cur_pt, bot.assigned_frontier);
                    bot.path_index = 0;
                    if (bot.planned_path.empty()) {
                        std::vector<Point> next_path;
                        Point best_vp = findBestReachableFrontier(bot, next_path);
                        if (best_vp.r >= 0 && !next_path.empty()) {
                            bot.tour.waypoints = {best_vp};
                            bot.tour.current_wp_idx = 0;
                            bot.assigned_frontier = best_vp;
                            bot.has_goal = true;
                            bot.planned_path = std::move(next_path);
                            bot.path_index = 0;
                            bot.next_waypoint = bot.planned_path[0];
                            bot.persistence_ticks = persistence_window_ticks;
                        } else {
                            bot.has_goal = false;
                            bot.assigned_frontier = {-1, -1};
                            bot.next_waypoint = {-1, -1};
                        }
                    } else {
                        bot.next_waypoint = bot.planned_path[0];
                    }
                } else {
                    // Shortcut lookahead: Find furthest waypoint with direct clear line-of-sight
                    size_t lookahead_steps = static_cast<size_t>(std::max(1, static_cast<int>(std::round(0.40 / cell_size_m))));
                    size_t max_idx = std::min(bot.planned_path.size() - 1, bot.path_index + lookahead_steps);
                    size_t target_idx = bot.path_index;
                    for (size_t idx = bot.path_index; idx <= max_idx; ++idx) {
                        auto [wpx, wpy] = cellToWorld(bot.planned_path[idx]);
                        if (isSegmentClear(bot.x, bot.y, wpx, wpy, eff_r)) {
                            target_idx = idx;
                        } else {
                            break;
                        }
                    }
                    if (target_idx > bot.path_index) {
                        bot.path_index = target_idx;
                    }
                    bot.next_waypoint = bot.planned_path[bot.path_index];
                }
            }
        } else {
            // Robot has no goal or goal out of bounds: immediately find reachable frontier
            std::vector<Point> fb_path;
            Point best_vp = findBestReachableFrontier(bot, fb_path);
            if (best_vp.r >= 0 && !fb_path.empty()) {
                bot.tour.waypoints = {best_vp};
                bot.tour.current_wp_idx = 0;
                bot.assigned_frontier = best_vp;
                bot.has_goal = true;
                bot.planned_path = std::move(fb_path);
                bot.path_index = 0;
                bot.next_waypoint = bot.planned_path[0];
                bot.persistence_ticks = persistence_window_ticks;
                bot.waypoint_stuck_count = 0;
            } else {
                bot.planned_path.clear();
                bot.path_index = 0;
                bot.next_waypoint = {-1, -1};
            }
        }
    }
}

/**
 * @brief Computes continuous forward linear velocity v and angular velocity omega.
 *
 * Implements a 3-layer guidance and control policy:
 * 1. Recovery State Machine:
 *    - BackUp: Reverse linear velocity (-0.35 * v_max) to clear corner/wall pinches.
 *    - RepulsiveSpin: In-place rotation away from the closest obstacle normal detected by LiDAR.
 *    - ReplanNudge: Re-computes A* corridor path from current pose to dislodge path followers.
 * 2. Obstacle-Repulsive Vector Fields:
 *    - Evaluates repulsive forces from all 360 LiDAR hits within range, steering the robot toward open corridors.
 * 3. Carrot-Following Path Tracking:
 *    - Steers toward a forward lookahead waypoint along planned_path with line-of-sight shortcutting.
 */
void MultiRobotSimulator::computeAutonomousControls(RobotState& robot, double dt) {
    (void)dt;
    if (robot.scan.empty()) return;

    // 0. Active Recovery Behaviors (Nav2-style BackUp & Repulsive Spin)
    if (robot.recovery_phase == RobotState::RecoveryPhase::BackUp) {
        robot.recovery_ticks--;
        robot.v = -0.35 * robot_config.v_max; // Smooth reverse away from corner/wall
        robot.omega = 0.0;
        if (robot.recovery_ticks <= 0) {
            // Next stage: Spin in place away from closest obstacle normal
            robot.recovery_phase = RobotState::RecoveryPhase::RepulsiveSpin;
            robot.recovery_ticks = 15; // ~0.75s spin
        }
        return;
    }

    if (robot.recovery_phase == RobotState::RecoveryPhase::RepulsiveSpin) {
        robot.recovery_ticks--;
        robot.v = 0.0;

        // Find angle to closest obstacle from 360 scan
        double min_obs_d = robot_config.lidar.range_m;
        double min_obs_ang = robot.theta;
        for (const auto& beam : robot.scan) {
            if (beam.distance_m < min_obs_d) {
                min_obs_d = beam.distance_m;
                min_obs_ang = beam.angle_rad;
            }
        }
        // Desired escape heading points directly away from closest obstacle
        double repulse_heading = normalizeAngle(min_obs_ang + M_PI);
        double err = normalizeAngle(repulse_heading - robot.theta);
        robot.omega = (err > 0.0 ? 1.0 : -1.0) * robot_config.omega_max;

        if (robot.recovery_ticks <= 0 || std::abs(err) < 0.20) {
            // Recovery complete! Re-plan fresh path from current open position
            robot.recovery_phase = RobotState::RecoveryPhase::None;
            robot.recovery_ticks = 0;
            robot.stuck_ticks = 0;
            if (robot.has_goal) {
                Point cur_pt = worldToCell(robot.x, robot.y);
                robot.planned_path = planPathAStar(cur_pt, robot.assigned_frontier);
                robot.path_index = 0;
                if (!robot.planned_path.empty()) {
                    double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
                    size_t lookahead_steps = static_cast<size_t>(std::max(1, static_cast<int>(std::round(0.35 / cell_size_m))));
                    size_t max_idx = std::min(robot.planned_path.size() - 1, lookahead_steps);
                    size_t target_idx = 0;
                    for (size_t idx = 0; idx <= max_idx; ++idx) {
                        auto [wpx, wpy] = cellToWorld(robot.planned_path[idx]);
                        if (isSegmentClear(robot.x, robot.y, wpx, wpy, eff_r)) target_idx = idx;
                        else break;
                    }
                    robot.next_waypoint = robot.planned_path[target_idx];
                }
            }
        }
        return;
    }

    if (robot.recovery_phase == RobotState::RecoveryPhase::ReplanNudge) {
        robot.recovery_ticks--;
        if (robot.recovery_ticks <= 0) {
            robot.recovery_phase = RobotState::RecoveryPhase::None;
        }
    }

    // Re-plan nudge ("rotayı tekrar bir dürt"):
    // If robot is pursuing a goal and has stalled at a corner for >= 6 ticks (~0.3s):
    if (robot.has_goal && robot.stuck_ticks >= 6 && robot.recovery_phase == RobotState::RecoveryPhase::None) {
        Point cur_pt = worldToCell(robot.x, robot.y);
        robot.planned_path = planPathAStar(cur_pt, robot.assigned_frontier);
        robot.path_index = 0;
        if (!robot.planned_path.empty()) {
            double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
            size_t lookahead_steps = static_cast<size_t>(std::max(1, static_cast<int>(std::round(0.35 / cell_size_m))));
            size_t max_idx = std::min(robot.planned_path.size() - 1, lookahead_steps);
            size_t target_idx = 0;
            for (size_t idx = 0; idx <= max_idx; ++idx) {
                auto [wpx, wpy] = cellToWorld(robot.planned_path[idx]);
                if (isSegmentClear(robot.x, robot.y, wpx, wpy, eff_r)) target_idx = idx;
                else break;
            }
            robot.next_waypoint = robot.planned_path[target_idx];
        }
        robot.recovery_phase = RobotState::RecoveryPhase::ReplanNudge;
        robot.recovery_ticks = 8;
    }

    // Fallback in-place rotation escape (if trail was empty)
    if (robot.escape_ticks > 0) {
        robot.escape_ticks--;
        robot.v = 0.0;
        robot.omega = robot.escape_turn_dir * robot_config.omega_max;
        return;
    }

    // 1. Path Following Controller (Zero Circling, Continuous Alignment)
    if (robot.has_goal && robot.next_waypoint.r >= 0 && robot.next_waypoint.c >= 0) {
        auto [wx, wy] = cellToWorld(robot.next_waypoint);
        double dx = wx - robot.x;
        double dy = wy - robot.y;
        double dist_to_target = std::hypot(dx, dy);
        double to_target_angle = std::atan2(dy, dx);
        double heading_err = normalizeAngle(to_target_angle - robot.theta);

        // Proportional angular steering
        robot.omega = std::clamp(3.5 * heading_err, -robot_config.omega_max, robot_config.omega_max);

        // Forward speed regulation to eliminate circling:
        // A. Heading alignment: If heading error is large (> 60 deg), stop forward drive (v = 0)
        //    to rotate cleanly in place and align with the path segment before driving forward.
        double abs_herr = std::abs(heading_err);
        double heading_align = 0.0;
        if (abs_herr < (M_PI * 0.35)) { // within ~63 degrees
            heading_align = std::cos(abs_herr * 1.4); // 1.0 at 0 deg, drops to 0 at 63 deg
            heading_align = std::max(0.0, heading_align);
        }

        // B. Proximity governor: Prevent unicycle orbiting singularity!
        // Turning radius R = v / omega. To ensure the robot turns inside dist_to_target without
        // orbiting around it, limit forward speed when approaching the target.
        double max_turn_v = std::max(0.08, dist_to_target * robot_config.omega_max * 0.85);
        double target_v = std::min(robot_config.v_max * heading_align, max_turn_v);

        // C. Forward clearance from LiDAR:
        double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
        double fwd_clearance = robot_config.lidar.range_m;
        for (const auto& beam : robot.scan) {
            double rel_ang = normalizeAngle(beam.angle_rad - robot.theta);
            if (std::abs(rel_ang) < (M_PI * 0.45)) {
                double fwd_d = beam.distance_m * std::cos(rel_ang);
                double lat_d = std::abs(beam.distance_m * std::sin(rel_ang));
                if (lat_d < (eff_r * 1.25) && fwd_d > 0.0) {
                    fwd_clearance = std::min(fwd_clearance, fwd_d);
                }
            }
        }
        for (const auto& other : robots) {
            if (other.id == robot.id) continue;
            double od = std::hypot(other.x - robot.x, other.y - robot.y);
            double o_ang = std::atan2(other.y - robot.y, other.x - robot.x);
            double rel_o = normalizeAngle(o_ang - robot.theta);
            if (std::abs(rel_o) < (M_PI * 0.45)) {
                double fwd_d = od * std::cos(rel_o);
                double lat_d = std::abs(od * std::sin(rel_o));
                if (lat_d < (eff_r * 1.5) && fwd_d > 0.0) {
                    fwd_clearance = std::min(fwd_clearance, std::max(0.0, fwd_d - eff_r * 2.0));
                }
            }
        }

        if (fwd_clearance < (cell_size_m * 1.2)) {
            robot.v = 0.0;
        } else if (fwd_clearance < (cell_size_m * 2.5)) {
            double clr_scale = (fwd_clearance - cell_size_m * 1.2) / (cell_size_m * 1.3);
            robot.v = target_v * clr_scale;
        } else {
            robot.v = target_v;
        }
        return;
    }

    // 2. If the robot has no active goal or valid waypoint, try to assign the best reachable frontier!
    if (!viewpoints.empty()) {
        std::vector<Point> fallback_path;
        Point best_vp = findBestReachableFrontier(robot, fallback_path);
        if (best_vp.r >= 0 && !fallback_path.empty()) {
            robot.tour.waypoints = {best_vp};
            robot.tour.current_wp_idx = 0;
            robot.assigned_frontier = best_vp;
            robot.has_goal = true;
            robot.planned_path = std::move(fallback_path);
            robot.path_index = 0;
            robot.next_waypoint = robot.planned_path[0];
            robot.persistence_ticks = persistence_window_ticks;
            robot.waypoint_stuck_count = 0;

            // Immediately steer along this path
            auto [wx, wy] = cellToWorld(robot.next_waypoint);
            double dx = wx - robot.x;
            double dy = wy - robot.y;
            double dist_to_target = std::hypot(dx, dy);
            double to_target_angle = std::atan2(dy, dx);
            double heading_err = normalizeAngle(to_target_angle - robot.theta);

            robot.omega = std::clamp(3.5 * heading_err, -robot_config.omega_max, robot_config.omega_max);
            double abs_herr = std::abs(heading_err);
            double heading_align = (abs_herr < (M_PI * 0.35)) ? std::max(0.0, std::cos(abs_herr * 1.4)) : 0.0;
            double max_turn_v = std::max(0.08, dist_to_target * robot_config.omega_max * 0.85);
            robot.v = std::min(robot_config.v_max * heading_align, max_turn_v);
            return;
        }
    }

    double target_x = robot.x + std::cos(robot.theta) * 1.0;
    double target_y = robot.y + std::sin(robot.theta) * 1.0;
    bool aimed_at_frontier = false;
    if (robot.has_goal && observed.inBounds(robot.assigned_frontier)) {
        const auto [gx, gy] = cellToWorld(robot.assigned_frontier);
        target_x = gx;
        target_y = gy;
        aimed_at_frontier = true;
    } else if (!viewpoints.empty()) {
        double best_d2 = 1e100;
        for (const auto& vp : viewpoints) {
            const auto [fx, fy] = cellToWorld(vp);
            const double dx = fx - robot.x;
            const double dy = fy - robot.y;
            const double d2 = dx * dx + dy * dy;
            if (d2 < best_d2) {
                best_d2 = d2;
                target_x = fx;
                target_y = fy;
                aimed_at_frontier = true;
            }
        }
    }
    if (!aimed_at_frontier) {
        double best_score_dist = 0.0;
        double best_beam_dist = 0.5;
        double best_beam_angle = robot.theta;
        for (const auto& beam : robot.scan) {
            double rel_ang = std::abs(normalizeAngle(beam.angle_rad - robot.theta));
            double dir_bias = (rel_ang < (M_PI * 0.45)) ? 1.3 : 0.8;
            double score_d = beam.distance_m * dir_bias;
            if (score_d > best_score_dist) {
                best_score_dist = score_d;
                best_beam_dist = beam.distance_m;
                best_beam_angle = beam.angle_rad;
            }
        }
        target_x = robot.x + std::cos(best_beam_angle) * std::max(0.5, best_beam_dist * 0.7);
        target_y = robot.y + std::sin(best_beam_angle) * std::max(0.5, best_beam_dist * 0.7);
    }

    // 2. DWA (Dynamic Window Approach) Velocity-Space Trajectory Evaluation
    // 5 forward linear velocity levels (strictly non-negative in normal mode; reverse is last resort)
    static const double v_ratios[5] = {0.0, 0.25, 0.50, 0.75, 1.00};
    // 11 angular velocity levels
    static const double omega_ratios[11] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);
    const double T_horiz = 0.60; // 0.60s prediction horizon
    const int N_sub = 4;
    const double dt_sub = T_horiz / N_sub;

    double best_score = -1e9;
    double best_v = 0.0;
    double best_omega = 0.0;
    bool found_valid = false;

    for (double v_ratio : v_ratios) {
        double cand_v = v_ratio * robot_config.v_max;

        for (double w_ratio : omega_ratios) {
            double cand_omega = w_ratio * robot_config.omega_max;

            // Trajectory forward rollout and collision check
            bool collision = false;
            double traj_theta = robot.theta;
            double traj_x = robot.x;
            double traj_y = robot.y;
            double min_traj_d = std::hypot(target_x - robot.x, target_y - robot.y);

            for (int step = 1; step <= N_sub; ++step) {
                double t = step * dt_sub;
                if (std::abs(cand_omega) < 1e-4) {
                    traj_theta = robot.theta;
                    traj_x = robot.x + cand_v * std::cos(robot.theta) * t;
                    traj_y = robot.y + cand_v * std::sin(robot.theta) * t;
                } else {
                    traj_theta = normalizeAngle(robot.theta + cand_omega * t);
                    double r_curve = cand_v / cand_omega;
                    traj_x = robot.x + r_curve * (std::sin(traj_theta) - std::sin(robot.theta));
                    traj_y = robot.y - r_curve * (std::cos(traj_theta) - std::cos(robot.theta));
                }

                double d_step = std::hypot(target_x - traj_x, target_y - traj_y);
                if (d_step < min_traj_d) min_traj_d = d_step;

                // Static obstacle collision check
                if (!checkRobotCollision(traj_x, traj_y, eff_r)) {
                    collision = true;
                    break;
                }

                // Inter-robot collision check
                for (const auto& other : robots) {
                    if (other.id == robot.id) continue;
                    double odist = std::hypot(traj_x - other.x, traj_y - other.y);
                    if (odist < (eff_r * 2.1)) {
                        collision = true;
                        break;
                    }
                }
                if (collision) break;
            }

            if (collision) continue;

            // Collision-free candidate: evaluate DWA objective score
            found_valid = true;

            // Heading alignment: evaluate alignment towards target
            double to_target_angle = std::atan2(target_y - robot.y, target_x - robot.x);
            double heading_diff = std::abs(normalizeAngle(to_target_angle - traj_theta));
            double score_heading = (M_PI - heading_diff) / M_PI; // [0, 1]

            // Progress towards target (reward closest approach along rollout)
            double cur_d = std::hypot(target_x - robot.x, target_y - robot.y);
            double score_progress = std::clamp((cur_d - min_traj_d) / (robot_config.v_max * T_horiz + 1e-4), -1.0, 1.0);

            // Obstacle clearance along trajectory endpoint from LiDAR
            double min_lidar_d = robot_config.lidar.range_m;
            for (const auto& beam : robot.scan) {
                double b_dist = std::hypot(beam.hit_x - traj_x, beam.hit_y - traj_y);
                if (b_dist < min_lidar_d) min_lidar_d = b_dist;
            }
            double score_clearance = std::clamp(min_lidar_d / (cell_size_m * 3.0), 0.0, 1.0);

            // Forward speed encouragement
            double score_speed = cand_v / robot_config.v_max;

            // Weighted composite score
            double score = 3.0 * score_heading + 2.5 * score_progress + 1.5 * score_clearance + 0.8 * score_speed;

            // Heavily penalize sitting completely idle (v=0, w=0) when moving or turning is available
            if (cand_v == 0.0 && std::abs(cand_omega) < 1e-3) {
                score -= 4.0;
            }
            // In-place rotation penalty if heading is already roughly aligned
            else if (cand_v == 0.0 && heading_diff < (M_PI * 0.35)) {
                score -= 2.0;
            }

            if (score > best_score) {
                best_score = score;
                best_v = cand_v;
                best_omega = cand_omega;
            }
        }
    }

    if (found_valid) {
        robot.v = best_v;
        robot.omega = best_omega;
    } else {
        // Last resort: all forward trajectories collide.
        // Check if in-place rotation is possible
        if (checkRobotCollision(robot.x, robot.y, eff_r)) {
            robot.v = 0.0;
            robot.omega = (robot.escape_turn_dir > 0.0 ? 1.0 : -1.0) * robot_config.omega_max;
        } else {
            // Truly wedged against wall: gentle reverse as absolute last resort
            robot.v = -0.25 * robot_config.v_max;
            robot.omega = 0.0;
        }
    }
}

/**
 * @brief Integrates continuous SE(2) unicycle kinematics with multi-stage collision resolution.
 *
 * Discrete Euler integration:
 *   next_theta = normalizeAngle(theta + omega * dt)
 *   dx = v * cos(theta) * dt
 *   dy = v * sin(theta) * dt
 *
 * Collision Resolution Hierarchy:
 * 1. Full Step: Evaluates circular footprint at (next_x, next_y). If clear, commits full step.
 * 2. X-Slide: If diagonal step collides, tests sliding along X axis only (next_x, y).
 * 3. Y-Slide: Tests sliding along Y axis only (x, next_y).
 * 4. Half-Step: Tests reduced half-displacement (x + 0.5*dx, y + 0.5*dy).
 * 5. Corner Repulsion: If completely blocked, rotates in place and applies a gentle repulsive nudge
 *    along the nearest obstacle corner normal to prevent wedge lock.
 */
void MultiRobotSimulator::integrateKinematics(RobotState& robot, double dt) {
    double old_x = robot.x;
    double old_y = robot.y;
    double old_theta = robot.theta;

    double next_theta = normalizeAngle(robot.theta + robot.omega * dt);
    double dx = robot.v * std::cos(robot.theta) * dt;
    double dy = robot.v * std::sin(robot.theta) * dt;
    double next_x = robot.x + dx;
    double next_y = robot.y + dy;

    double eff_r = std::min(robot_config.radius_m, cell_size_m * 0.38);

    // 1. Try full diagonal step (supports forward and reverse)
    if (checkRobotCollision(next_x, next_y, eff_r)) {
        robot.x = next_x;
        robot.y = next_y;
        robot.theta = next_theta;
    }
    // 2. Try sliding in X only
    else if (std::abs(dx) > 1e-6 && checkRobotCollision(next_x, robot.y, eff_r)) {
        robot.x = next_x;
        robot.theta = next_theta;
    }
    // 3. Try sliding in Y only
    else if (std::abs(dy) > 1e-6 && checkRobotCollision(robot.x, next_y, eff_r)) {
        robot.y = next_y;
        robot.theta = next_theta;
    }
    // 4. Try half-step
    else if (checkRobotCollision(robot.x + 0.5 * dx, robot.y + 0.5 * dy, eff_r)) {
        robot.x += 0.5 * dx;
        robot.y += 0.5 * dy;
        robot.theta = next_theta;
    }
    // 5. Blocked: rotate in place, and if touching an obstacle corner, gently push away along vertex normal
    else {
        robot.theta = next_theta;
        double push_dx = 0.0;
        double push_dy = 0.0;
        int min_c = static_cast<int>(std::floor((robot.x - eff_r * 1.5) / cell_size_m));
        int max_c = static_cast<int>(std::floor((robot.x + eff_r * 1.5) / cell_size_m));
        int min_r = static_cast<int>(std::floor((robot.y - eff_r * 1.5) / cell_size_m));
        int max_r = static_cast<int>(std::floor((robot.y + eff_r * 1.5) / cell_size_m));
        for (int r_idx = min_r; r_idx <= max_r; ++r_idx) {
            for (int c_idx = min_c; c_idx <= max_c; ++c_idx) {
                if (ground_truth.inBounds(r_idx, c_idx) && ground_truth.getObserved(r_idx, c_idx) == CellState::Obstacle) {
                    double cell_min_x = c_idx * cell_size_m;
                    double cell_max_x = (c_idx + 1) * cell_size_m;
                    double cell_min_y = r_idx * cell_size_m;
                    double cell_max_y = (r_idx + 1) * cell_size_m;

                    double closest_x = std::clamp(robot.x, cell_min_x, cell_max_x);
                    double closest_y = std::clamp(robot.y, cell_min_y, cell_max_y);
                    double ox = robot.x - closest_x;
                    double oy = robot.y - closest_y;
                    double odist = std::hypot(ox, oy);
                    if (odist < (eff_r * 1.35) && odist > 1e-5) {
                        push_dx += (ox / odist);
                        push_dy += (oy / odist);
                    }
                }
            }
        }
        if (std::hypot(push_dx, push_dy) > 1e-3) {
            double push_mag = std::hypot(push_dx, push_dy);
            double ux = push_dx / push_mag;
            double uy = push_dy / push_mag;
            for (double step_scale : {0.15, 0.08, 0.04}) {
                double nx = robot.x + ux * (cell_size_m * step_scale);
                double ny = robot.y + uy * (cell_size_m * step_scale);
                if (checkRobotCollision(nx, ny, eff_r)) {
                    robot.x = nx;
                    robot.y = ny;
                    break;
                }
            }
        }
    }

    // Stuck detection: position hasn't made net translation progress despite commanded translation
    double dist_moved = std::hypot(robot.x - old_x, robot.y - old_y);
    double angle_turned = std::abs(normalizeAngle(robot.theta - old_theta));
    bool is_rotating_in_place = (std::abs(robot.omega) > 0.2) || (angle_turned > 0.015);
    bool commanded_translation = (std::abs(robot.v) > 0.04) || (robot.has_goal && !is_rotating_in_place);

    if (commanded_translation && dist_moved < (0.15 * robot_config.v_max * dt)) {
        robot.stuck_ticks++;
        // If stuck for >= 14 ticks (~0.70s), trigger Nav2 BackUp or RepulsiveSpin recovery maneuver!
        if (robot.stuck_ticks >= 14 && robot.recovery_phase != RobotState::RecoveryPhase::BackUp &&
            robot.recovery_phase != RobotState::RecoveryPhase::RepulsiveSpin) {
            // Check rear clearance from 360 LiDAR before backing up
            double rear_clearance = robot_config.lidar.range_m;
            for (const auto& beam : robot.scan) {
                double rel_ang = normalizeAngle(beam.angle_rad - robot.theta);
                if (std::abs(rel_ang) > (M_PI * 0.55)) {
                    double bwd_d = beam.distance_m * std::abs(std::cos(rel_ang));
                    double lat_d = beam.distance_m * std::abs(std::sin(rel_ang));
                    if (lat_d < (eff_r * 1.25) && bwd_d > 0.0) {
                        rear_clearance = std::min(rear_clearance, bwd_d);
                    }
                }
            }
            if (rear_clearance > (cell_size_m * 1.0)) {
                robot.recovery_phase = RobotState::RecoveryPhase::BackUp;
                robot.recovery_ticks = 14; // ~0.70s of backing up
            } else {
                robot.recovery_phase = RobotState::RecoveryPhase::RepulsiveSpin;
                robot.recovery_ticks = 16; // Direct spin away from obstacle normal if rear is blocked
            }

            robot.waypoint_stuck_count++;
            if (robot.waypoint_stuck_count > 3 && robot.has_goal) {
                robot.tour.advance(); // Skip this repeatedly stuck goal!
                if (robot.tour.empty()) {
                    robot.has_goal = false;
                    robot.assigned_frontier = {-1, -1};
                } else {
                    robot.assigned_frontier = robot.tour.currentGoal();
                }
                robot.next_waypoint = {-1, -1};
                robot.waypoint_stuck_count = 0;
            }
        }
    } else if (dist_moved >= (0.15 * robot_config.v_max * dt)) {
        if (robot.stuck_ticks > 0) robot.stuck_ticks--;
        if (robot.recovery_phase == RobotState::RecoveryPhase::ReplanNudge) {
            robot.recovery_phase = RobotState::RecoveryPhase::None;
        }
    } else if (is_rotating_in_place) {
        // Deliberate in-place rotation: making active rotational progress, NOT stuck!
        if (robot.stuck_ticks > 0) robot.stuck_ticks--;
    } else {
        if (robot.stuck_ticks > 0) robot.stuck_ticks--;
    }

    // Net translation tracking for auto-respawn policy (50 steps without significant translation)
    double dist_from_stuck_ref = std::hypot(robot.x - robot.stuck_ref_x, robot.y - robot.stuck_ref_y);
    double threshold_m = (stuck_dist_threshold_m > 0.0) ? stuck_dist_threshold_m : std::max(0.20, cell_size_m * 1.5);
    if (dist_from_stuck_ref < threshold_m) {
        robot.stuck_translation_ticks++;
        if (robot.stuck_translation_ticks >= 30 && robot.has_goal && robot.recovery_phase == RobotState::RecoveryPhase::None) {
            robot.stuck_ticks = 14;
        }
    } else {
        robot.stuck_ref_x = robot.x;
        robot.stuck_ref_y = robot.y;
        robot.stuck_translation_ticks = 0;
    }

    // Record trail point every few steps
    if (total_ticks % 3 == 0) {
        robot.trail.push_back({static_cast<float>(robot.x), static_cast<float>(robot.y)});
        if (robot.trail.size() > 250) {
            robot.trail.erase(robot.trail.begin(), robot.trail.begin() + 25);
        }
    }
}

/**
 * @brief Advances the complete continuous multi-robot exploration simulation by dt seconds.
 *
 * Full step sequence:
 * 1. Commits cells indexed by background worker into indexed_map.
 * 2. Measures simulation frequency and updates tile visual decay timers.
 * 3. Detects frontiers on observed (or indexed) map and groups into spatial tile buckets.
 * 4. Dispatches m-TSP fleet assignment and all-pairs shortest-path queries.
 * 5. Evaluates autonomous controls (carrot following, obstacle repulsion, recovery).
 * 6. Integrates SE(2) kinematics and checks for stationary/stuck robots.
 * 7. Sweeps planar LiDAR rays to reveal fog-of-war.
 * 8. Flushes newly opened free cells into H-TRIP hierarchy (sync or async).
 * 9. Streams telemetry traces to mission logger if active.
 */
void MultiRobotSimulator::step(double dt) {
    total_ticks++;

    // Commit any cells whose H-TRIP tree indexing completed on background worker
    {
        std::lock_guard<std::mutex> lock(commit_mutex);
        if (!completed_indexed_cells.empty()) {
            for (const auto& pt : completed_indexed_cells) {
                indexed_map.setObserved(pt, CellState::Free);
            }
            completed_indexed_cells.clear();
        }
    }

    // Track real-time simulation step frequency (Hz) independently of map updates
    auto now = std::chrono::high_resolution_clock::now();
    if (last_step_time.time_since_epoch().count() > 0) {
        double dt_actual_s = std::chrono::duration<double>(now - last_step_time).count();
        if (dt_actual_s > 0.0001 && dt_actual_s < 2.0) {
            double inst_sim_hz = 1.0 / dt_actual_s;
            sim_step_freq_hz = (sim_step_freq_hz <= 0.0)
                ? inst_sim_hz
                : (0.90 * sim_step_freq_hz + 0.10 * inst_sim_hz);
        }
    }
    last_step_time = now;

    // Decay highlight timers for recently updated tiles
    {
        std::lock_guard<std::mutex> lock(decay_mutex);
        for (auto& decay : active_tile_decay) {
            if (decay > 0) decay--;
        }
        for (auto& decay : active_index_decay) {
            if (decay > 0) decay--;
        }
    }
    if (total_ticks - last_update_tick > 40) {
        map_update_freq_hz *= 0.96;
    }

    for (auto& bot : robots) {
        if (bot.persistence_ticks > 0) {
            bot.persistence_ticks--;
        }
    }

    // 1. Detect active frontiers / patrol targets
    if (!enable_fog_of_war) {
        for (auto it = patrol_targets.begin(); it != patrol_targets.end(); ) {
            bool reached = false;
            for (const auto& bot : robots) {
                Point bcell = worldToCell(bot.x, bot.y);
                double d = std::hypot(bcell.r - it->r, bcell.c - it->c);
                if (d <= 6.0) { // within ~1.2m
                    reached = true;
                    break;
                }
            }
            if (reached) {
                targets_visited_count++;
                it = patrol_targets.erase(it);
            } else {
                ++it;
            }
        }
        replenishPatrolTargets(static_cast<size_t>(max_frontiers_eval > 0 ? max_frontiers_eval : 25));

        frontiers.clear();
        for (size_t i = 0; i < patrol_targets.size(); ++i) {
            FrontierCluster cl;
            cl.id = static_cast<int>(i + 1);
            cl.representative = patrol_targets[i];
            cl.cells = {patrol_targets[i]};
            frontiers.push_back(std::move(cl));
        }
        viewpoints = patrol_targets;
    } else if (enable_live_frontier_projection) {
        // Mode B: Detect frontiers on the live observed grid as LiDAR sweeps in real time.
        // Newly opened frontier viewpoints are projected (snapped) to the nearest indexed cell
        // during H-TRIP queries so exploration planning never lags behind perception.
        frontiers = frontier_manager.detectAndCluster(observed);
    } else {
        // Mode A: Strictly detect frontiers on the double-buffered indexed map.
        // Frontiers only emerge after background worker finishes H-TRIP indexing and commits cells.
        frontiers = frontier_manager.detectAndCluster(indexed_map);
    }

    if (enable_fog_of_war) {
        viewpoints.clear();
        int max_vp = std::min(static_cast<int>(frontiers.size()), max_frontiers_eval);
        viewpoints.reserve(static_cast<size_t>(max_vp));

        if (static_cast<int>(frontiers.size()) <= max_vp) {
            for (const auto& f : frontiers) {
                viewpoints.push_back(f.representative);
            }
        } else {
            // Tile Binning: Homogenize frontier selection across H-TRIP spatial tiles
            int tile_b = (oracle.config.tile_size_b > 0) ? oracle.config.tile_size_b : 16;
            int num_tile_cols = (observed.cols + tile_b - 1) / tile_b;
            int num_tile_rows = (observed.rows + tile_b - 1) / tile_b;
            int total_tiles = num_tile_rows * num_tile_cols;

            std::vector<std::vector<size_t>> tile_buckets(static_cast<size_t>(total_tiles));
            for (size_t i = 0; i < frontiers.size(); ++i) {
                const auto& p = frontiers[i].representative;
                int tr = std::clamp(p.r / tile_b, 0, num_tile_rows - 1);
                int tc = std::clamp(p.c / tile_b, 0, num_tile_cols - 1);
                int tid = tr * num_tile_cols + tc;
                tile_buckets[static_cast<size_t>(tid)].push_back(i);
            }

            std::vector<size_t> active_tiles;
            active_tiles.reserve(static_cast<size_t>(total_tiles));
            for (size_t t = 0; t < tile_buckets.size(); ++t) {
                if (!tile_buckets[t].empty()) {
                    active_tiles.push_back(t);
                    // Sort clusters in descending order of information gain (cells count)
                    std::sort(tile_buckets[t].begin(), tile_buckets[t].end(),
                              [&](size_t a, size_t b) {
                                  return frontiers[a].cells.size() > frontiers[b].cells.size();
                              });
                }
            }

            if (!active_tiles.empty()) {
                size_t round = 0;
                size_t start_idx = 0;
                while (viewpoints.size() < static_cast<size_t>(max_vp)) {
                    bool added_in_round = false;
                    for (size_t i = 0; i < active_tiles.size(); ++i) {
                        size_t t_idx = active_tiles[(start_idx + i) % active_tiles.size()];
                        auto& bucket = tile_buckets[t_idx];
                        if (round < bucket.size()) {
                            viewpoints.push_back(frontiers[bucket[round]].representative);
                            added_in_round = true;
                            if (viewpoints.size() == static_cast<size_t>(max_vp)) break;
                        }
                    }
                    if (!added_in_round) break;
                    round++;
                }
            }
        }
    }

    // 2. Perform intelligent multi-robot m-TSP frontier tour allocation
    assignFrontiersToFleet();

    // 3. Update robot controls & continuous kinematics
    for (auto& bot : robots) {
        computeAutonomousControls(bot, dt);
        integrateKinematics(bot, dt);
    }

    // Auto-respawn robots that have remained in place (translation < threshold for >= 50 steps)
    checkAndRespawnStuckRobots();

    // 4. Cast LiDAR scans and reveal fog-of-war
    for (auto& bot : robots) {
        castLidarRays(bot);
    }

    // 5. Newly revealed cells are flushed incrementally into H-TRIP
    if (!newly_opened_cells.empty() || !async_pending_cells.empty()) {
        flushPendingUpdates();
    } else {
        if (!is_updating_htrip.load()) {
            last_htrip_update_us = 0.0;
            last_revealed_cells_count = 0;
            updated_leaf_indices.clear();
        }
    }

    if (is_updating_htrip.load()) {
        replica_lag_ticks.store(total_ticks > last_update_tick ? total_ticks - last_update_tick : 0);
    } else {
        replica_lag_ticks.store(0);
    }
    pending_openings.store(async_pending_cells.size());

    // 6. Record high-resolution academic telemetry if active
    if (mission_logger.isLogging()) {
        mission_logger.logStep(*this);
    }
}

void MultiRobotSimulator::runSyncBaselines(const GridMap& map_snapshot, const std::vector<Point>& pts, double htrip_us, uint64_t tick) {
    if (pts.size() < 2) return;

    // 1. Full exhaustive all-pairs BFS
    double bfs_us = 0.0;
    if (enable_bfs_baseline) {
        auto t_bfs0 = std::chrono::high_resolution_clock::now();
        auto bfs_mat = BfsOracle::queryAllPairs(map_snapshot, pts);
        auto t_bfs1 = std::chrono::high_resolution_clock::now();
        bfs_us = std::chrono::duration<double, std::micro>(t_bfs1 - t_bfs0).count();
        (void)bfs_mat;
    }

    // 2. Representative A* sampling to measure actual per-pair search time and extrapolate to all pairs
    const size_t k_pts = pts.size();
    const size_t total_pairs = k_pts * (k_pts - 1) / 2;
    double astar_us = 0.0;

    if (enable_astar_baseline && total_pairs > 0) {
        const size_t sample_target = std::min(size_t{16}, total_pairs);
        auto t_astar0 = std::chrono::high_resolution_clock::now();
        size_t evaluated = 0;
        size_t stride = std::max(size_t{1}, total_pairs / sample_target);
        size_t pair_idx = 0;

        for (size_t i = 0; i < k_pts && evaluated < sample_target; ++i) {
            for (size_t j = i + 1; j < k_pts && evaluated < sample_target; ++j) {
                if ((pair_idx % stride) == 0) {
                    (void)AStarOracle::queryDistance(map_snapshot, pts[i], pts[j]);
                    evaluated++;
                }
                pair_idx++;
            }
        }
        auto t_astar1 = std::chrono::high_resolution_clock::now();
        double sample_time_us = std::chrono::duration<double, std::micro>(t_astar1 - t_astar0).count();
        if (evaluated > 0) {
            astar_us = (sample_time_us / static_cast<double>(evaluated)) * static_cast<double>(total_pairs);
        }
    }

    uint64_t job_id = bfs_job_seq.fetch_add(1, std::memory_order_relaxed) + 1;
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        commitBfsCompareLocked(job_id, tick, pts.size(), htrip_us, bfs_us, astar_us);
    }
    if (mission_logger.isLogging()) {
        mission_logger.logBfsCompare(job_id, tick, pts.size(), htrip_us, bfs_us, astar_us);
    }
}

void MultiRobotSimulator::runBfsBaseline() {
    if (last_query_pts.size() < 2) return;
    GridMap map_copy;
    {
        std::lock_guard<std::mutex> lock(query_oracle_mutex);
        map_copy = *active_query_map_ptr.load(std::memory_order_acquire);
    }
    runBfsBaseline(std::move(map_copy), last_query_pts);
}

void MultiRobotSimulator::runBfsBaseline(GridMap map_snapshot, std::vector<Point> pts) {
    if (pts.size() < 2) return;
    createBfsCompareJob(std::move(map_snapshot), std::move(pts));
}

uint64_t MultiRobotSimulator::createBfsCompareJob(GridMap map_snapshot, std::vector<Point> pts) {
    const uint64_t job_id = bfs_job_seq.fetch_add(1, std::memory_order_relaxed) + 1;
    BfsCompareJob job;
    job.job_id = job_id;
    job.tick = total_ticks;
    job.num_points = pts.size();
    job.map = std::move(map_snapshot);
    job.pts = std::move(pts);
    job.occupied = true;

    bool launch_now = false;
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        if (is_running_bfs.load()) {
            bfs_pending_job = std::move(job);
        } else {
            is_running_bfs.store(true);
            launch_now = true;
        }
    }
    if (launch_now) {
        launchBfsWorker(std::move(job.map), std::move(job.pts), job.job_id, job.tick,
                        job.num_points, job.htrip_query_us, job.has_htrip);
    }
    return job_id;
}

void MultiRobotSimulator::launchBfsWorker(GridMap map_snapshot, std::vector<Point> pts, uint64_t job_id,
                                         uint64_t tick, size_t num_points, double htrip_query_us, bool has_htrip) {
    if (bfs_worker_thread.joinable()) {
        bfs_worker_thread.join();
    }
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        bfs_inflight_job = {};
        bfs_inflight_job.occupied = true;
        bfs_inflight_job.job_id = job_id;
        bfs_inflight_job.tick = tick;
        bfs_inflight_job.num_points = num_points;
        bfs_inflight_job.htrip_query_us = htrip_query_us;
        bfs_inflight_job.has_htrip = has_htrip;
    }
    bfs_worker_thread = std::thread([this, job_id, obs = std::move(map_snapshot), q_pts = std::move(pts)]() {
        double bfs_us = 0.0;
        if (this->enable_bfs_baseline) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto bfs_matrix = BfsOracle::queryAllPairs(obs, q_pts);
            auto t1 = std::chrono::high_resolution_clock::now();
            bfs_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
            (void)bfs_matrix;
        }

        double astar_us = 0.0;
        if (this->enable_astar_baseline) {
            auto t2 = std::chrono::high_resolution_clock::now();
            auto astar_matrix = AStarOracle::queryAllPairs(obs, q_pts);
            auto t3 = std::chrono::high_resolution_clock::now();
            astar_us = std::chrono::duration<double, std::micro>(t3 - t2).count();
            (void)astar_matrix;
        }

        this->onBfsWorkerFinished(job_id, bfs_us, astar_us);
    });
}

void MultiRobotSimulator::pumpBfsPending() {
    BfsCompareJob job;
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        if (is_running_bfs.load() || !bfs_pending_job.occupied) return;
        is_running_bfs.store(true);
        job = std::move(bfs_pending_job);
        bfs_pending_job = {};
    }
    launchBfsWorker(std::move(job.map), std::move(job.pts), job.job_id, job.tick,
                    job.num_points, job.htrip_query_us, job.has_htrip);
}

void MultiRobotSimulator::noteHtripQueryTime(uint64_t job_id, double htrip_us) {
    if (job_id == 0 || htrip_us <= 0.0) return;
    bool emit = false;
    BfsCompareSample sample;
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        if (bfs_inflight_job.occupied && bfs_inflight_job.job_id == job_id) {
            bfs_inflight_job.htrip_query_us = htrip_us;
            bfs_inflight_job.has_htrip = true;
            if (bfs_inflight_job.has_bfs) {
                commitBfsCompareLocked(job_id, bfs_inflight_job.tick, bfs_inflight_job.num_points,
                                       htrip_us, bfs_inflight_job.bfs_query_us, bfs_inflight_job.astar_query_us);
                bfs_inflight_job.occupied = false;
                emit = true;
                sample = bfs_compare_samples.back();
            }
        } else if (bfs_pending_job.occupied && bfs_pending_job.job_id == job_id) {
            bfs_pending_job.htrip_query_us = htrip_us;
            bfs_pending_job.has_htrip = true;
        } else if (bfs_await_htrip_job.occupied && bfs_await_htrip_job.job_id == job_id) {
            commitBfsCompareLocked(job_id, bfs_await_htrip_job.tick, bfs_await_htrip_job.num_points,
                                   htrip_us, bfs_await_htrip_job.bfs_query_us, bfs_await_htrip_job.astar_query_us);
            bfs_await_htrip_job.occupied = false;
            emit = true;
            sample = bfs_compare_samples.back();
        }
    }
    if (emit && mission_logger.isLogging()) {
        mission_logger.logBfsCompare(sample.job_id, sample.tick, sample.num_points,
                                     sample.htrip_query_us, sample.bfs_query_us, sample.astar_query_us);
    }
}

void MultiRobotSimulator::onBfsWorkerFinished(uint64_t job_id, double bfs_us, double astar_us) {
    bool emit = false;
    BfsCompareSample sample;
    {
        std::lock_guard<std::mutex> lock(bfs_mutex);
        if (bfs_inflight_job.occupied && bfs_inflight_job.job_id == job_id) {
            bfs_inflight_job.bfs_query_us = bfs_us;
            bfs_inflight_job.astar_query_us = astar_us;
            bfs_inflight_job.has_bfs = true;
            if (bfs_inflight_job.has_htrip) {
                commitBfsCompareLocked(job_id, bfs_inflight_job.tick, bfs_inflight_job.num_points,
                                       bfs_inflight_job.htrip_query_us, bfs_us, astar_us);
                bfs_inflight_job.occupied = false;
                emit = true;
                sample = bfs_compare_samples.back();
            } else {
                bfs_await_htrip_job = bfs_inflight_job;
                bfs_await_htrip_job.map = {};
                bfs_await_htrip_job.pts.clear();
                bfs_inflight_job.occupied = false;
            }
        }
        is_running_bfs.store(false);
    }
    if (emit && mission_logger.isLogging()) {
        mission_logger.logBfsCompare(sample.job_id, sample.tick, sample.num_points,
                                     sample.htrip_query_us, sample.bfs_query_us, sample.astar_query_us);
    }
}

void MultiRobotSimulator::commitBfsCompareLocked(uint64_t job_id, uint64_t tick, size_t num_points,
                                                double htrip_us, double bfs_us, double astar_us) {
    BfsCompareSample sample;
    sample.job_id = job_id;
    sample.tick = tick;
    sample.num_points = num_points;
    sample.htrip_query_us = htrip_us;
    sample.bfs_query_us = bfs_us;
    sample.astar_query_us = astar_us;
    sample.speedup = (htrip_us > 0.0) ? (bfs_us / htrip_us) : 0.0;
    sample.astar_speedup = (htrip_us > 0.0) ? (astar_us / htrip_us) : 0.0;
    bfs_compare_samples.push_back(sample);
    last_bfs_query_us.store(bfs_us);
    last_astar_query_us.store(astar_us);
    last_paired_htrip_query_us.store(htrip_us);
    last_paired_speedup.store(sample.speedup);
    last_paired_astar_speedup.store(sample.astar_speedup);
    bfs_compare_count.store(bfs_compare_samples.size());
}

void MultiRobotSimulator::drainBfsComparisons() {
    if (!enable_bfs_baseline && !enable_astar_baseline) return;

    while (is_querying_pairwise.load()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    if (query_worker_thread.joinable()) {
        query_worker_thread.join();
    }

    for (;;) {
        pumpBfsPending();
        if (is_running_bfs.load()) {
            if (bfs_worker_thread.joinable()) {
                bfs_worker_thread.join();
            }
            continue;
        }
        bool pending = false;
        {
            std::lock_guard<std::mutex> lock(bfs_mutex);
            pending = bfs_pending_job.occupied;
        }
        if (!pending) break;
    }
}

} // namespace htrip
