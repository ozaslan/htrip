#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include "minplus_oracle.hpp"
#include "frontier.hpp"
#include "tour_planner.hpp"
#include "mission_logger.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <thread>
#include <memory>

namespace htrip {

/**
 * @brief Configuration parameters for a 2D planar range-finding LiDAR sensor.
 *
 * Models a planar scanning LiDAR emitting rays uniformly spaced across an angular sector
 * or a full 360-degree sweep. In the simulation loop, raymarching is executed over the
 * continuous metric domain and converted to grid coordinates to simulate realistic line-of-sight
 * visibility, obstacle detection, and fog-of-war revelation.
 */
struct LidarConfig {
    double range_m = 8.0;      ///< Maximum sensing horizon in meters; rays terminate here if unobstructed.
    double span_deg = 360.0;   ///< Total angular field-of-view (FOV) in degrees (e.g., 360 for omnidirectional, 180 or 90 for forward sector).
    int num_rays = 360;        ///< Total number of discrete angular beams evaluated in each scan.

    /**
     * @brief Computes the field-of-view span in radians.
     */
    [[nodiscard]] double spanRad() const noexcept {
        return span_deg * (M_PI / 180.0);
    }

    /**
     * @brief Computes the angular increment between adjacent beams in degrees.
     *
     * For a full 360-degree sweep, the increment is span / num_rays (periodic wrap).
     * For a restricted sector (span < 360), the increment is span / (num_rays - 1) to cover both boundary beams.
     */
    [[nodiscard]] double angularResolutionDeg() const noexcept {
        if (num_rays <= 1) return span_deg;
        return span_deg / static_cast<double>(span_deg >= 360.0 ? num_rays : (num_rays - 1));
    }

    /**
     * @brief Computes the angular increment between adjacent beams in radians.
     */
    [[nodiscard]] double angularResolutionRad() const noexcept {
        return angularResolutionDeg() * (M_PI / 180.0);
    }
};

/**
 * @brief Physical, kinematic, and sensor configuration for a differential-drive / unicycle robot.
 *
 * Enforces non-holonomic motion constraints and physical footprint dimensions for collision checking.
 */
struct RobotConfig {
    double radius_m = 0.20;    ///< Robot circular collision footprint radius in meters.
    double v_max = 0.60;       ///< Maximum allowable forward linear velocity (m/s).
    double omega_max = 1.80;   ///< Maximum allowable angular rotational velocity (rad/s).
    LidarConfig lidar;         ///< Mounted planar LiDAR configuration.
};

/**
 * @brief Measurement payload resulting from a single simulated LiDAR beam.
 *
 * Stores both polar and world Cartesian endpoints, along with obstacle termination status.
 */
struct LidarHit {
    double angle_rad = 0.0;    ///< Global polar heading of the ray in world coordinates.
    double distance_m = 0.0;   ///< Measured radial distance to first obstacle or max sensor range (meters).
    double hit_x = 0.0;        ///< Cartesian X coordinate of the beam endpoint in the world frame.
    double hit_y = 0.0;        ///< Cartesian Y coordinate of the beam endpoint in the world frame.
    bool hit_obstacle = false; ///< True if the ray terminated upon striking an obstacle cell; false if max range was reached.
};

/**
 * @brief Operational execution policy for baseline metric comparisons (BFS and A*).
 *
 * Enables fair runtime benchmarking against classical graph search without degrading real-time simulation interactivity.
 */
enum class BaselineExecutionMode {
    Asynchronous, ///< Interactive exploration mode. Baselines execute concurrently on background worker threads without blocking simulation ticks.
    Synchronous   ///< Exact benchmarking mode. Every query tick stalls until H-TRIP, BFS, and A* complete on the exact same snapshot.
};

/**
 * @brief Continuous SE(2) pose, kinematic state, and navigation controller context for an individual robot.
 *
 * Encapsulates the complete runtime state of a robot, including:
 * 1. Physical kinematics: planar position (x, y), heading theta, forward speed v, and angular rate omega.
 * 2. Visual telemetry: position trail and latest LiDAR scan hits for GUI/rendering.
 * 3. Exploration hierarchy: assigned m-TSP tour, current target frontier, and planned corridor path.
 * 4. Recovery state machine: multi-stage un-stucking behavior (ReplanNudge, BackUp, RepulsiveSpin).
 */
struct RobotState {
    int id = 0;
    double x = 0.0;          ///< World X position in continuous Cartesian space (meters).
    double y = 0.0;          ///< World Y position in continuous Cartesian space (meters).
    double theta = 0.0;      ///< Planar heading angle in radians, normalized to [-pi, +pi].
    double v = 0.0;          ///< Current forward linear velocity (m/s).
    double omega = 0.0;      ///< Current angular rotational velocity (rad/s).

    uint32_t color = 0xFF3388EE; ///< RGBA color identifier for rendering.
    std::vector<std::pair<float, float>> trail; ///< Historic trajectory positions for visualization trails.
    std::vector<LidarHit> scan; ///< Most recent LiDAR scan hits.

    // Frontier assignment and corridor navigation state
    RobotTour tour;                  ///< Pre-planned m-TSP exploration tour assigned to this robot.
    Point assigned_frontier{-1, -1}; ///< Immediate target frontier viewpoint (grid coordinates, matching tour.currentGoal()).
    int assigned_frontier_id = -1;   ///< Cluster ID of the assigned frontier (-1 if unassigned).
    bool has_goal = false;           ///< True if the robot is actively navigating toward an assigned frontier.
    Point next_waypoint{-1, -1};     ///< Immediate next intermediate grid cell waypoint along the planned path.
    std::vector<Point> planned_path; ///< Obstacle-cleared A* grid path connecting robot pose to assigned frontier.
    size_t path_index = 0;           ///< Index of the currently pursued waypoint within planned_path.

    /**
     * @brief Discrete phases of the autonomous recovery state machine.
     *
     * Mimics professional mobile robot recovery behaviors (e.g., ROS 2 Nav2 Recovery Server):
     * - ReplanNudge: Force immediate A* replanning from the current pose to break local deadlocks.
     * - BackUp: Execute negative linear velocity to disengage from wall or corner contact.
     * - RepulsiveSpin: Execute in-place rotation away from the obstacle gradient to find open corridor headings.
     */
    enum class RecoveryPhase : uint8_t {
        None = 0,
        ReplanNudge = 1,   ///< Re-plan path immediately from current pose.
        BackUp = 2,        ///< Reverse linear motion to clear obstacle contact.
        RepulsiveSpin = 3  ///< Pure in-place rotation away from obstacle normal.
    };
    RecoveryPhase recovery_phase = RecoveryPhase::None;
    int recovery_ticks = 0;          ///< Remaining simulation ticks allocated to the active recovery phase.

    int stuck_ticks = 0;             ///< Consecutive ticks where forward progress has been impeded.
    double stuck_ref_x = 0.0;        ///< Reference X position for measuring net translational displacement.
    double stuck_ref_y = 0.0;        ///< Reference Y position for measuring net translational displacement.
    int stuck_translation_ticks = 0; ///< Consecutive ticks where net displacement remained below threshold.
    int total_respawns = 0;          ///< Cumulative number of emergency respawn operations triggered for this robot.
    int persistence_ticks = 0;       ///< Remaining ticks of persistent commitment to the current assigned frontier.
    int escape_ticks = 0;            ///< Ticks dedicated to collision escape maneuver.
    int backtrack_ticks = 0;         ///< Ticks dedicated to backtracking toward open space.
    double backtrack_target_x = 0.0; ///< World X coordinate of backtrack destination.
    double backtrack_target_y = 0.0; ///< World Y coordinate of backtrack destination.
    int waypoint_stuck_count = 0;    ///< Counter of failed attempts to reach the immediate waypoint.
    double escape_turn_dir = 1.0;    ///< Direction of rotational escape (+1.0 counter-clockwise, -1.0 clockwise).
    double wander_angle_offset = 0.0;///< Offset applied during stochastic corridor wandering.
};

/**
 * @brief Continuous SE(2) multi-robot autonomous exploration simulator and H-TRIP benchmarking harness.
 *
 * MultiRobotSimulator coordinates the full exploration lifecycle across an arbitrary fleet of mobile robots:
 * 1. Physical Grid & Sensor Simulation:
 *    Maintains ground-truth and observed grid maps, simulates uniform polar LiDAR raycasting,
 *    and reveals unknown cells in real time (fog-of-war exploration).
 * 2. Double-Buffered Dynamic Hierarchical Indexing:
 *    Propagates newly revealed free cells into the H-TRIP hierarchy (`MinPlusOracle`).
 *    To maintain smooth, real-time simulation rates (e.g. 20 Hz) on large maps, dynamic index updates
 *    and pairwise all-pairs distance queries can run asynchronously on dedicated worker threads
 *    (`update_worker_thread`, `query_worker_thread`) using double-buffered grid and oracle instances.
 * 3. Frontier Extraction & Stand-off Viewpoint Placement:
 *    Detects frontier boundaries between explored free cells and unknown space, partitions them into
 *    connected clusters, and computes clearance-safe stand-off viewpoints (clearance >= 2).
 * 4. All-Pairs Shortest Path Queries & Fleet Tour Planning:
 *    Executes all-pairs distance queries across all robot positions and candidate frontier viewpoints
 *    using SIMD-accelerated H-TRIP or classical baselines (k x BFS, A*). Solves the multi-robot traveling
 *    salesperson problem (m-TSP) using greedy allocation and 2-opt local search.
 * 5. Autonomous Kinematic Control & Collision Avoidance:
 *    Drives unicycle robots along planned paths using obstacle-repulsive vector fields and a multi-tier
 *    recovery state machine (Replan, BackUp, RepulsiveSpin, and auto-respawn).
 */
class MultiRobotSimulator {
public:
    GridMap ground_truth;    ///< Complete ground-truth world representation loaded from MovingAI benchmark.
    GridMap observed;        ///< Live observed occupancy grid updated in real time by LiDAR raycasting.
    GridMap indexed_map;     ///< Double-buffered grid strictly matching the state indexed in H-TRIP.
    MinPlusOracle oracle;    ///< Primary H-TRIP dynamic reachability index on indexed_map.
    FrontierManager frontier_manager{2, 16}; ///< Frontier detection and stand-off viewpoint generation engine.

    std::vector<FrontierCluster> frontiers;
    std::vector<Point> viewpoints;
    std::vector<Point> last_query_pts;
    std::vector<dist_t> all_pairs_matrix;

    std::vector<RobotState> robots;
    RobotConfig robot_config;

    int scale_factor = 1;         ///< 1x, 2x, or 4x resolution scaling
    double cell_size_m = 0.05;    ///< Metric size of a single cell in the scaled grid
    double world_w_m = 0.0;       ///< Total physical world width in meters
    double world_h_m = 0.0;       ///< Total physical world height in meters

    int total_free_cells_gt = 0;
    int explored_free_cells = 0;

    std::vector<Point> newly_opened_cells; ///< Free cells revealed in the current tick
    int last_revealed_cells_count = 0;     ///< Count of cells processed in the latest update batch
    std::vector<int> updated_leaf_indices; ///< Leaf tile indices modified in the current tick
    std::vector<uint8_t> active_tile_decay;  ///< Decay ticks for base tiles being updated (warm amber glow)
    std::vector<uint8_t> active_index_decay; ///< Decay ticks for tiles updated across H-TRIP index (electric emerald glow)
    std::vector<uint32_t> dirty_leaf_stamp;  ///< Epoch stamps for dirty leaf tile deduplication (zero allocation)
    uint32_t dirty_stamp_epoch = 0;

    uint64_t total_map_updates = 0;         ///< Cumulative count of H-TRIP dynamic index batch updates
    uint64_t last_update_tick = 0;          ///< Simulation tick when the last batch update was triggered
    std::atomic<uint64_t> replica_lag_ticks{0}; ///< Ticks the query replica trails live occupancy while an update is in flight.
    std::atomic<uint64_t> pending_openings{0};  ///< Newly free cells waiting because the update worker is busy.
    double map_update_freq_hz = 0.0;        ///< Measured map update rate in Hz
    double sim_step_freq_hz = 0.0;          ///< Measured simulation step loop rate in Hz
    std::chrono::high_resolution_clock::time_point last_step_time{}; ///< Timestamp of previous step() call

    // Timing and performance metrics (microseconds, thread-safe atomics)
    std::atomic<double> last_htrip_update_us{0.0};
    std::atomic<double> last_htrip_query_us{0.0};
    std::atomic<double> last_mtsp_plan_us{0.0};
    std::atomic<double> last_bfs_query_us{0.0};
    std::atomic<double> last_astar_query_us{0.0};
    std::atomic<double> last_paired_htrip_query_us{0.0}; ///< H-TRIP time from the same snapshot as last baseline
    std::atomic<double> last_paired_speedup{0.0};        ///< last_bfs_query_us / last_paired_htrip_query_us
    std::atomic<double> last_paired_astar_speedup{0.0};  ///< last_astar_query_us / last_paired_htrip_query_us
    std::atomic<uint64_t> bfs_compare_count{0};          ///< Completed same-input H-TRIP vs baseline pairs

    struct BfsCompareSample {
        uint64_t job_id = 0;
        uint64_t tick = 0;
        size_t num_points = 0;
        double htrip_query_us = 0.0;
        double bfs_query_us = 0.0;
        double astar_query_us = 0.0;
        double speedup = 0.0;
        double astar_speedup = 0.0;
    };
    std::vector<BfsCompareSample> bfs_compare_samples; ///< Filled when both kernels finish; copy under bfs_mutex

    // Multi-robot coordination and benchmarking control flags
    uint64_t total_ticks = 0;        ///< Monotonically increasing simulation tick counter.
    int max_frontiers_eval = 64;     ///< Maximum frontier clusters included in m-TSP (Default: 64).
    int replan_period_ticks = 20;    ///< Periodic m-TSP re-planning interval (Default: 20 ticks = 1.0s at 20 Hz).
    int persistence_window_ticks = 50; ///< Minimum ticks (at 20Hz: 50 ticks = 2.5s) to persist on an assigned frontier before allowing preemption.
    int bfs_compare_interval = 0;    ///< BFS comparison cadence in H-TRIP queries (0 = never, N = evaluate every N queries).
    uint64_t last_replan_tick = 0;   ///< Simulation tick when m-TSP re-planning was last initiated.
    bool enable_bfs_baseline = false; ///< Enables concurrent k x BFS on a dedicated worker; never blocks H-TRIP/m-TSP.
    bool enable_astar_baseline = false; ///< Enables A* baseline comparison (disabled by default for speed).
    bool enable_fog_of_war = true;    ///< True if perception is unrevealed (fog-of-war); false if the floor is fully revealed.
    std::vector<Point> patrol_targets; ///< Active target waypoints when fog-of-war is disabled.
    size_t targets_visited_count = 0;  ///< Cumulative count of patrol/inspection targets reached by fleet in no-fog mode.
    BaselineExecutionMode baseline_mode = BaselineExecutionMode::Asynchronous; ///< Mode 1: Asynchronous (fluid), Mode 2: Synchronous (lock-step log).
    bool enable_async_updates = true; ///< Offloads H-TRIP dynamic index updates to a background worker thread.
    bool enable_live_frontier_projection = true; ///< Detects frontiers on live observed map and snaps to indexed map for H-TRIP queries.
    std::atomic<bool> is_updating_htrip{false}; ///< Thread-safe flag indicating background H-TRIP update worker is actively processing.
    std::atomic<bool> is_running_bfs{false};    ///< Thread-safe flag indicating background BFS comparison worker is active.
    std::atomic<bool> is_querying_pairwise{false}; ///< Thread-safe flag indicating background pairwise H-TRIP query worker is active.
    std::atomic<bool> new_tours_ready{false};       ///< Set to true when background query worker produces fresh staged tours.

    // Thread synchronization mutexes and double-buffering instances
    mutable std::mutex oracle_mutex;  ///< Synchronizes modifications to primary oracle and indexed_map during updates.
    mutable std::mutex commit_mutex;  ///< Protects completed_indexed_cells queue transferred from worker to main thread.
    mutable std::mutex query_mutex;   ///< Protects staged query distance matrix, fleet tours, and viewpoints.
    mutable std::mutex query_oracle_mutex; ///< Protects isolated double-buffered query oracle and query map snapshot.
    mutable std::mutex bfs_mutex;          ///< Serializes baseline job queues, worker dispatch, and timing aggregation.
    mutable std::mutex decay_mutex;        ///< Synchronizes visual decay counters for active tiles and index updates.

    MinPlusOracle query_oracle;            ///< Double-buffered H-TRIP oracle instance used by background query worker.
    GridMap query_map;                     ///< Double-buffered grid map snapshot corresponding to query_oracle.
    std::atomic<MinPlusOracle*> active_query_oracle_ptr{&query_oracle}; ///< Atomic pointer to active query oracle buffer.
    std::atomic<GridMap*> active_query_map_ptr{&query_map};             ///< Atomic pointer to active query map buffer.
    std::atomic<MinPlusOracle*> active_update_oracle_ptr{&oracle};      ///< Atomic pointer to background update oracle buffer.
    std::atomic<GridMap*> active_update_map_ptr{&indexed_map};          ///< Atomic pointer to background update map buffer.
    bool double_buffer_dirty = false;      ///< Legacy synchronization flag (kept for API compatibility).
    int grid_tile_size_b = 16;             ///< Cached tile dimension b for zero-lock dirty leaf index computation.
    int grid_num_leaf_tiles_x = 0;         ///< Cached count of leaf tiles along horizontal axis.
    int grid_num_leaf_tiles_y = 0;         ///< Cached count of leaf tiles along vertical axis.
    size_t grid_num_leaves = 0;            ///< Cached total count of leaf tiles in the hierarchy.

    // Pre-allocated zero-heap scratch buffers for hot query paths
    mutable MinPlusOracle::HTripQueryScratch async_query_scratch; ///< Caller-owned scratch buffer for async background queries.
    mutable MinPlusOracle::HTripQueryScratch sync_query_scratch;  ///< Caller-owned scratch buffer for synchronous main-thread queries.
    mutable std::vector<dist_t> async_query_matrix;///< Pre-allocated flat buffer storing pairwise distance matrix results.

    // Dedicated worker threads for non-blocking concurrent execution
    std::thread query_worker_thread;       ///< Single worker executing H-TRIP all-pairs queries and m-TSP tour optimization.
    std::thread update_worker_thread;      ///< Single worker executing dynamic H-TRIP min-plus leaf and internal updates.
    std::thread bfs_worker_thread;         ///< Single worker executing classical k x BFS and A* baseline benchmarks.
    std::atomic<uint64_t> bfs_job_seq{0};  ///< Monotonic sequence identifier for baseline benchmarking jobs.

    /**
     * @brief Context descriptor for a paired H-TRIP vs BFS/A* benchmark job.
     */
    struct BfsCompareJob {
        uint64_t job_id = 0;               ///< Unique sequential identifier for this job.
        uint64_t tick = 0;                 ///< Simulation tick at which the snapshot was captured.
        size_t num_points = 0;             ///< Number of viewpoints k evaluated in the k x k query.
        GridMap map;                       ///< Complete grid map snapshot used for graph search.
        std::vector<Point> pts;            ///< Set of k viewpoints queried.
        double htrip_query_us = 0.0;       ///< Wall-clock execution time for H-TRIP on this snapshot (microseconds).
        double bfs_query_us = 0.0;         ///< Wall-clock execution time for k x BFS on this snapshot (microseconds).
        double astar_query_us = 0.0;       ///< Wall-clock execution time for k(k-1)/2 A* searches (microseconds).
        bool occupied = false;             ///< True if this job slot holds an active benchmark task.
        bool has_htrip = false;            ///< True if H-TRIP timing has been recorded.
        bool has_bfs = false;              ///< True if BFS/A* timing has been recorded.
    };
    BfsCompareJob bfs_inflight_job;        ///< Baseline job currently being executed on the worker thread.
    BfsCompareJob bfs_pending_job;         ///< Queued baseline job waiting for worker availability.
    BfsCompareJob bfs_await_htrip_job;     ///< Job where BFS finished early, awaiting matching H-TRIP timing from query thread.

    std::vector<Point> completed_indexed_cells; ///< Free cells processed by update worker, pending commit to indexed_map.
    std::vector<Point> async_pending_cells;     ///< Batch of newly opened cells queued for the next background update run.
    std::vector<dist_t> staged_all_pairs_matrix;///< Pairwise distance matrix produced by query worker, pending main thread adoption.
    std::vector<RobotTour> staged_fleet_tours;  ///< Multi-robot exploration tours optimized by query worker.
    std::vector<Point> staged_viewpoints;       ///< Set of viewpoints corresponding to staged_fleet_tours.

    double htrip_query_freq_hz = 0.0;           ///< Measured operational frequency of H-TRIP all-pairs queries (Hz).
    double target_query_hz = 20.0;              ///< Target continuous query frequency (Hz, default: 20 Hz).
    uint64_t total_htrip_queries = 0;           ///< Total count of completed pairwise distance queries across mission.
    uint64_t last_query_tick = 0;               ///< Simulation tick when pairwise query was last dispatched.
    std::chrono::steady_clock::time_point last_pairwise_dispatch_time{}; ///< Timestamp when pairwise query was dispatched.
    std::chrono::steady_clock::time_point last_query_completed_time{};   ///< Timestamp when pairwise query finished.
    uint32_t robot_seed = 0;                    ///< Deterministic pseudorandom seed for robot pose and spawn generation.

    // Stationary / stuck robot detection and automatic recovery policy
    bool enable_auto_respawn = true;         ///< Enables automatic emergency respawning for immobilized robots.
    int auto_respawn_ticks = 100;            ///< Consecutive stuck ticks required before triggering respawn (default: 100 ticks = 5.0s).
    double stuck_dist_threshold_m = 0.20;    ///< Maximum net displacement threshold (meters) defining lack of translational progress.
    uint32_t total_robot_respawns = 0;       ///< Cumulative count of emergency respawns executed across entire fleet.
    MissionLogger mission_logger;            ///< High-precision telemetry and experiment trace logger.

    /**
     * @brief Thread-safe snapshot of recently updated tile decay counters for UI and video rendering.
     * @param out_tile Output vector receiving base tile decay counters (amber glow).
     * @param out_index Output vector receiving H-TRIP hierarchical update decay counters (emerald glow).
     */
    void getTileDecaySnapshots(std::vector<uint8_t>& out_tile, std::vector<uint8_t>& out_index) const {
        std::lock_guard<std::mutex> lock(decay_mutex);
        out_tile = active_tile_decay;
        out_index = active_index_decay;
    }

    /**
     * @brief Triggers an immediate one-shot k x BFS baseline benchmark on the current observed map.
     */
    void runBfsBaseline();

    /**
     * @brief Triggers an immediate one-shot k x BFS baseline benchmark on a provided map snapshot.
     */
    void runBfsBaseline(GridMap map_snapshot, std::vector<Point> pts);

    /**
     * @brief Runs synchronous H-TRIP, k x BFS, and A* baselines in lock-step on the calling thread (Mode 2).
     */
    void runSyncBaselines(const GridMap& map_snapshot, const std::vector<Point>& pts, double htrip_us, uint64_t tick);

    /**
     * @brief Flushes and synchronizes all pending baseline benchmarking jobs before session shutdown.
     */
    void drainBfsComparisons();

    /**
     * @brief Starts high-resolution telemetry recording for experiment benchmarking and paper figures.
     */
    bool startMissionLogging(const std::string& output_dir = "recordings/telemetry",
                             const std::string& session_prefix = "mission",
                             bool timestamp_subdir = true) {
        return mission_logger.startSession(output_dir, session_prefix, timestamp_subdir);
    }

    /**
     * @brief Terminates telemetry recording and generates comprehensive academic experiment summary JSON.
     */
    void stopMissionLogging() {
        drainBfsComparisons();
        mission_logger.endSession(*this);
    }

    /**
     * @brief Checks if mission telemetry logging is actively recording.
     */
    [[nodiscard]] bool isMissionLogging() const noexcept {
        return mission_logger.isLogging();
    }

    MultiRobotSimulator() = default;
    ~MultiRobotSimulator();
    MultiRobotSimulator(const MultiRobotSimulator&) = delete;
    MultiRobotSimulator& operator=(const MultiRobotSimulator&) = delete;

    /**
     * @brief Gracefully terminates and joins all background worker threads (update, query, BFS).
     */
    void stopBackgroundThreads() noexcept;

    /**
     * @brief Identifies leaf tiles affected by newly opened cells and updates dirty tile decay arrays.
     */
    void updateDirtyTileIndices(const std::vector<Point>& cells);

    /**
     * @brief Configures pseudorandom seed for deterministic initial pose generation.
     */
    void setSeed(uint32_t seed) noexcept { robot_seed = seed; }

    /**
     * @brief Retrieves active pseudorandom seed.
     */
    [[nodiscard]] uint32_t getSeed() const noexcept { return robot_seed; }

    /**
     * @brief Asynchronously dispatches H-TRIP pairwise distance evaluation and m-TSP tour optimization.
     */
    void startAsyncPairwiseQuery();

    /**
     * @brief Adopts newly computed m-TSP fleet tours while respecting persistence windows of active goals.
     */
    void applyFleetTours(std::vector<RobotTour>& tours, const std::vector<Point>& current_vps);

    /**
     * @brief Generates distributed patrol targets across the free space when fog-of-war is disabled.
     */
    void initPatrolTargets(size_t count);

    /**
     * @brief Replenishes visited patrol targets to maintain active target count.
     */
    void replenishPatrolTargets(size_t target_count);

    /**
     * @brief Initializes the simulator from a raw ground-truth grid map.
     *
     * Upscales grid by scale factor, allocates observed fog-of-war grid, initializes robots,
     * and constructs the initial H-TRIP hierarchical index.
     */
    void init(const GridMap& raw_gt,
              double raw_cell_size_m,
              int scale,
              int num_robots,
              const RobotConfig& cfg,
              int b = 16,
              int g = 2,
              uint32_t seed = 0,
              bool randomize_poses = true,
              bool enable_fog = true);

    /**
     * @brief Advances the continuous multi-robot simulation by dt seconds.
     *
     * 1. Evaluates autonomous guidance, obstacle repulsion, and recovery controllers.
     * 2. Integrates continuous SE(2) unicycle kinematics with continuous collision avoidance.
     * 3. Executes uniform polar LiDAR raymarching to reveal unknown cells into free/obstacle.
     * 4. Propagates revealed cells into H-TRIP via dynamic rank-1 Conway updates.
     * 5. Detects active frontier clusters and dispatches all-pairs shortest-path queries.
     */
    void step(double dt = 0.05);

    /**
     * @brief Resets fleet poses and fog-of-war while preserving loaded map geometry.
     */
    void reset(bool randomize_poses = true, uint32_t new_seed = 0);

    /**
     * @brief Dynamically resizes the active robot fleet.
     */
    void setRobotCount(int n, bool randomize_poses = false, uint32_t seed = 0);

    /**
     * @brief Updates LiDAR hardware parameters across all robots.
     */
    void updateLidarConfig(const LidarConfig& lcfg);

    /**
     * @brief Synchronously flushes all pending revealed free cells into the H-TRIP hierarchical index.
     */
    void flushPendingUpdates();

    /**
     * @brief Teleports an immobilized or stuck robot to a safe, clearance-verified passable cell.
     */
    bool respawnRobot(int robot_id);

    /**
     * @brief Monitors fleet progress and triggers emergency respawns for stationary robots.
     */
    int checkAndRespawnStuckRobots();

    /**
     * @brief Calculates cumulative exploration progress as a fraction in [0.0, 1.0].
     */
    [[nodiscard]] double explorationFraction() const noexcept {
        if (total_free_cells_gt <= 0) return 0.0;
        return static_cast<double>(explored_free_cells) / static_cast<double>(total_free_cells_gt);
    }

    /**
     * @brief Converts continuous world coordinates (meters) to discrete grid coordinates (row, col).
     */
    [[nodiscard]] Point worldToCell(double x, double y) const noexcept {
        int c = static_cast<int>(std::floor(x / cell_size_m));
        int r = static_cast<int>(std::floor(y / cell_size_m));
        c = std::clamp(c, 0, observed.cols - 1);
        r = std::clamp(r, 0, observed.rows - 1);
        return {r, c};
    }

    /**
     * @brief Converts discrete grid cell coordinates to continuous world center coordinates (meters).
     */
    [[nodiscard]] std::pair<double, double> cellToWorld(Point p) const noexcept {
        double x = (static_cast<double>(p.c) + 0.5) * cell_size_m;
        double y = (static_cast<double>(p.r) + 0.5) * cell_size_m;
        return {x, y};
    }

    /**
     * @brief Plans an obstacle-cleared, clearance-weighted A* shortest path on the observed grid.
     */
    [[nodiscard]] std::vector<Point> planPathAStar(Point start, Point goal) const;

    /**
     * @brief Determines whether all cells within a specified radius around point p are fully revealed.
     */
    [[nodiscard]] bool isFrontierAreaExplored(Point p, int radius_cells = 4) const noexcept;

private:
    /**
     * @brief Simulates planar LiDAR raymarching from robot pose across the ground-truth map.
     * Updates the observed grid map and records revealed free cells.
     */
    void castLidarRays(RobotState& robot);

    /**
     * @brief Computes forward linear velocity v and angular velocity omega via obstacle-repulsive carrot following.
     */
    void computeAutonomousControls(RobotState& robot, double dt);

    /**
     * @brief Integrates non-holonomic unicycle kinematics with swept-circle collision checks.
     */
    void integrateKinematics(RobotState& robot, double dt);

    /**
     * @brief Tests whether a circular robot footprint at (x, y) with radius r intersects any obstacle.
     */
    bool checkRobotCollision(double x, double y, double r) const noexcept;

    /**
     * @brief Tests whether a swept circular footprint along the segment from (x1, y1) to (x2, y2) is collision-free.
     */
    bool isSegmentClear(double x1, double y1, double x2, double y2, double r) const noexcept;

    /**
     * @brief Selects the closest reachable frontier for a single robot when fleet allocation is unassigned.
     */
    [[nodiscard]] Point findBestReachableFrontier(const RobotState& bot, std::vector<Point>& out_path) const;

    /**
     * @brief Dispatches m-TSP fleet planning or local heuristic assignment across active robots.
     */
    void assignFrontiersToFleet();

    /**
     * @brief Allocates and queues a paired benchmark job comparing H-TRIP against BFS/A*.
     */
    uint64_t createBfsCompareJob(GridMap map_snapshot, std::vector<Point> pts);

    /**
     * @brief Dispatches the background baseline worker thread for a queued comparison job.
     */
    void launchBfsWorker(GridMap map_snapshot, std::vector<Point> pts, uint64_t job_id,
                         uint64_t tick, size_t num_points, double htrip_query_us, bool has_htrip);

    /**
     * @brief Advances pending baseline jobs when the worker thread becomes idle.
     */
    void pumpBfsPending();

    /**
     * @brief Records H-TRIP query execution time for a pending paired comparison job.
     */
    void noteHtripQueryTime(uint64_t job_id, double htrip_us);

    /**
     * @brief Callback invoked when the baseline worker thread completes BFS and A* evaluations.
     */
    void onBfsWorkerFinished(uint64_t job_id, double bfs_us, double astar_us);

    /**
     * @brief Commits paired benchmarking timings to the sample history under mutex protection.
     */
    void commitBfsCompareLocked(uint64_t job_id, uint64_t tick, size_t num_points,
                                double htrip_us, double bfs_us, double astar_us);

    /**
     * @brief Computes the immediate next corridor cell along the shortest path toward goal.
     */
    [[nodiscard]] Point findNextStepToGoal(Point start, Point goal) const;

    /**
     * @brief Tests line-of-sight visibility between two grid cells via Bresenham traversal.
     */
    [[nodiscard]] bool hasLineOfSight(Point p1, Point p2) const noexcept;

    std::vector<int> tile_cell_counts;       ///< Internal workspace tracking cell allocations per tile.
    mutable std::vector<dist_t> corridor_dist;    ///< Scratch distance buffer for corridor Dijkstra searches.
    mutable std::vector<int> corridor_visited;    ///< Scratch visited buffer for corridor Dijkstra searches.
};

} // namespace htrip
