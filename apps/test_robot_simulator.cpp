#include "htrip/robot_simulator.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

int main() {
    std::cout << "========================================================\n";
    std::cout << "  RUNNING MULTI-ROBOT CONTINUOUS SIMULATOR UNIT TESTS\n";
    std::cout << "========================================================\n\n";

    // -------------------------------------------------------------------------
    // TEST 1: Resolution Scaling & Physical Metric Invariance
    // -------------------------------------------------------------------------
    std::cout << "[Test 1] Resolution Scaling & Metric Invariance...\n";
    {
        int raw_r = 16;
        int raw_c = 16;
        double raw_cell_size = 0.10; // 10 cm -> 1.6m x 1.6m
        htrip::GridMap raw_gt(raw_r, raw_c, htrip::CellState::Free);
        // Place a 4x4 obstacle block in the center (from 6 to 9)
        for (int r = 6; r < 10; ++r) {
            for (int c = 6; c < 10; ++c) {
                raw_gt.setObserved(r, c, htrip::CellState::Obstacle);
            }
        }

        for (int scale : {1, 2, 4}) {
            htrip::RobotConfig rcfg;
            htrip::MultiRobotSimulator sim;
            sim.init(raw_gt, raw_cell_size, scale, 1, rcfg);

            int expected_r = raw_r * scale;
            int expected_c = raw_c * scale;
            double expected_cell_size = raw_cell_size / static_cast<double>(scale);

            if (sim.ground_truth.rows != expected_r || sim.ground_truth.cols != expected_c) {
                std::cerr << "FAIL: Scale " << scale << " wrong rows/cols: "
                          << sim.ground_truth.rows << "x" << sim.ground_truth.cols << "\n";
                return 1;
            }
            if (std::abs(sim.cell_size_m - expected_cell_size) > 1e-6) {
                std::cerr << "FAIL: Scale " << scale << " wrong cell size: "
                          << sim.cell_size_m << " vs " << expected_cell_size << "\n";
                return 1;
            }
            // Invariance check
            if (std::abs(sim.world_w_m - 1.60) > 1e-6 || std::abs(sim.world_h_m - 1.60) > 1e-6) {
                std::cerr << "FAIL: Scale " << scale << " physical world dimension altered: "
                          << sim.world_w_m << "x" << sim.world_h_m << " (expected 1.60x1.60)\n";
                return 1;
            }

            // Verify obstacle block scaled properly: center cell must be Obstacle
            htrip::Point center_pt = sim.worldToCell(0.80, 0.80);
            if (sim.ground_truth.getObserved(center_pt) != htrip::CellState::Obstacle) {
                std::cerr << "FAIL: Scale " << scale << " center obstacle missing!\n";
                return 1;
            }
        }
        std::cout << "  -> PASSED: Metric dimensions strictly invariant across 1x, 2x, 4x scaling.\n";
    }

    // -------------------------------------------------------------------------
    // TEST 2: LiDAR Ray Uniform Polar Distribution & Beam Visibility
    // -------------------------------------------------------------------------
    std::cout << "\n[Test 2] Uniform Polar LiDAR Beam Distribution...\n";
    {
        htrip::LidarConfig lcfg_360;
        lcfg_360.span_deg = 360.0;
        lcfg_360.num_rays = 360;
        if (std::abs(lcfg_360.angularResolutionDeg() - 1.0) > 1e-5) {
            std::cerr << "FAIL: 360 deg / 360 rays resolution != 1.0 deg\n";
            return 1;
        }

        htrip::LidarConfig lcfg_180;
        lcfg_180.span_deg = 180.0;
        lcfg_180.num_rays = 181;
        if (std::abs(lcfg_180.angularResolutionDeg() - 1.0) > 1e-5) {
            std::cerr << "FAIL: 180 deg / 181 rays resolution != 1.0 deg\n";
            return 1;
        }

        // Test raymarching hitting wall at known distance
        htrip::GridMap box(32, 32, htrip::CellState::Free);
        // Wall at column 20
        for (int r = 0; r < 32; ++r) {
            box.setObserved(r, 20, htrip::CellState::Obstacle);
        }

        htrip::RobotConfig rcfg;
        rcfg.lidar.range_m = 10.0;
        rcfg.lidar.span_deg = 360.0;
        rcfg.lidar.num_rays = 360;

        htrip::MultiRobotSimulator sim;
        sim.init(box, 0.10, 1, 1, rcfg);

        // Place robot at (1.0m, 1.6m) facing east (theta = 0 rad, towards wall at x = 2.0m)
        sim.robots[0].x = 1.0;
        sim.robots[0].y = 1.6;
        sim.robots[0].theta = 0.0;
        sim.step(0.01);

        // Find beam along theta = 0
        bool found_wall_hit = false;
        for (const auto& beam : sim.robots[0].scan) {
            if (std::abs(beam.angle_rad) < 0.02) {
                found_wall_hit = true;
                // Wall is at 2.0m, robot is at 1.0m -> distance ~ 1.0m
                if (!beam.hit_obstacle) {
                    std::cerr << "FAIL: Eastward LiDAR ray missed obstacle wall!\n";
                    return 1;
                }
                if (std::abs(beam.distance_m - 1.0) > 0.15) {
                    std::cerr << "FAIL: LiDAR measured dist " << beam.distance_m << " != 1.0m\n";
                    return 1;
                }
            }
        }
        if (!found_wall_hit) {
            std::cerr << "FAIL: Zero-degree beam not found in scan!\n";
            return 1;
        }
        std::cout << "  -> PASSED: 360 uniformly spaced rays correctly measured wall distance.\n";
    }

    // -------------------------------------------------------------------------
    // TEST 3: Multi-Robot Exploration & Progressive Fog-of-War Revelation
    // -------------------------------------------------------------------------
    std::cout << "\n[Test 3] Multi-Robot Progressive Exploration & Collision Safety...\n";
    {
        // 32x32 maze with rooms
        htrip::GridMap maze(32, 32, htrip::CellState::Free);
        // Outer boundaries
        for (int i = 0; i < 32; ++i) {
            maze.setObserved(0, i, htrip::CellState::Obstacle);
            maze.setObserved(31, i, htrip::CellState::Obstacle);
            maze.setObserved(i, 0, htrip::CellState::Obstacle);
            maze.setObserved(i, 31, htrip::CellState::Obstacle);
        }
        // Central wall with doorway
        for (int r = 0; r < 32; ++r) {
            if (r != 16) maze.setObserved(r, 16, htrip::CellState::Obstacle);
        }

        htrip::RobotConfig rcfg;
        rcfg.radius_m = 0.15;
        rcfg.v_max = 0.5;
        rcfg.omega_max = 1.5;
        rcfg.lidar.range_m = 3.0;
        rcfg.lidar.span_deg = 360.0;
        rcfg.lidar.num_rays = 120;

        htrip::MultiRobotSimulator sim;
        // Scale 2x -> 64x64
        sim.init(maze, 0.10, 2, 3, rcfg, 16, 2);

        int init_explored = sim.explored_free_cells;
        if (init_explored <= 0) {
            std::cerr << "FAIL: Initial scan revealed 0 cells!\n";
            return 1;
        }

        // Run 40 steps of simulation
        for (int step = 0; step < 40; ++step) {
            sim.step(0.05);

            // Verify no robot is in obstacle collision
            for (const auto& bot : sim.robots) {
                htrip::Point cp = sim.worldToCell(bot.x, bot.y);
                if (sim.ground_truth.getObserved(cp) == htrip::CellState::Obstacle) {
                    std::cerr << "FAIL: Robot " << bot.id << " collided into obstacle at ("
                              << bot.x << ", " << bot.y << ")!\n";
                    return 1;
                }
            }
        }

        int final_explored = sim.explored_free_cells;
        if (final_explored <= init_explored) {
            std::cerr << "FAIL: Explored cell count did not expand during wandering: "
                      << init_explored << " -> " << final_explored << "\n";
            return 1;
        }

        std::cout << "  -> Initial Explored Cells: " << init_explored << "\n";
        std::cout << "  -> Explored Cells after 40 Steps: " << final_explored
                  << " (" << (sim.explorationFraction() * 100.0) << "% of map)\n";
        std::cout << "  -> Active Frontiers Detected: " << sim.frontiers.size() << "\n";
        std::cout << "  -> PASSED: Robots wandered naturally without collisions and progressively cleared fog-of-war.\n";
    }

    // -------------------------------------------------------------------------
    // TEST 4: Stuck Robot Detection & Auto-Respawn (100 Steps)
    // -------------------------------------------------------------------------
    std::cout << "\n[Test 4] Stuck Robot Detection & Auto-Respawn (100 Steps)...\n";
    {
        htrip::GridMap box(32, 32, htrip::CellState::Free);
        for (int i = 0; i < 32; ++i) {
            box.setObserved(0, i, htrip::CellState::Obstacle);
            box.setObserved(31, i, htrip::CellState::Obstacle);
            box.setObserved(i, 0, htrip::CellState::Obstacle);
            box.setObserved(i, 31, htrip::CellState::Obstacle);
        }

        htrip::RobotConfig rcfg;
        rcfg.radius_m = 0.15;
        rcfg.v_max = 0.5;
        rcfg.omega_max = 1.5;

        htrip::MultiRobotSimulator sim;
        sim.init(box, 0.10, 1, 1, rcfg, 16, 2);

        sim.enable_auto_respawn = true;
        sim.auto_respawn_ticks = 100;

        double initial_x = sim.robots[0].x;
        double initial_y = sim.robots[0].y;

        // Run 99 steps holding robot near initial position
        for (int step = 0; step < 99; ++step) {
            sim.robots[0].x = initial_x;
            sim.robots[0].y = initial_y;
            sim.robots[0].v = 0.0;
            sim.step(0.05);
        }

        if (sim.total_robot_respawns != 0) {
            std::cerr << "FAIL: Robot respawned prematurely before 100 steps!\n";
            return 1;
        }
        if (sim.robots[0].stuck_translation_ticks < 95) {
            std::cerr << "FAIL: stuck_translation_ticks (" << sim.robots[0].stuck_translation_ticks
                      << ") did not increment properly!\n";
            return 1;
        }

        // Run the 100th step (forcing it to remain in place during integration)
        sim.robots[0].x = initial_x;
        sim.robots[0].y = initial_y;
        sim.robots[0].v = 0.0;
        sim.step(0.05);

        if (sim.total_robot_respawns < 1) {
            std::cerr << "FAIL: Robot was NOT respawned after 100 stuck steps!\n";
            return 1;
        }
        if (sim.robots[0].stuck_translation_ticks != 0) {
            std::cerr << "FAIL: stuck_translation_ticks was not reset to 0 after respawn!\n";
            return 1;
        }

        double new_x = sim.robots[0].x;
        double new_y = sim.robots[0].y;
        double dist = std::hypot(new_x - initial_x, new_y - initial_y);
        if (dist < 0.50) {
            std::cerr << "FAIL: Respawned position is too close to original stuck position: " << dist << " m\n";
            return 1;
        }

        std::cout << "  -> Initial Stuck Position: (" << initial_x << ", " << initial_y << ")\n";
        std::cout << "  -> Respawned Position: (" << new_x << ", " << new_y << ") - Relocation dist: " << dist << " m\n";
        std::cout << "  -> Total Respawns: " << sim.total_robot_respawns << "\n";
        std::cout << "  -> PASSED: Robot successfully detected as stuck and relocated to a safe position.\n";
    }

    // -------------------------------------------------------------------------
    // TEST 5: Procedural Room Maze Single-Component Connectivity
    // -------------------------------------------------------------------------
    std::cout << "\n[Test 5] Procedural Room Maze Single-Component Connectivity...\n";
    {
        htrip::GridMap grid(64, 64, htrip::CellState::Obstacle);
        grid.generateProceduralRooms(64, 64, 16, 42);

        int total_free = 0;
        htrip::Point start_free = {-1, -1};
        for (int r = 0; r < 64; ++r) {
            for (int c = 0; c < 64; ++c) {
                if (grid.getObserved(r, c) == htrip::CellState::Free) {
                    total_free++;
                    if (start_free.r < 0) start_free = {r, c};
                }
            }
        }

        if (total_free == 0 || start_free.r < 0) {
            std::cerr << "FAIL: Procedural rooms generated 0 free cells!\n";
            return 1;
        }

        // BFS flood fill to verify all rooms are connected into 1 component
        std::vector<bool> visited(64 * 64, false);
        std::vector<htrip::Point> q;
        q.push_back(start_free);
        visited[start_free.r * 64 + start_free.c] = true;
        int reachable_free = 0;

        size_t head = 0;
        const int dr[4] = {-1, 1, 0, 0};
        const int dc[4] = {0, 0, -1, 1};

        while (head < q.size()) {
            auto p = q[head++];
            reachable_free++;
            for (int i = 0; i < 4; ++i) {
                int nr = p.r + dr[i];
                int nc = p.c + dc[i];
                if (grid.inBounds(nr, nc) && grid.getObserved(nr, nc) == htrip::CellState::Free) {
                    int idx = nr * 64 + nc;
                    if (!visited[idx]) {
                        visited[idx] = true;
                        q.push_back({nr, nc});
                    }
                }
            }
        }

        if (reachable_free != total_free) {
            std::cerr << "FAIL: Procedural rooms are DISCONNECTED! Reachable: "
                      << reachable_free << " vs Total: " << total_free << "\n";
            return 1;
        }

        std::cout << "  -> Total Free Cells in 4x4 Room Grid: " << total_free << "\n";
        std::cout << "  -> Reachable from Room (0,0): " << reachable_free << " (100% connected)\n";
        std::cout << "  -> PASSED: All doorways cleanly pierced; exactly 1 connected component verified.\n";
    }

    // -------------------------------------------------------------------------
    // TEST 6: Idle robots still receive a frontier when K < N or all K sit near one robot
    // -------------------------------------------------------------------------
    std::cout << "\n[Test 6] Fleet coverage fill for starved robots...\n";
    {
        using htrip::Point;
        using htrip::dist_t;
        using htrip::INF;
        using htrip::MTSPPlanner;
        using htrip::RobotTour;

        auto run_case = [](const std::vector<Point>& robots,
                           const std::vector<Point>& frontiers,
                           const std::vector<dist_t>& D,
                           const char* name) -> int {
            std::vector<RobotTour> tours;
            MTSPPlanner::planFleetTours(robots, frontiers, D, tours, 32);
            if (tours.size() != robots.size()) {
                std::cerr << "FAIL: " << name << " tour count " << tours.size() << "\n";
                return 1;
            }
            for (size_t i = 0; i < tours.size(); ++i) {
                if (tours[i].empty()) {
                    std::cerr << "FAIL: " << name << " robot " << i
                              << " has an empty tour while frontiers remain\n";
                    return 1;
                }
            }
            return 0;
        };

        // 8 robots, 3 frontiers all nearest to robot 0 in geodesic distance.
        {
            const size_t N = 8, K = 3, M = N + K;
            std::vector<Point> robots(N), frontiers(K);
            for (size_t i = 0; i < N; ++i) robots[i] = {0, static_cast<int>(i) * 10};
            for (size_t f = 0; f < K; ++f) frontiers[f] = {0, 100 + static_cast<int>(f)};
            std::vector<dist_t> D(M * M, INF);
            for (size_t i = 0; i < M; ++i) D[i * M + i] = 0;
            for (size_t f = 0; f < K; ++f) {
                D[0 * M + (N + f)] = 5;
                D[(N + f) * M + 0] = 5;
                for (size_t r = 1; r < N; ++r) {
                    D[r * M + (N + f)] = static_cast<dist_t>(50 + r);
                    D[(N + f) * M + r] = static_cast<dist_t>(50 + r);
                }
            }
            if (run_case(robots, frontiers, D, "clustered-K<N") != 0) return 1;
        }

        // 4 robots, 10 frontiers all closest to robot 0; remaining robots must still be tasked.
        {
            const size_t N = 4, K = 10, M = N + K;
            std::vector<Point> robots(N), frontiers(K);
            for (size_t i = 0; i < N; ++i) robots[i] = {static_cast<int>(i), 0};
            for (size_t f = 0; f < K; ++f) frontiers[f] = {0, 20 + static_cast<int>(f)};
            std::vector<dist_t> D(M * M, INF);
            for (size_t i = 0; i < M; ++i) D[i * M + i] = 0;
            for (size_t f = 0; f < K; ++f) {
                D[0 * M + (N + f)] = 3;
                D[(N + f) * M + 0] = 3;
                for (size_t r = 1; r < N; ++r) {
                    D[r * M + (N + f)] = static_cast<dist_t>(30 + r);
                    D[(N + f) * M + r] = static_cast<dist_t>(30 + r);
                }
            }
            if (run_case(robots, frontiers, D, "clustered-K>N") != 0) return 1;
        }

        std::cout << "  -> PASSED: Every robot receives a frontier when any remain.\n";
    }

    std::cout << "\n========================================================\n";
    std::cout << "  ALL MULTI-ROBOT SIMULATOR UNIT TESTS PASSED!\n";
    std::cout << "========================================================\n";
    return 0;
}
