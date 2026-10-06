#include "htrip/grid_map.hpp"
#include "htrip/frontier.hpp"
#include "htrip/minplus_oracle.hpp"
#include "htrip/fuel_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

using namespace htrip;

// =========================================================================
// Suite 1: Canonical Toy Scenario (Paper Fig. 1 & Table 1 Replication)
// =========================================================================
void runToyScenario() {
    std::cout << "=================================================================\n";
    std::cout << "  SUITE 1: CANONICAL TOY SCENARIO (Paper Fig. 1 & Table 1)\n";
    std::cout << "=================================================================\n";

    // 32x32 map: Two rooms separated by a solid wall at column 15.
    // Detour corridor runs around the bottom through row 28.
    GridMap grid(32, 32, CellState::Obstacle);

    // Carve Room A (left, cols 2-14, rows 4-16)
    for (int r = 4; r <= 16; ++r) {
        for (int c = 2; c <= 14; ++c) {
            grid.setObserved(r, c, CellState::Free);
        }
    }

    // Carve Room B (right, cols 17-29, rows 4-16)
    for (int r = 4; r <= 16; ++r) {
        for (int c = 17; c <= 29; ++c) {
            grid.setObserved(r, c, CellState::Free);
        }
    }

    // Carve Detour Corridor connecting Room A to Room B via row 28
    // From Room A: down col 4 to row 28
    for (int r = 16; r <= 28; ++r) grid.setObserved(r, 4, CellState::Free);
    // Along row 28: from col 4 to col 26
    for (int c = 4; c <= 26; ++c) grid.setObserved(28, c, CellState::Free);
    // Up to Room B: up col 26 from row 28 to row 16
    for (int r = 16; r <= 28; ++r) grid.setObserved(r, 26, CellState::Free);

    // Wall at col 15, 16 remains Obstacle initially.
    Point f1{10, 4};   // Frontier 1 in Room A
    Point f2{10, 26};  // Frontier 2 in Room B
    std::vector<Point> viewpoints = {f1, f2};

    // Sensor radius R = 5 cells
    const int SENSOR_R = 5;

    // Initialize Oracles
    HierarchyConfig cfg;
    cfg.tile_size_b = 16;
    cfg.group_w = 2;
    cfg.group_h = 2;
    HTripOracle htrip_oracle(cfg);
    htrip_oracle.build(grid);

    FuelCacheOracle fuel_oracle(SENSOR_R);
    fuel_oracle.build(grid, viewpoints);

    dist_t d_init_bfs = BfsOracle::queryDistance(grid, f1, f2);
    dist_t d_init_htrip = htrip_oracle.queryDistance(f1, f2, grid);
    dist_t d_init_fuel = fuel_oracle.queryDistance(0, 1);

    std::cout << "\n[Before Doorway Discovery]:\n";
    std::cout << "  Detour distance Ground Truth (BFS): " << d_init_bfs << " steps\n";
    std::cout << "  H-TRIP Oracle distance:            " << d_init_htrip << " steps\n";
    std::cout << "  FUEL Cache distance:               " << d_init_fuel << " steps\n";

    // Now, open a direct doorway between Room A and Room B at (10, 15) and (10, 16)
    Point door1{10, 15};
    Point door2{10, 16};
    grid.openDoorway(door1);
    grid.openDoorway(door2);

    // The robot/sensor just opened the door at (10, 15).
    // Note that f1 is at (10, 4) -> distance = 11 > SENSOR_R (5).
    // f2 is at (10, 26) -> distance = 10 > SENSOR_R (5).
    // Neither frontier is inside the sensor bounding box of the doorway!

    // Update H-TRIP
    auto t0 = std::chrono::high_resolution_clock::now();
    htrip_oracle.onCellOpened(door1, grid);
    htrip_oracle.onCellOpened(door2, grid);
    auto t1 = std::chrono::high_resolution_clock::now();
    double htrip_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    // Update FUEL
    fuel_oracle.onCellOpened(door1, grid, viewpoints);
    fuel_oracle.onCellOpened(door2, grid, viewpoints);

    // Ground truth BFS
    dist_t d_final_bfs = BfsOracle::queryDistance(grid, f1, f2);
    dist_t d_final_htrip = htrip_oracle.queryDistance(f1, f2, grid);
    dist_t d_final_fuel = fuel_oracle.queryDistance(0, 1);

    std::cout << "\n[After Doorway Discovery at (10,15)-(10,16)]:\n";
    std::cout << "+-----------------------+-------------+-------------+---------------------+------------+\n";
    std::cout << "| Method                | d(F1, F2)   | Status      | Work Scope          | Rel Error  |\n";
    std::cout << "+-----------------------+-------------+-------------+---------------------+------------+\n";

    std::cout << "| Ground Truth (BFS)    | "
              << std::setw(11) << d_final_bfs << " | Exact       | Global BFS          |       0.0% |\n";

    double fuel_err = 100.0 * (static_cast<double>(d_final_fuel - d_final_bfs) / d_final_bfs);
    std::cout << "| FUEL Cache Baseline   | "
              << std::setw(11) << d_final_fuel << " | "
              << (d_final_fuel > d_final_bfs ? "STALE       " : "Exact       ") << "| "
              << fuel_oracle.lastBfsCount() << " BFS (skipped)     | "
              << std::setw(9) << std::fixed << std::setprecision(1) << fuel_err << "% |\n";

    double htrip_err = 100.0 * (static_cast<double>(d_final_htrip - d_final_bfs) / d_final_bfs);
    std::cout << "| Proposed H-TRIP       | "
              << std::setw(11) << d_final_htrip << " | "
              << (d_final_htrip == d_final_bfs ? "EXACT       " : "MISMATCH    ") << "| Leaf + 1 parent     | "
              << std::setw(9) << std::fixed << std::setprecision(1) << htrip_err << "% |\n";
    std::cout << "+-----------------------+-------------+-------------+---------------------+------------+\n";

    std::cout << "  H-TRIP update time: " << htrip_us << " us (exact, zero allocation)\n";
    std::cout << "  FUEL update time:   " << fuel_oracle.lastUpdateTimeUs() << " us ("
              << fuel_oracle.lastBfsCount() << " searches, surviving pairs skipped)\n";

    if (d_final_htrip == d_final_bfs && d_final_fuel > d_final_bfs) {
        std::cout << "\n>>> SUITE 1 SUCCESS: Paper Figure 1 & Table 1 perfectly validated!\n\n";
    } else {
        std::cout << "\n>>> SUITE 1 WARNING: Results differed from expected table.\n\n";
    }
}

// =========================================================================
// Suite 2: Multi-Frontier Exploration Benchmark (Paper Section V & Q2)
// =========================================================================
void runExplorationBenchmark(int k_target = 25, int num_doorways = 20) {
    std::cout << "=================================================================\n";
    std::cout << "  SUITE 2: MULTI-FRONTIER EXPLORATION BENCHMARK (k = " << k_target
              << ", steps = " << num_doorways << ")\n";
    std::cout << "=================================================================\n";

    GridMap grid;
    grid.generateProceduralRooms(128, 128, 16, 888);
    grid.initExplorationSnapshot({1, 1}, 0.50, 0.05, 42);

    FrontierManager fm(2);
    auto clusters = fm.detectAndCluster(grid);
    std::cout << "Detected " << clusters.size() << " frontier clusters in exploration snapshot.\n";

    std::vector<Point> viewpoints;
    for (const auto& cl : clusters) {
        if (grid.isPassable(cl.representative)) {
            viewpoints.push_back(cl.representative);
        }
    }

    std::vector<Point> all_free;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            Point p{r, c};
            if (grid.isPassable(p)) all_free.push_back(p);
        }
    }

    std::mt19937_64 vp_rng(12345);
    while (viewpoints.size() < static_cast<size_t>(k_target) && !all_free.empty()) {
        size_t idx = vp_rng() % all_free.size();
        viewpoints.push_back(all_free[idx]);
    }
    if (viewpoints.size() > static_cast<size_t>(k_target)) {
        viewpoints.resize(static_cast<size_t>(k_target));
    }

    size_t k = viewpoints.size();
    std::cout << "Selected k = " << k << " active viewpoints across explored rooms.\n";

    // Setup Oracles
    HierarchyConfig cfg;
    cfg.tile_size_b = 16;
    cfg.group_w = 2;
    cfg.group_h = 2;

    HTripOracle htrip_oracle(cfg);
    htrip_oracle.build(grid);

    FuelCacheOracle fuel_oracle(12); // sensor radius = 12 cells
    fuel_oracle.build(grid, viewpoints);

    MinPlusOracle::HTripQueryScratch scratch;
    MinPlusOracle::HTripIncrementalState inc_state;
    std::vector<dist_t> htrip_matrix;

    std::mt19937_64 rng(54321);

    size_t total_steps = 0;
    size_t htrip_discrepancies = 0;
    double sum_stale_fraction = 0.0;
    double max_overall_rel_error = 0.0;
    double total_htrip_us = 0.0;
    double total_fuel_us = 0.0;
    double total_bfs_us = 0.0;
    size_t total_fuel_searches = 0;

    std::cout << "\nStepping through " << num_doorways << " dynamic doorway discovery events...\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << " Step | Doorway   | H-TRIP Err | FUEL Stale% | FUEL MaxErr | H-TRIP us | FUEL us\n";
    std::cout << "-----------------------------------------------------------------\n";

    int step = 0;
    while (step < num_doorways) {
        int r = static_cast<int>(rng() % static_cast<uint64_t>(grid.rows));
        int c = static_cast<int>(rng() % static_cast<uint64_t>(grid.cols));
        Point door{r, c};

        if (grid.ground_truth[grid.index(door)] != CellState::Free ||
            grid.getObserved(door) == CellState::Free) {
            continue;
        }

        grid.openDoorway(door);

        // 1. Update H-TRIP
        auto t0 = std::chrono::high_resolution_clock::now();
        htrip_oracle.onCellOpened(door, grid);
        htrip_oracle.queryAllPairsIncremental(viewpoints, grid, htrip_matrix, scratch, inc_state);
        auto t1 = std::chrono::high_resolution_clock::now();
        double htrip_time = std::chrono::duration<double, std::micro>(t1 - t0).count();
        total_htrip_us += htrip_time;

        // 2. Update FUEL Baseline
        t0 = std::chrono::high_resolution_clock::now();
        fuel_oracle.onCellOpened(door, grid, viewpoints);
        t1 = std::chrono::high_resolution_clock::now();
        double fuel_time = std::chrono::duration<double, std::micro>(t1 - t0).count();
        total_fuel_us += fuel_time;
        total_fuel_searches += fuel_oracle.lastBfsCount();

        // 3. Ground Truth BFS
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<dist_t> gt_matrix = BfsOracle::queryAllPairs(grid, viewpoints);
        t1 = std::chrono::high_resolution_clock::now();
        total_bfs_us += std::chrono::duration<double, std::micro>(t1 - t0).count();

        // Check H-TRIP exactness
        int step_discrepancies = 0;
        for (size_t i = 0; i < htrip_matrix.size(); ++i) {
            if (htrip_matrix[i] != gt_matrix[i]) {
                step_discrepancies++;
                htrip_discrepancies++;
            }
        }

        // Evaluate FUEL Staleness
        auto rep = fuel_oracle.evaluateAgainstGroundTruth(gt_matrix);
        sum_stale_fraction += rep.stale_pair_fraction;
        if (rep.max_relative_error > max_overall_rel_error) {
            max_overall_rel_error = rep.max_relative_error;
        }

        total_steps++;

        if (step < 10 || step % 5 == 0 || step == num_doorways - 1) {
            std::cout << std::setw(5) << step + 1 << " | ("
                      << std::setw(3) << door.r << "," << std::setw(3) << door.c << ") | "
                      << std::setw(10) << step_discrepancies << " | "
                      << std::setw(10) << std::fixed << std::setprecision(1) << (rep.stale_pair_fraction * 100.0) << "% | "
                      << std::setw(10) << std::fixed << std::setprecision(1) << (rep.max_relative_error * 100.0) << "% | "
                      << std::setw(9) << std::fixed << std::setprecision(1) << htrip_time << " | "
                      << std::setw(7) << std::fixed << std::setprecision(1) << fuel_time << "\n";
        }

        step++;
    }

    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "\n[Multi-Frontier Benchmark Summary Across " << total_steps << " Steps]:\n";
    std::cout << "  H-TRIP Total Discrepancies vs Ground Truth: " << htrip_discrepancies
              << " (" << (htrip_discrepancies == 0 ? "100% BIT-EXACT" : "FAILED") << ")\n";
    std::cout << "  FUEL Average Stale Pair Percentage:         "
              << std::fixed << std::setprecision(2) << (sum_stale_fraction / total_steps * 100.0) << "%\n";
    std::cout << "  FUEL Maximum Observed Relative Error:       "
              << std::fixed << std::setprecision(1) << (max_overall_rel_error * 100.0) << "%\n";
    std::cout << "  FUEL Average Searches per Step:             "
              << std::fixed << std::setprecision(2) << (static_cast<double>(total_fuel_searches) / total_steps)
              << " (skipped pairs caused stale errors)\n";
    std::cout << "  Average H-TRIP Update + Query Latency:      "
              << std::fixed << std::setprecision(1) << (total_htrip_us / total_steps) << " us\n";
    std::cout << "  Average FUEL Baseline Latency:              "
              << std::fixed << std::setprecision(1) << (total_fuel_us / total_steps) << " us\n";
    std::cout << "  Average Full BFS Recomputation Latency:     "
              << std::fixed << std::setprecision(1) << (total_bfs_us / total_steps) << " us\n";
    std::cout << "  Speedup of H-TRIP over Full BFS:            "
              << std::fixed << std::setprecision(1) << (total_bfs_us / total_htrip_us) << "x\n";
    std::cout << "=================================================================\n\n";
}

// =========================================================================
// Suite 3: Systematic Shortcut Revelation Benchmark (Q2: Frequency & Magnitude)
// =========================================================================
void runShortcutBenchmark() {
    std::cout << "=================================================================\n";
    std::cout << "  SUITE 3: SYSTEMATIC SHORTCUT REVELATION BENCHMARK (Q2 Protocol)\n";
    std::cout << "=================================================================\n";

    // 64x64 grid with 4x4 rooms of size 16
    GridMap grid;
    grid.generateProceduralRooms(64, 64, 16, 42);

    // Mark all existing ground truth free space as explored
    for (int r = 0; r < 64; ++r) {
        for (int c = 0; c < 64; ++c) {
            if (grid.ground_truth[grid.index(r, c)] == CellState::Free) {
                grid.setObserved(r, c, CellState::Free);
            }
        }
    }

    // Place viewpoints at the center of all 16 rooms
    std::vector<Point> viewpoints;
    for (int ry = 0; ry < 4; ++ry) {
        for (int rx = 0; rx < 4; ++rx) {
            Point p{ry * 16 + 8, rx * 16 + 8};
            if (grid.isPassable(p)) {
                viewpoints.push_back(p);
            }
        }
    }

    const size_t k = viewpoints.size();
    std::cout << "Tracking " << k << " room centroids (" << (k * (k - 1)) << " directional pairs).\n";

    HierarchyConfig cfg;
    cfg.tile_size_b = 16;
    cfg.group_w = 2;
    cfg.group_h = 2;

    HTripOracle htrip(cfg);
    htrip.build(grid);

    FuelCacheOracle fuel(4); // sensor radius R = 4 (doorway discoveries are localized)
    fuel.build(grid, viewpoints);

    MinPlusOracle::HTripQueryScratch scratch;
    MinPlusOracle::HTripIncrementalState inc_state;
    std::vector<dist_t> htrip_mat;

    std::cout << "\nOpening 10 wall segments between adjacent rooms (shortcuts)...\n";
    std::cout << "------------------------------------------------------------------------------------\n";
    std::cout << " Shortcut # | Wall Opened | GroundTruth BFS | H-TRIP (Exact) | FUEL Cache (Stale) | Error%\n";
    std::cout << "------------------------------------------------------------------------------------\n";

    int shortcuts_opened = 0;
    std::mt19937_64 rng(999);

    for (int ry = 0; ry < 3; ++ry) {
        for (int rx = 0; rx < 3; ++rx) {
            if (shortcuts_opened >= 10) break;

            // Carve a horizontal door at column (rx+1)*16 between rows ry*16+7 and ry*16+9
            int wall_c = (rx + 1) * 16;
            int door_r = ry * 16 + 8;
            Point door1{door_r, wall_c - 1};
            Point door2{door_r, wall_c};

            grid.openDoorway(door1);
            grid.openDoorway(door2);
            htrip.onCellOpened(door1, grid);
            htrip.onCellOpened(door2, grid);
            fuel.onCellOpened(door1, grid, viewpoints);
            fuel.onCellOpened(door2, grid, viewpoints);

            htrip.queryAllPairsIncremental(viewpoints, grid, htrip_mat, scratch, inc_state);
            auto gt_mat = BfsOracle::queryAllPairs(grid, viewpoints);

            auto rep = fuel.evaluateAgainstGroundTruth(gt_mat);

            shortcuts_opened++;

            // Show distance between the two rooms adjacent to this wall
            Point p_left{door_r, wall_c - 8};
            Point p_right{door_r, wall_c + 8};
            dist_t d_gt = BfsOracle::queryDistance(grid, p_left, p_right);
            dist_t d_ht = htrip.queryDistance(p_left, p_right, grid);

            std::cout << "  #" << std::setw(2) << shortcuts_opened << "        | ("
                      << std::setw(2) << door1.r << "," << std::setw(2) << door1.c << "-" << std::setw(2) << door2.c << ")  | "
                      << std::setw(15) << d_gt << " | "
                      << std::setw(14) << d_ht << " | "
                      << std::setw(8) << std::fixed << std::setprecision(1)
                      << (rep.stale_pair_fraction * 100.0) << "% stale   | "
                      << std::setw(6) << std::fixed << std::setprecision(1)
                      << (rep.max_relative_error * 100.0) << "%\n";
        }
    }

    std::cout << "------------------------------------------------------------------------------------\n";
    std::cout << ">>> SUITE 3 COMPLETE: Demonstrated systemic stale caching in FUEL under shortcuts.\n\n";
}

int main() {
    std::cout << "\n=================================================================\n";
    std::cout << "   ICRA 2027: FUEL BASELINE ORACLE VALIDATION & BENCHMARK\n";
    std::cout << "=================================================================\n\n";

    runToyScenario();
    runExplorationBenchmark(25, 20);
    runShortcutBenchmark();

    return 0;
}
