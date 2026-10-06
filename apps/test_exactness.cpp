#include "htrip/grid_map.hpp"
#include "htrip/frontier.hpp"
#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include <iostream>
#include <cassert>
#include <random>
#include <algorithm>
#include <iomanip>
#include <vector>

using namespace htrip;

// Helper to run full matrix verification and report stats
bool verifyAllPairsMatrix(const std::string& test_name,
                          const MinPlusOracle& oracle,
                          const GridMap& grid,
                          const std::vector<Point>& viewpoints) {
    size_t k = viewpoints.size();
    if (k == 0) return true;

    auto prop_mat = oracle.queryAllPairs(viewpoints, grid);
    auto bfs_mat = BfsOracle::queryAllPairs(grid, viewpoints);

    int total_queries = static_cast<int>(k * k);
    int exact_matches = 0;
    int max_diff = 0;

    for (size_t i = 0; i < k * k; ++i) {
        dist_t d_prop = prop_mat[i];
        dist_t d_bfs = bfs_mat[i];

        int diff = std::abs(static_cast<int>(d_prop) - static_cast<int>(d_bfs));
        if (diff > max_diff) max_diff = diff;
        if (d_prop == d_bfs) {
            exact_matches++;
        } else {
            size_t r = i / k;
            size_t c = i % k;
            std::cerr << "  MISMATCH at pair (" << viewpoints[r].r << "," << viewpoints[r].c
                      << ") -> (" << viewpoints[c].r << "," << viewpoints[c].c << "): "
                      << "prop=" << d_prop << " vs bfs=" << d_bfs << "\n";
        }
    }

    std::cout << "  [" << test_name << "] Queries: " << total_queries
              << ", Matches: " << exact_matches << " / " << total_queries
              << ", MaxDiff: " << max_diff << "\n";

    return exact_matches == total_queries;
}

// Helper: incremental matrix refresh must be bit-identical to the full query.
bool verifyIncrementalMatches(const std::string& test_name,
                              const MinPlusOracle& oracle,
                              const GridMap& grid,
                              const std::vector<Point>& viewpoints,
                              MinPlusOracle::HTripQueryScratch& scratch,
                              MinPlusOracle::HTripIncrementalState& state) {
    std::vector<dist_t> inc_mat;
    std::vector<dist_t> full_mat;
    oracle.queryAllPairsIncremental(viewpoints, grid, inc_mat, scratch, state);
    oracle.queryAllPairs(viewpoints, grid, full_mat, scratch);

    size_t mismatch = 0;
    if (inc_mat.size() != full_mat.size()) {
        std::cerr << "  INCREMENTAL size mismatch: inc=" << inc_mat.size()
                  << " full=" << full_mat.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < inc_mat.size(); ++i) {
        if (inc_mat[i] != full_mat[i]) {
            if (mismatch < 5) {
                std::cerr << "  INCREMENTAL MISMATCH at " << i << ": inc=" << inc_mat[i]
                          << " full=" << full_mat[i] << "\n";
            }
            mismatch++;
        }
    }
    if (mismatch != 0) {
        std::cout << "  [" << test_name << "] incremental refresh FAILED ("
                  << mismatch << " mismatches)\n";
    }
    return mismatch == 0;
}

int main(int argc, char** argv) {
    int seed = 12345;
    if (argc > 1) {
        seed = std::atoi(argv[1]);
    }

    std::cout << "=========================================================\n";
    std::cout << "  ICRA 2027: H-TRIP Rigorous Exactness Verification Suite\n";
    std::cout << "=========================================================\n\n";

    // -------------------------------------------------------------
    // Suite 1: Deterministic Multi-Level Seam Crossing Regression (P1)
    // -------------------------------------------------------------
    std::cout << "[Suite 1: Deterministic Multi-Level Seam Crossing Regression]\n";
    {
        GridMap seam_grid(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                seam_grid.blockCorridor({r, c});
            }
        }
        // Open cell (4, 16) in leaf (2, 0)
        seam_grid.openDoorway({4, 16});

        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle seam_oracle(cfg);
        seam_oracle.build(seam_grid);

        // Dynamically open cell (4, 15) in leaf (1, 0)
        Point p1{4, 15};
        Point p2{4, 16};
        seam_grid.openDoorway(p1);
        seam_oracle.onCellOpened(p1, seam_grid);

        dist_t d_prop = seam_oracle.queryDistance(p1, p2, seam_grid);
        dist_t d_bfs = BfsOracle::queryDistance(seam_grid, p1, p2);

        std::cout << "  Seam (4,15)->(4,16): prop=" << d_prop << " vs bfs=" << d_bfs << "\n";
        if (d_prop != 1 || d_bfs != 1) {
            std::cerr << "FAIL: Suite 1 seam crossing failed! prop=" << d_prop << ", bfs=" << d_bfs << "\n";
            return 1;
        }

        // Open a corridor connecting to the rest of the map
        for (int c = 10; c <= 20; ++c) {
            Point p{4, c};
            seam_grid.openDoorway(p);
            seam_oracle.onCellOpened(p, seam_grid);
        }

        std::vector<Point> seam_viewpoints = {{4, 10}, {4, 13}, {4, 15}, {4, 16}, {4, 18}, {4, 20}};
        if (!verifyAllPairsMatrix("Multi-tile seam corridor", seam_oracle, seam_grid, seam_viewpoints)) {
            std::cerr << "FAIL: Suite 1 corridor all-pairs verification failed!\n";
            return 1;
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 1b: Deterministic Level-2+ Seam Crossing Regression (N3)
    // -------------------------------------------------------------
    std::cout << "[Suite 1b: Deterministic Level-2+ Seam Crossing Regression]\n";
    {
        // 32x32 grid with b=8, g=2.
        // Leaf size is 8x8. Level 1 nodes span 16x16. Level 2 node spans 32x32 (Root).
        // (15, 8) is in leaf (tx=1, ty=1) -> Level 1 node (0, 0).
        // (16, 8) is in leaf (tx=1, ty=2) -> Level 1 node (0, 1).
        // Opening (15,8)-(16,8) is a vertical seam crossing between two DIFFERENT level-1 nodes!
        // Its LCA level is exactly 2 (root level).
        GridMap seam_grid(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                seam_grid.blockCorridor({r, c});
            }
        }
        Point p2{16, 8};
        seam_grid.openDoorway(p2);

        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle seam_oracle(cfg);
        seam_oracle.build(seam_grid);

        Point p1{15, 8};
        seam_grid.openDoorway(p1);
        seam_oracle.onCellOpened(p1, seam_grid);

        dist_t d_prop = seam_oracle.queryDistance(p1, p2, seam_grid);
        dist_t d_bfs = BfsOracle::queryDistance(seam_grid, p1, p2);

        std::cout << "  Level-2 Seam (15,8)->(16,8): prop=" << d_prop << " vs bfs=" << d_bfs << "\n";
        if (d_prop != 1 || d_bfs != 1) {
            std::cerr << "FAIL: Suite 1b level-2 seam crossing failed! prop=" << d_prop << ", bfs=" << d_bfs << "\n";
            return 1;
        }

        // Extend vertical corridor across the level-2 boundary: rows 10 to 22 at c=8
        for (int r = 10; r <= 22; ++r) {
            Point p{r, 8};
            seam_grid.openDoorway(p);
            seam_oracle.onCellOpened(p, seam_grid);
        }

        std::vector<Point> test_vps = {{10, 8}, {12, 8}, {15, 8}, {16, 8}, {18, 8}, {22, 8}};
        if (!verifyAllPairsMatrix("Level-2 vertical seam corridor", seam_oracle, seam_grid, test_vps)) {
            std::cerr << "FAIL: Suite 1b corridor verification failed!\n";
            return 1;
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 1c: child_to_ext staleness regression (interior-seam shortcut)
    // -------------------------------------------------------------
    std::cout << "[Suite 1c: Interior-Seam Shortcut Regression (child_to_ext freshness)]\n";
    {
        // 32x32, b=8, g=2. Node (0,0) spans rows 0-15, cols 0-15.
        // A U-shaped corridor inside node (0,0) connects interior ports
        // (8,8) -> (8,13) -> (14,13) -> (14,8) and touches the exterior
        // ring at (15,8) only. Opening column-8 cells shortcuts
        // (8,8) <-> (14,8) from 16 down to 6 WITHOUT changing any
        // exterior-exterior summary S_U of node (0,0).
        // A stale child_to_ext cache would corrupt the lifted distance
        // from (8,9) to exterior port (15,8) and produce wrong answers.
        GridMap g(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                g.blockCorridor({r, c});
            }
        }
        // Exterior ring of node (0,0): rows 0/15, cols 0/15 (+ row 15 eastward)
        for (int c = 0; c <= 15; ++c) g.openDoorway({0, c});
        for (int c = 0; c <= 31; ++c) g.openDoorway({15, c});
        for (int r = 1; r <= 14; ++r) g.openDoorway({r, 0});
        for (int r = 1; r <= 14; ++r) g.openDoorway({r, 15});
        // U corridor inside node (0,0)
        for (int c = 8; c <= 13; ++c) g.openDoorway({8, c});
        for (int r = 9; r <= 14; ++r) g.openDoorway({r, 13});
        for (int c = 8; c <= 12; ++c) g.openDoorway({14, c});
        // East corridor + target point
        for (int r = 15; r <= 20; ++r) g.openDoorway({r, 20});
        g.openDoorway({20, 20});

        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        std::vector<Point> vps = {{8, 9}, {14, 9}, {9, 13}, {20, 20}, {15, 10}};
        if (!verifyAllPairsMatrix("Suite 1c initial", oracle, g, vps)) {
            std::cerr << "FAIL: Suite 1c initial verification failed!\n";
            return 1;
        }

        // Open the interior column-8 shortcut cells one at a time
        const std::vector<Point> shortcut = {{9, 8}, {10, 8}, {11, 8}, {12, 8}};
        for (const Point& p : shortcut) {
            g.openDoorway(p);
            oracle.onCellOpened(p, g);
            std::string label = "Suite 1c after (" + std::to_string(p.r) + "," + std::to_string(p.c) + ")";
            if (!verifyAllPairsMatrix(label, oracle, g, vps)) {
                std::cerr << "FAIL: Suite 1c stale child_to_ext regression!\n";
                return 1;
            }
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 2: Multi-Branching 16-ary Hierarchy (g = 4) Exactness
    // -------------------------------------------------------------
    std::cout << "[Suite 2: Multi-Branching 16-ary Hierarchy (g = 4) Exactness]\n";
    {
        GridMap g4_grid(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                g4_grid.openDoorway({r, c});
            }
        }

        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 4;
        cfg.group_h = 4; // 16-ary branching
        HTripOracle g4_oracle(cfg);
        g4_oracle.build(g4_grid);

        int total_q = 0;
        int matches = 0;
        for (int r1 = 2; r1 < 30; r1 += 7) {
            for (int c1 = 2; c1 < 30; c1 += 7) {
                Point p1{r1, c1};
                for (int r2 = 3; r2 < 30; r2 += 7) {
                    for (int c2 = 3; c2 < 30; c2 += 7) {
                        Point p2{r2, c2};
                        dist_t d_prop = g4_oracle.queryDistance(p1, p2, g4_grid);
                        dist_t d_bfs = BfsOracle::queryDistance(g4_grid, p1, p2);
                        total_q++;
                        if (d_prop == d_bfs) matches++;
                        else {
                            std::cerr << "FAIL: Suite 2 mismatch for (" << r1 << "," << c1
                                      << ") -> (" << r2 << "," << c2 << "): prop="
                                      << d_prop << " vs bfs=" << d_bfs << "\n";
                            return 1;
                        }
                    }
                }
            }
        }
        std::cout << "  [16-ary Tree] Tested " << total_q << " queries | Matches: "
                  << matches << " / " << total_q << "\n";
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 3: Dynamic Root Expansion (liftRoot) (B1 & N1)
    // -------------------------------------------------------------
    std::cout << "[Suite 3: Dynamic Root Expansion (liftRoot)]\n";
    {
        // 1. Standard expansion: 16x16 (H=1) -> 32x32 (H=2)
        GridMap grid16(16, 16);
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                grid16.openDoorway({r, c});
            }
        }

        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle lift_oracle(cfg);
        lift_oracle.build(grid16);

        // Expand to 32x32
        GridMap grid32(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                grid32.openDoorway({r, c});
            }
        }
        lift_oracle.liftRoot(grid32);

        std::vector<Point> test_pts = {{2, 2}, {5, 5}, {2, 20}, {20, 20}, {25, 10}};
        if (!verifyAllPairsMatrix("Lifted root expanded world", lift_oracle, grid32, test_pts)) {
            std::cerr << "FAIL: Suite 3 liftRoot all-pairs mismatch!\n";
            return 1;
        }

        // 2. N1 regression test: liftRoot when old_H == 0 (single-leaf initial oracle)
        GridMap grid8(8, 8);
        for (int r = 0; r < 8; ++r) {
            for (int c = 0; c < 8; ++c) {
                grid8.openDoorway({r, c});
            }
        }
        HierarchyConfig cfg0;
        cfg0.tile_size_b = 8;
        cfg0.group_w = 2;
        cfg0.group_h = 2;
        MinPlusOracle lift_oracle0(cfg0);
        lift_oracle0.build(grid8);
        assert(lift_oracle0.depth() == 0); // old_H == 0!

        // Expand 8x8 to 16x16 (liftRoot from old_H == 0)
        lift_oracle0.liftRoot(grid16);
        assert(lift_oracle0.depth() == 1);
        std::vector<Point> test_pts0 = {{1, 1}, {3, 3}, {6, 6}, {10, 10}, {14, 14}};
        if (!verifyAllPairsMatrix("Lifted root from old_H=0", lift_oracle0, grid16, test_pts0)) {
            std::cerr << "FAIL: Suite 3 liftRoot from old_H=0 mismatch!\n";
            return 1;
        }

        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Large Grid Shared Setup for Suites 4, 5, 6
    // -------------------------------------------------------------
    int large_size = 128;
    int tile_b = 16;
    GridMap large_grid;
    large_grid.generateProceduralRooms(large_size, large_size, 16, seed);

    // Exploration snapshot with 60% observed
    Point start{2, 2};
    large_grid.initExplorationSnapshot(start, 0.60, 0.15, seed);

    FrontierManager fm(2);
    auto clusters = fm.detectAndCluster(large_grid);
    std::vector<Point> viewpoints;
    for (const auto& cl : clusters) {
        if (large_grid.isPassable(cl.representative)) {
            viewpoints.push_back(cl.representative);
        }
    }

    std::mt19937_64 vp_rng(seed + 42);
    std::vector<Point> all_free_cells;
    for (int r = 0; r < large_size; ++r) {
        for (int c = 0; c < large_size; ++c) {
            Point pt{r, c};
            if (large_grid.isPassable(pt)) {
                all_free_cells.push_back(pt);
            }
        }
    }
    while (viewpoints.size() < 64 && !all_free_cells.empty()) {
        size_t idx = vp_rng() % all_free_cells.size();
        viewpoints.push_back(all_free_cells[idx]);
    }
    if (viewpoints.size() > 64) {
        viewpoints.resize(64);
    }

    HierarchyConfig cfg_large;
    cfg_large.tile_size_b = tile_b;
    cfg_large.group_w = 2;
    cfg_large.group_h = 2;

    HTripOracle large_oracle(cfg_large);
    large_oracle.build(large_grid);

    // -------------------------------------------------------------
    // Suite 4: Large Grid from Scratch with Thousands of Pairwise Queries
    // -------------------------------------------------------------
    std::cout << "[Suite 4: Large Grid All-Pairs Exactness from Scratch]\n";
    {
        std::cout << "  Testing " << viewpoints.size() << " viewpoints ("
                  << (viewpoints.size() * viewpoints.size()) << " all-pairs queries)...\n";
        std::cout << "  Oracle depth H: " << large_oracle.depth()
                  << ", Total memory: " << (large_oracle.totalMemoryBytes() / 1024) << " KB\n";

        if (!verifyAllPairsMatrix("Large Grid 4,096 Queries", large_oracle, large_grid, viewpoints)) {
            std::cerr << "FAIL: Suite 4 large grid verification failed!\n";
            return 1;
        }
        std::cout << "  STATUS: PASSED (100% bit-exact across thousands of queries)\n\n";
    }

    // -------------------------------------------------------------
    // Suite 5: Dynamic Edge / Cell Revelations (Conway Rank-1 Updates)
    // -------------------------------------------------------------
    std::cout << "[Suite 5: Dynamic Doorway Revelations & Conway Rank-1 Perturbations]\n";
    {
        std::mt19937_64 rng(seed + 99);
        int dynamic_openings = 40;
        int opened_count = 0;

        for (int step = 0; step < dynamic_openings; ++step) {
            int r = static_cast<int>(rng() % static_cast<uint64_t>(large_size));
            int c = static_cast<int>(rng() % static_cast<uint64_t>(large_size));
            Point pt{r, c};

            if (large_grid.ground_truth[large_grid.index(pt)] == CellState::Free &&
                large_grid.getObserved(pt) != CellState::Free) {

                large_grid.openDoorway(pt);
                large_oracle.onCellOpened(pt, large_grid);
                opened_count++;
                all_free_cells.push_back(pt);

                // Verify all-pairs matrix periodically
                if (opened_count % 10 == 0) {
                    std::string label = "Dynamic step " + std::to_string(opened_count);
                    if (!verifyAllPairsMatrix(label, large_oracle, large_grid, viewpoints)) {
                        std::cerr << "FAIL: Suite 5 mismatch after opening " << opened_count << " doorways!\n";
                        return 1;
                    }
                }
            }
        }
        std::cout << "  Successfully tested " << opened_count << " dynamic Conway star perturbations.\n";
        std::cout << "  STATUS: PASSED (100% bit-exact across all dynamic openings)\n\n";
    }

    // -------------------------------------------------------------
    // Suite 6: Dynamic Obstacle Blockages
    // -------------------------------------------------------------
    std::cout << "[Suite 6: Dynamic Corridor Blockages & Rebuilds]\n";
    {
        int blocked_count = 0;
        for (int step = 0; step < 20; ++step) {
            if (all_free_cells.empty()) break;
            size_t idx = vp_rng() % all_free_cells.size();
            Point pt = all_free_cells[idx];

            if (large_grid.getObserved(pt) == CellState::Free) {
                large_grid.blockCorridor(pt);
                large_oracle.onCellBlocked(pt, large_grid);
                blocked_count++;

                if (blocked_count % 5 == 0) {
                    std::string label = "Blockage step " + std::to_string(blocked_count);
                    if (!verifyAllPairsMatrix(label, large_oracle, large_grid, viewpoints)) {
                        std::cerr << "FAIL: Suite 6 mismatch after blocking " << blocked_count << " corridors!\n";
                        return 1;
                    }
                }
            }
        }
        std::cout << "  Successfully tested " << blocked_count << " dynamic corridor blockages.\n";
        std::cout << "  STATUS: PASSED (100% bit-exact across all dynamic blockages)\n\n";
    }

    // -------------------------------------------------------------
    // Suite 7: Version-Stamped Incremental Matrix Refresh (exact vs full)
    // -------------------------------------------------------------
    std::cout << "[Suite 7: Version-Stamped Incremental Matrix Refresh]\n";
    {
        GridMap inc_grid;
        inc_grid.generateProceduralRooms(128, 128, 16, seed + 777);
        inc_grid.initExplorationSnapshot({2, 2}, 0.60, 0.15, seed + 7);

        std::vector<Point> inc_free;
        for (int r = 0; r < 128; ++r) {
            for (int c = 0; c < 128; ++c) {
                if (inc_grid.getObserved(r, c) == CellState::Free) {
                    inc_free.push_back({r, c});
                }
            }
        }
        std::mt19937_64 rng7(seed + 1234);
        std::shuffle(inc_free.begin(), inc_free.end(), rng7);
        std::vector<Point> inc_vps;
        for (size_t i = 0; i < 48 && i < inc_free.size(); ++i) {
            inc_vps.push_back(inc_free[i]);
        }

        HierarchyConfig cfg7;
        cfg7.tile_size_b = 16;
        cfg7.group_w = 2;
        cfg7.group_h = 2;
        HTripOracle oracle7(cfg7);
        oracle7.build(inc_grid);

        HTripOracle::HTripQueryScratch scratch7;
        HTripOracle::HTripIncrementalState state7;
        std::vector<dist_t> inc_mat;
        std::vector<dist_t> full_mat;

        auto matrices_equal = [](const std::vector<dist_t>& a, const std::vector<dist_t>& b) {
            if (a.size() != b.size()) return false;
            for (size_t t = 0; t < a.size(); ++t) {
                if (a[t] != b[t]) return false;
            }
            return true;
        };

        // Initial incremental call (internally forces a full recomputation)
        oracle7.queryAllPairsIncremental(inc_vps, inc_grid, inc_mat, scratch7, state7);
        full_mat = oracle7.queryAllPairs(inc_vps, inc_grid);
        if (!matrices_equal(inc_mat, full_mat)) {
            std::cerr << "FAIL: Suite 7 initial incremental matrix differs from full!\n";
            return 1;
        }

        // Dynamic openings: incremental refresh must stay bit-exact against the
        // full recomputation and against the BFS oracle after every opening.
        std::mt19937_64 rng8(seed + 99);
        int opened = 0;
        for (int step = 0; step < 60 && opened < 25; ++step) {
            const int r = static_cast<int>(rng8() % 128);
            const int c = static_cast<int>(rng8() % 128);
            const Point p{r, c};
            if (inc_grid.ground_truth[inc_grid.index(p)] != CellState::Free) continue;
            if (inc_grid.getObserved(p) == CellState::Free) continue;

            inc_grid.openDoorway(p);
            oracle7.onCellOpened(p, inc_grid);
            opened++;

            oracle7.queryAllPairsIncremental(inc_vps, inc_grid, inc_mat, scratch7, state7);

            // Hazard check: sharing the same scratch with an unrelated full query
            // (different viewpoint list!) must not corrupt the incremental
            // state's cached lifts.
            if (opened % 5 == 0) {
                std::vector<dist_t> throwaway;
                std::vector<Point> other_vps(inc_vps.rbegin(), inc_vps.rbegin() + 10);
                oracle7.queryAllPairs(other_vps, inc_grid, throwaway, scratch7);
            }

            full_mat = oracle7.queryAllPairs(inc_vps, inc_grid);
            if (!matrices_equal(inc_mat, full_mat)) {
                std::cerr << "FAIL: Suite 7 incremental vs full mismatch after opening " << opened << "\n";
                return 1;
            }
        }

        auto bfs_mat = BfsOracle::queryAllPairs(inc_grid, inc_vps);
        if (!matrices_equal(inc_mat, bfs_mat)) {
            std::cerr << "FAIL: Suite 7 incremental matrix differs from BFS oracle!\n";
            return 1;
        }

        std::cout << "  Incremental refresh bit-exact over " << opened << " dynamic openings.\n";
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 8: Dense AVX2 Meet Tail Regression (>=32-port target)
    // -------------------------------------------------------------
    std::cout << "[Suite 8: Dense Meet Tail Regression (AVX2 port tail)]\n";
    {
        // 32x32, b=16, g=2 (H=1). Leaf (0,1) (rows 16-31, cols 0-15) is fully
        // free and contains the target T=(25,10): all 60 perimeter ports are
        // active (>= 32 -> dense AVX2 meet path). Its ONLY connection to the
        // rest of the map is the seam (20,15)-(20,16); port (20,15) is a
        // right-column port with index 49, i.e. inside the 16-lane tail
        // [48..59]. The northern (checked) ports lead through a long detour,
        // so the true optimum uses the tail port. A dense meet that drops the
        // scalar tail overestimates d((2,24),(25,10)) (true value 37).
        GridMap g(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                g.blockCorridor({r, c});
            }
        }
        // Source-side corridor in leaf (1,0): row 2 and column 24
        for (int c = 16; c <= 31; ++c) g.openDoorway({2, c});
        for (int r = 2; r <= 15; ++r) g.openDoorway({r, 24});
        // Highway through leaf (1,1): column 24 (rows 16..20) and row 20
        for (int r = 16; r <= 20; ++r) g.openDoorway({r, 24});
        for (int c = 16; c <= 31; ++c) g.openDoorway({20, c});
        // Fully open target leaf (0,1)
        for (int r = 16; r <= 31; ++r) {
            for (int c = 0; c <= 15; ++c) {
                g.openDoorway({r, c});
            }
        }

        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        const Point S{2, 24};
        const Point T{25, 10};
        std::vector<Point> vps = {S, T, {20, 20}, {8, 24}};

        dist_t d_bfs = BfsOracle::queryDistance(g, S, T);
        dist_t d_point = oracle.queryDistance(S, T, g);
        std::cout << "  Tail-port pair (2,24)->(25,10): bfs=" << d_bfs
                  << " point=" << d_point << "\n";
        if (d_bfs != 37) {
            std::cerr << "FAIL: Suite 8 scenario drift, expected bfs=37 got " << d_bfs << "\n";
            return 1;
        }
        if (d_point != d_bfs) {
            std::cerr << "FAIL: Suite 8 point query tail-port regression!\n";
            return 1;
        }
        if (!verifyAllPairsMatrix("Suite 8 dense meet tail", oracle, g, vps)) {
            std::cerr << "FAIL: Suite 8 dense AVX2 meet tail regression!\n";
            return 1;
        }

        // Incremental path must equal the full recomputation on this geometry too
        MinPlusOracle::HTripQueryScratch scratch8;
        MinPlusOracle::HTripIncrementalState state8;
        std::vector<dist_t> inc8, full8;
        oracle.queryAllPairsIncremental(vps, g, inc8, scratch8, state8);
        full8 = oracle.queryAllPairs(vps, g);
        for (size_t t = 0; t < full8.size(); ++t) {
            if (inc8[t] != full8[t]) {
                std::cerr << "FAIL: Suite 8 incremental mismatch on tail geometry!\n";
                return 1;
            }
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 9: Incremental Edge Cases (hybrid fallback, reopened viewpoint)
    // -------------------------------------------------------------
    std::cout << "[Suite 9: Incremental Edge Cases (hybrid fallback, reopened viewpoint)]\n";
    {
        // (a) Hybrid fallback must be safe for the incremental API, even with a
        // fresh scratch that the full query never populated.
        GridMap hg(32, 32);
        hg.generateProceduralRooms(32, 32, 16, seed + 55);
        hg.initExplorationSnapshot({1, 1}, 0.60, 0.10, seed + 6);
        HierarchyConfig hcfg;
        hcfg.tile_size_b = 16;
        hcfg.group_w = 2;
        hcfg.group_h = 2;
        hcfg.hybrid_k_threshold = 4; // k = 3 below threshold -> BFS fallback
        MinPlusOracle h_oracle(hcfg);
        h_oracle.build(hg);

        std::vector<Point> hvps;
        for (int r = 0; r < 32 && hvps.size() < 3; ++r) {
            for (int c = 0; c < 32 && hvps.size() < 3; ++c) {
                if (hg.getObserved(r, c) == CellState::Free) hvps.push_back({r, c});
            }
        }
        MinPlusOracle::HTripQueryScratch h_scratch;
        MinPlusOracle::HTripIncrementalState h_state;
        std::vector<dist_t> h_mat;
        h_oracle.queryAllPairsIncremental(hvps, hg, h_mat, h_scratch, h_state);
        std::vector<dist_t> h_ref = BfsOracle::queryAllPairs(hg, hvps);
        for (size_t t = 0; t < h_ref.size(); ++t) {
            if (h_mat[t] != h_ref[t]) {
                std::cerr << "FAIL: Suite 9a hybrid incremental mismatch at " << t << "\n";
                return 1;
            }
        }
        h_oracle.queryAllPairsIncremental(hvps, hg, h_mat, h_scratch, h_state);
        for (size_t t = 0; t < h_ref.size(); ++t) {
            if (h_mat[t] != h_ref[t]) {
                std::cerr << "FAIL: Suite 9a hybrid incremental mismatch on repeat at " << t << "\n";
                return 1;
            }
        }

        // (b) A viewpoint blocked during a reset and later reopened must regain a
        // zero self-distance in the incremental matrix.
        GridMap rg(32, 32);
        rg.generateProceduralRooms(32, 32, 16, seed + 77);
        rg.initExplorationSnapshot({1, 1}, 0.70, 0.10, seed + 8);
        HierarchyConfig rcfg;
        rcfg.tile_size_b = 16;
        rcfg.group_w = 2;
        rcfg.group_h = 2;
        MinPlusOracle r_oracle(rcfg);
        r_oracle.build(rg);

        std::vector<Point> rvps;
        for (int r = 0; r < 32 && rvps.size() < 4; ++r) {
            for (int c = 0; c < 32 && rvps.size() < 4; ++c) {
                if (rg.getObserved(r, c) == CellState::Free) rvps.push_back({r, c});
            }
        }
        MinPlusOracle::HTripQueryScratch r_scratch;
        MinPlusOracle::HTripIncrementalState r_state;
        std::vector<dist_t> r_mat;
        r_oracle.queryAllPairsIncremental(rvps, rg, r_mat, r_scratch, r_state);

        const size_t blocked_idx = 1;
        rg.blockCorridor(rvps[blocked_idx]);
        r_oracle.onCellBlocked(rvps[blocked_idx], rg);
        r_oracle.queryAllPairsIncremental(rvps, rg, r_mat, r_scratch, r_state); // reset while blocked

        rg.openDoorway(rvps[blocked_idx]);
        r_oracle.onCellOpened(rvps[blocked_idx], rg);
        r_oracle.queryAllPairsIncremental(rvps, rg, r_mat, r_scratch, r_state);

        const size_t k9 = rvps.size();
        if (r_mat[blocked_idx * k9 + blocked_idx] != 0) {
            std::cerr << "FAIL: Suite 9b reopened viewpoint diagonal is "
                      << r_mat[blocked_idx * k9 + blocked_idx] << " (expected 0)\n";
            return 1;
        }
        std::vector<dist_t> r_full = r_oracle.queryAllPairs(rvps, rg);
        for (size_t t = 0; t < r_full.size(); ++t) {
            if (r_mat[t] != r_full[t]) {
                std::cerr << "FAIL: Suite 9b incremental vs full mismatch at " << t << "\n";
                return 1;
            }
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 10: Single-Leaf Hierarchy (H = 0) Dynamic Exactness
    // -------------------------------------------------------------
    std::cout << "[Suite 10: Single-Leaf Hierarchy (H=0) Dynamic Exactness]\n";
    {
        GridMap g(16, 16);
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                g.blockCorridor({r, c});
            }
        }
        for (int c = 2; c <= 13; ++c) g.openDoorway({8, c});

        HierarchyConfig cfg;
        cfg.tile_size_b = 16; // S == b => H == 0 (single leaf, no internal levels)
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        std::vector<Point> vps = {{8, 2}, {8, 6}, {8, 10}, {8, 13}};
        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        if (!verifyAllPairsMatrix("Suite 10 initial", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 10 initial", oracle, g, vps, scratch, state)) return 1;

        const std::vector<Point> branches = {{7, 6}, {6, 6}, {5, 6}, {4, 6},
                                             {9, 10}, {10, 10}, {11, 10}, {12, 10}};
        for (const Point& p : branches) {
            g.openDoorway(p);
            oracle.onCellOpened(p, g);
        }
        if (!verifyAllPairsMatrix("Suite 10 after branches", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 10 after branches", oracle, g, vps, scratch, state)) return 1;
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 11: Non-Square Workspace Exactness
    // -------------------------------------------------------------
    std::cout << "[Suite 11: Non-Square Workspace (20x44) Exactness]\n";
    {
        GridMap g;
        g.generateProceduralRooms(20, 44, 8, seed + 411);
        g.initExplorationSnapshot({2, 2}, 0.65, 0.10, seed + 12);

        std::vector<Point> vps;
        for (int r = 0; r < 20 && vps.size() < 10; ++r) {
            for (int c = 0; c < 44 && vps.size() < 10; ++c) {
                if (g.getObserved(r, c) == CellState::Free) vps.push_back({r, c});
            }
        }
        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2; // S pads to 64 => H = 3
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        if (!verifyAllPairsMatrix("Suite 11 initial", oracle, g, vps)) return 1;
        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        if (!verifyIncrementalMatches("Suite 11 initial", oracle, g, vps, scratch, state)) return 1;

        std::mt19937_64 rng(seed + 13);
        int opened = 0;
        for (int step = 0; step < 120 && opened < 15; ++step) {
            const int r = static_cast<int>(rng() % 20);
            const int c = static_cast<int>(rng() % 44);
            const Point p{r, c};
            if (g.ground_truth[g.index(p)] != CellState::Free) continue;
            if (g.getObserved(p) == CellState::Free) continue;
            g.openDoorway(p);
            oracle.onCellOpened(p, g);
            opened++;
        }
        if (!verifyAllPairsMatrix("Suite 11 dynamic", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 11 dynamic", oracle, g, vps, scratch, state)) return 1;
        std::cout << "  STATUS: PASSED\n";
        std::cout << "  Openings applied: " << opened << "\n\n";
    }

    // -------------------------------------------------------------
    // Suite 12: Disconnected Components and Non-Passable Viewpoint
    // -------------------------------------------------------------
    std::cout << "[Suite 12: Disconnected Halves and Blocked Viewpoint]\n";
    {
        GridMap g(32, 32);
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                g.blockCorridor({r, c});
            }
        }
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                if (c != 16) g.openDoorway({r, c});
            }
        }
        const Point blocked_viewpoint{20, 16}; // stays impassable for the whole suite
        std::vector<Point> vps = {{4, 4}, {8, 12}, {4, 20}, {8, 28}, blocked_viewpoint};

        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        if (!verifyAllPairsMatrix("Suite 12 separated", oracle, g, vps)) return 1;
        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        if (!verifyIncrementalMatches("Suite 12 separated", oracle, g, vps, scratch, state)) return 1;

        // Blocked viewpoint must stay INF on its whole row and diagonal
        auto separated = oracle.queryAllPairs(vps, g);
        const size_t k = vps.size();
        if (separated[4 * k + 4] != INF) {
            std::cerr << "FAIL: Suite 12 blocked viewpoint diagonal is "
                      << separated[4 * k + 4] << " (expected INF)\n";
            return 1;
        }

        // Open a door in the wall and re-verify (cross pairs become finite)
        g.openDoorway({8, 16});
        oracle.onCellOpened({8, 16}, g);
        if (!verifyAllPairsMatrix("Suite 12 door opened", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 12 door opened", oracle, g, vps, scratch, state)) return 1;
        auto connected = oracle.queryAllPairs(vps, g);
        if (connected[0 * k + 2] >= INF) {
            std::cerr << "FAIL: Suite 12 cross-half distance stayed INF after opening the door\n";
            return 1;
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 13: 16-ary (g=4) Incremental Dynamic Exactness
    // -------------------------------------------------------------
    std::cout << "[Suite 13: 16-ary (g=4) Incremental Dynamic Exactness]\n";
    {
        GridMap g;
        g.generateProceduralRooms(64, 64, 16, seed + 613);
        g.initExplorationSnapshot({2, 2}, 0.62, 0.10, seed + 14);

        std::vector<Point> vps;
        for (int r = 0; r < 64 && vps.size() < 20; ++r) {
            for (int c = 0; c < 64 && vps.size() < 20; ++c) {
                if (g.getObserved(r, c) == CellState::Free) vps.push_back({r, c});
            }
        }
        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 4;
        cfg.group_h = 4; // S == 64 => H = 2
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        if (!verifyAllPairsMatrix("Suite 13 initial", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 13 initial", oracle, g, vps, scratch, state)) return 1;

        std::mt19937_64 rng(seed + 15);
        int opened = 0;
        for (int step = 0; step < 200 && opened < 20; ++step) {
            const int r = static_cast<int>(rng() % 64);
            const int c = static_cast<int>(rng() % 64);
            const Point p{r, c};
            if (g.ground_truth[g.index(p)] != CellState::Free) continue;
            if (g.getObserved(p) == CellState::Free) continue;
            g.openDoorway(p);
            oracle.onCellOpened(p, g);
            opened++;
            if (opened % 5 == 0) {
                if (!verifyAllPairsMatrix("Suite 13 step", oracle, g, vps)) return 1;
                if (!verifyIncrementalMatches("Suite 13 step", oracle, g, vps, scratch, state)) return 1;
            }
        }
        if (!verifyAllPairsMatrix("Suite 13 final", oracle, g, vps)) return 1;
        if (!verifyIncrementalMatches("Suite 13 final", oracle, g, vps, scratch, state)) return 1;
        std::cout << "  STATUS: PASSED\n";
        std::cout << "  Openings applied: " << opened << "\n\n";
    }

    // -------------------------------------------------------------
    // Suite 14: Viewpoint List Changes (k and order) Reset State
    // -------------------------------------------------------------
    std::cout << "[Suite 14: Viewpoint List Changes Across Incremental Calls]\n";
    {
        GridMap g;
        g.generateProceduralRooms(32, 32, 16, seed + 814);
        g.initExplorationSnapshot({1, 1}, 0.70, 0.10, seed + 16);

        std::vector<Point> free_cells;
        for (int r = 0; r < 32; ++r) {
            for (int c = 0; c < 32; ++c) {
                if (g.getObserved(r, c) == CellState::Free) free_cells.push_back({r, c});
            }
        }
        std::mt19937_64 rng(seed + 17);
        std::shuffle(free_cells.begin(), free_cells.end(), rng);

        std::vector<Point> list_a(free_cells.begin(), free_cells.begin() + 5);
        std::vector<Point> list_b(free_cells.begin(), free_cells.begin() + 8);
        std::reverse(list_b.begin(), list_b.end());

        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        if (!verifyIncrementalMatches("Suite 14 list A", oracle, g, list_a, scratch, state)) return 1;
        if (!verifyIncrementalMatches("Suite 14 list B", oracle, g, list_b, scratch, state)) return 1;
        if (!verifyIncrementalMatches("Suite 14 list A again", oracle, g, list_a, scratch, state)) return 1;

        // Same list, different order: must also reset and stay exact
        std::vector<Point> list_a_shuffled = list_a;
        std::reverse(list_a_shuffled.begin(), list_a_shuffled.end());
        if (!verifyIncrementalMatches("Suite 14 list A reordered", oracle, g, list_a_shuffled, scratch, state)) return 1;
        std::cout << "  STATUS: PASSED\n\n";
    }

    // -------------------------------------------------------------
    // Suite 15: queryDistance vs Batched Matrix Consistency
    // -------------------------------------------------------------
    std::cout << "[Suite 15: Single-Pair Query vs Batched Matrix Consistency]\n";
    {
        GridMap g;
        g.generateProceduralRooms(40, 40, 16, seed + 915);
        g.initExplorationSnapshot({2, 2}, 0.68, 0.10, seed + 18);

        std::vector<Point> vps;
        for (int r = 0; r < 40 && vps.size() < 10; ++r) {
            for (int c = 0; c < 40 && vps.size() < 10; ++c) {
                if (g.getObserved(r, c) == CellState::Free) vps.push_back({r, c});
            }
        }
        HierarchyConfig cfg;
        cfg.tile_size_b = 8;
        cfg.group_w = 2;
        cfg.group_h = 2;
        MinPlusOracle oracle(cfg);
        oracle.build(g);

        auto mat = oracle.queryAllPairs(vps, g);
        const size_t k = vps.size();
        size_t checked = 0;
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                const dist_t single = oracle.queryDistance(vps[i], vps[j], g);
                if (single != mat[i * k + j]) {
                    std::cerr << "FAIL: Suite 15 queryDistance(" << i << "," << j << ")=" << single
                              << " but matrix=" << mat[i * k + j] << "\n";
                    return 1;
                }
                checked++;
            }
        }
        std::cout << "  STATUS: PASSED (" << checked << " single-pair queries matched the matrix)\n\n";
    }

    // -------------------------------------------------------------
    // Suite 16: Batched Cell Openings & Hierarchy Boundary Exactness
    // (Multi-Tile, Multi-Level & H=0 Leaf Guard)
    // -------------------------------------------------------------
    std::cout << "[Suite 16: Batched Cell Openings (Multi-Tile & H=0 Guard)]\n";
    {
        // Sub-test 16.1: H=0 single tile (no internal hierarchy)
        {
            GridMap g(16, 16);
            for (int r = 0; r < 16; ++r) {
                for (int c = 0; c < 16; ++c) {
                    g.blockCorridor({r, c});
                }
            }
            // Open initial pathway
            for (int i = 0; i < 16; ++i) {
                g.openDoorway({i, 0});
            }

            HierarchyConfig cfg;
            cfg.tile_size_b = 16;
            cfg.group_w = 2;
            cfg.group_h = 2;
            MinPlusOracle oracle(cfg);
            oracle.build(g);

            std::vector<Point> vps = {{0, 0}, {5, 0}, {10, 0}, {15, 0}, {5, 10}, {10, 10}};
            std::vector<Point> batch;
            for (int c = 1; c <= 10; ++c) {
                batch.push_back({5, c});
                batch.push_back({10, c});
            }
            for (Point p : batch) {
                g.openDoorway(p);
            }
            oracle.onBatchCellsOpened(batch, g);

            auto mat = oracle.queryAllPairs(vps, g);
            auto bfs_mat = BfsOracle::queryAllPairs(g, vps);
            if (mat != bfs_mat) {
                std::cerr << "FAIL: Suite 16.1 H=0 batch update differs from BFS oracle!\n";
                return 1;
            }
            std::cout << "  Sub-test 16.1 (H=0 leaf-only batch): PASSED\n";
        }

        // Sub-test 16.2: Multi-tile batch updates with both small and large batches
        {
            GridMap g;
            g.generateProceduralRooms(64, 64, 16, seed + 916);
            // Snapshot with ~50% unknown
            g.initExplorationSnapshot({4, 4}, 0.50, 0.15, seed + 33);

            HierarchyConfig cfg;
            cfg.tile_size_b = 8;
            cfg.group_w = 2;
            cfg.group_h = 2;
            MinPlusOracle oracle(cfg);
            oracle.build(g);

            std::vector<Point> vps;
            for (int r = 4; r < 60 && vps.size() < 12; r += 5) {
                for (int c = 4; c < 60 && vps.size() < 12; c += 5) {
                    if (g.getObserved(r, c) == CellState::Free) {
                        vps.push_back({r, c});
                    }
                }
            }

            std::mt19937_64 rng(seed + 42);
            for (int round = 0; round < 10; ++round) {
                std::vector<Point> batch;
                int batch_target = (round % 2 == 0) ? 8 : 28; // test both small SIMD loop and large fallback
                for (int tries = 0; tries < 200 && static_cast<int>(batch.size()) < batch_target; ++tries) {
                    int r = static_cast<int>(rng() % 64);
                    int c = static_cast<int>(rng() % 64);
                    Point p{r, c};
                    if (g.ground_truth[g.index(p)] == CellState::Free && g.getObserved(p) != CellState::Free) {
                        batch.push_back(p);
                        g.openDoorway(p);
                    }
                }

                if (!batch.empty()) {
                    oracle.onBatchCellsOpened(batch, g);
                    auto mat = oracle.queryAllPairs(vps, g);
                    auto bfs_mat = BfsOracle::queryAllPairs(g, vps);
                    if (mat != bfs_mat) {
                        std::cerr << "FAIL: Suite 16.2 round " << round << " (batch size "
                                  << batch.size() << ") differs from BFS oracle!\n";
                        for (size_t i = 0; i < vps.size(); ++i) {
                            for (size_t j = 0; j < vps.size(); ++j) {
                                size_t idx = i * vps.size() + j;
                                if (mat[idx] != bfs_mat[idx]) {
                                    std::cerr << "  vp[" << i << "](" << vps[i].r << "," << vps[i].c << ") -> vp["
                                              << j << "](" << vps[j].r << "," << vps[j].c << "): HTRIP="
                                              << mat[idx] << " BFS=" << bfs_mat[idx] << "\n";
                                }
                            }
                        }
                        return 1;
                    }
                }
            }
            std::cout << "  Sub-test 16.2 (Multi-tile batch Conway): PASSED\n";
        }

        // Sub-test 16.3: Large single-tile batch (> 16 cells in one tile, K > 24) on H > 0 hierarchy
        {
            GridMap g(64, 64);
            for (int r = 0; r < 64; ++r) {
                for (int c = 0; c < 64; ++c) {
                    g.setObserved(r, c, (c == 31) ? CellState::Obstacle : CellState::Free);
                }
            }
            HierarchyConfig cfg;
            cfg.tile_size_b = 16;
            cfg.group_w = 2;
            cfg.group_h = 2;
            MinPlusOracle oracle(cfg);
            oracle.build(g);

            std::vector<Point> batch;
            for (int r = 0; r < 64; ++r) {
                g.setObserved(r, 31, CellState::Free);
                batch.push_back({r, 31});
            }
            oracle.onBatchCellsOpened(batch, g);

            Point p1{32, 10};
            Point p2{32, 45};
            dist_t d_prop = oracle.queryDistance(p1, p2, g);
            dist_t d_bfs = BfsOracle::queryDistance(g, p1, p2);
            if (d_prop != d_bfs || d_prop == INF) {
                std::cerr << "FAIL: Suite 16.3 large batch mismatch: prop=" << d_prop << " vs bfs=" << d_bfs << "\n";
                return 1;
            }
            std::cout << "  Sub-test 16.3 (Large single-tile batch >16 cells): PASSED\n";
        }

        // Sub-test 16.4: liftRoot after batch openings (verifying root S_U freshness)
        {
            GridMap g(64, 64);
            for (int r = 0; r < 64; ++r) {
                for (int c = 0; c < 64; ++c) {
                    g.setObserved(r, c, (r == 0) ? CellState::Obstacle : CellState::Free);
                }
            }
            HierarchyConfig cfg;
            cfg.tile_size_b = 8;
            cfg.group_w = 2;
            cfg.group_h = 2;
            MinPlusOracle oracle(cfg);
            oracle.build(g);

            std::vector<Point> batch_r0;
            for (int c = 0; c < 64; ++c) {
                g.setObserved(0, c, CellState::Free);
                batch_r0.push_back({0, c});
            }
            oracle.onBatchCellsOpened(batch_r0, g);

            GridMap g128(128, 128);
            for (int r = 0; r < 128; ++r) {
                for (int c = 0; c < 128; ++c) {
                    if (r < 64 && c < 64) {
                        g128.setObserved(r, c, g.isPassable(r, c) ? CellState::Free : CellState::Obstacle);
                    } else {
                        g128.setObserved(r, c, CellState::Free);
                    }
                }
            }
            oracle.liftRoot(g128);

            Point p_old{0, 0};
            Point p_new{0, 100};
            dist_t d_prop = oracle.queryDistance(p_old, p_new, g128);
            dist_t d_bfs = BfsOracle::queryDistance(g128, p_old, p_new);
            if (d_prop != d_bfs || d_prop == INF) {
                std::cerr << "FAIL: Suite 16.4 liftRoot after batch mismatch: prop=" << d_prop << " vs bfs=" << d_bfs << "\n";
                return 1;
            }
            std::cout << "  Sub-test 16.4 (liftRoot after batch openings): PASSED\n";
        }
        std::cout << "  STATUS: PASSED\n\n";
    }

    std::cout << "=========================================================\n";
    std::cout << "  ALL EXACTNESS VERIFICATION SUITES PASSED PERFECTLY!\n";
    std::cout << "=========================================================\n";
    return 0;
}
