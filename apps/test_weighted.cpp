#include "htrip/minplus_oracle.hpp"
#include "htrip/dijkstra_oracle.hpp"
#include "htrip/grid_map.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace htrip;

static const char* GREEN = "\033[32m";
static const char* RED   = "\033[31m";
static const char* RESET = "\033[0m";

static void createWeightedMaze(GridMap& grid, uint32_t seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<dist_t> weight_dist(1, 20);

    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            grid.setObserved(r, c, CellState::Free);
            grid.setCellWeight(r, c, weight_dist(rng));
        }
    }

    // Add obstacles (corridor walls)
    for (int r = 4; r < grid.rows; r += 8) {
        for (int c = 0; c < grid.cols; ++c) {
            if (c % 10 != 5) {
                grid.setObserved(r, c, CellState::Obstacle);
            }
        }
    }
}

static std::vector<Point> pickViewpoints(const GridMap& grid, size_t count, uint32_t seed = 123) {
    std::vector<Point> all_passable;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            if (grid.isPassable(r, c)) {
                all_passable.emplace_back(r, c);
            }
        }
    }
    std::mt19937 rng(seed);
    std::shuffle(all_passable.begin(), all_passable.end(), rng);
    if (all_passable.size() > count) {
        all_passable.resize(count);
    }
    return all_passable;
}

static bool verifyExactness(const MinPlusOracle& oracle,
                            const GridMap& grid,
                            const std::vector<Point>& viewpoints,
                            const std::string& label,
                            const std::vector<dist_t>* prev_matrix = nullptr,
                            bool expect_changes = false) {
    size_t k = viewpoints.size();
    std::vector<dist_t> htrip_mat = oracle.queryAllPairs(viewpoints, grid);
    std::vector<dist_t> dijkstra_mat = DijkstraOracle::queryAllPairs(grid, viewpoints);

    size_t matches = 0;
    size_t total = k * k;
    int max_diff = 0;
    size_t changed_from_prev = 0;

    for (size_t i = 0; i < total; ++i) {
        dist_t d_h = htrip_mat[i];
        dist_t d_g = dijkstra_mat[i];
        if (d_h == d_g) {
            matches++;
        } else {
            size_t src_idx = i / k;
            size_t dst_idx = i % k;
            Point ps = viewpoints[src_idx];
            Point pt = viewpoints[dst_idx];
            int diff = std::abs(static_cast<int>(d_h) - static_cast<int>(d_g));
            if (diff > max_diff) max_diff = diff;
            std::cout << "    Mismatch: (" << ps.r << "," << ps.c << ") -> ("
                      << pt.r << "," << pt.c << ") H-TRIP=" << d_h << ", Dijkstra=" << d_g
                      << ", diff=" << diff << "\n";
        }
        if (prev_matrix && (*prev_matrix)[i] != d_h) {
            changed_from_prev++;
        }
    }

    bool ok = (matches == total);
    if (expect_changes && changed_from_prev == 0) {
        std::cout << "    [CRITICAL ERROR] Update had no effect on any viewpoint pair!\n";
        ok = false;
    }

    std::cout << "  [" << label << "] Queries: " << total
              << ", Matches: " << matches << " / " << total
              << ", MaxDiff: " << max_diff;
    if (prev_matrix) {
        std::cout << ", Changed Pairs: " << changed_from_prev;
    }
    std::cout << " -> " << (ok ? GREEN : RED) << (ok ? "PASSED" : "FAILED") << RESET << "\n";
    return ok;
}

int main() {
    std::cout << "=================================================================\n";
    std::cout << "  ICRA 2027: WEIGHTED GRAPH EXACTNESS & DYNAMIC UNIT TESTS\n";
    std::cout << "=================================================================\n";

    bool all_ok = true;

    // -------------------------------------------------------------
    // Suite 1: Exactness on Weighted Graph Across Multiple Configurations
    // -------------------------------------------------------------
    std::cout << "\n[Suite 1: Build Exactness on Weighted Graph Across (b, g) Knobs]\n";
    {
        GridMap grid(64, 64);
        createWeightedMaze(grid, 42);
        auto vps = pickViewpoints(grid, 25, 101);

        struct Cfg { int b; int g; };
        std::vector<Cfg> configs = {
            {4, 2},
            {8, 2},
            {16, 2},
            {4, 4},
            {8, 4}
        };

        for (const auto& c : configs) {
            HierarchyConfig hc;
            hc.tile_size_b = c.b;
            hc.group_w = c.g;
            hc.group_h = c.g;
            MinPlusOracle oracle(hc);
            bool built = oracle.build(grid);
            std::string label = "b=" + std::to_string(c.b) + ", g=" + std::to_string(c.g);
            if (!built) {
                std::cout << "  [" << label << "] Build failed! -> " << RED << "FAILED" << RESET << "\n";
                all_ok = false;
                continue;
            }
            if (!verifyExactness(oracle, grid, vps, label)) {
                all_ok = false;
            }
        }
    }

    // -------------------------------------------------------------
    // Suite 2: Highly Discriminative Single-Cell Weight Updates (LCA >= 2)
    // -------------------------------------------------------------
    std::cout << "\n[Suite 2: Discriminative Single-Cell Weight Dynamic Updates (LCA >= 2)]\n";
    {
        GridMap grid(64, 64);
        createWeightedMaze(grid, 55);

        // Place a chokepoint corridor at row 12, col 25 (connects upper and lower sectors)
        Point p_choke(12, 25);
        grid.setObserved(p_choke, CellState::Free);
        grid.setCellWeight(p_choke, 30); // Start with heavy weight

        // Viewpoints positioned on opposite sides of p_choke (LCA >= 2)
        std::vector<Point> vps = {
            {2, 25},   // Quadrant 0 (top)
            {20, 25},  // Quadrant 2 (bottom)
            {2, 10},
            {20, 10},
            {5, 50},
            {50, 50}
        };

        HierarchyConfig hc;
        hc.tile_size_b = 8;
        hc.group_w = 2;
        hc.group_h = 2;
        MinPlusOracle oracle(hc);
        if (!oracle.build(grid)) {
            std::cerr << "Failed to build oracle for Suite 2\n";
            return 1;
        }

        // Verify baseline
        std::vector<dist_t> baseline_mat = oracle.queryAllPairs(vps, grid);
        if (!verifyExactness(oracle, grid, vps, "Step 0: Baseline")) all_ok = false;

        // Step 1: Weight decrease on p_choke (mud paved into superhighway: 30 -> 1)
        oracle.onCellWeightChanged(p_choke, 1, grid);
        std::vector<dist_t> dec_mat = oracle.queryAllPairs(vps, grid);
        if (!verifyExactness(oracle, grid, vps, "Step 1: Weight Decrease (30 -> 1)", &baseline_mat, true)) {
            all_ok = false;
        }

        // Verify strict distance reduction between opposite viewpoints
        dist_t d_base = baseline_mat[0 * vps.size() + 1];
        dist_t d_dec = dec_mat[0 * vps.size() + 1];
        if (d_dec >= d_base) {
            std::cout << "  [CRITICAL] d(P0, P1) did not strictly decrease (" << d_base << " -> " << d_dec << ")\n";
            all_ok = false;
        } else {
            std::cout << "  [Confirmed] d(P0, P1) strictly decreased: " << d_base << " -> " << d_dec << "\n";
        }

        // Step 2: Weight increase on p_choke (hazard added: 1 -> 60)
        oracle.onCellWeightChanged(p_choke, 60, grid);
        std::vector<dist_t> inc_mat = oracle.queryAllPairs(vps, grid);
        if (!verifyExactness(oracle, grid, vps, "Step 2: Weight Increase (1 -> 60)", &dec_mat, true)) {
            all_ok = false;
        }
        dist_t d_inc = inc_mat[0 * vps.size() + 1];
        if (d_inc <= d_dec) {
            std::cout << "  [CRITICAL] d(P0, P1) did not strictly increase (" << d_dec << " -> " << d_inc << ")\n";
            all_ok = false;
        } else {
            std::cout << "  [Confirmed] d(P0, P1) strictly increased: " << d_dec << " -> " << d_inc << "\n";
        }
    }

    // -------------------------------------------------------------
    // Suite 3: Batch Update of 10 Cells Inside a Single Tile
    // -------------------------------------------------------------
    std::cout << "\n[Suite 3: Batch Update of 10 Cells Within a Single Tile]\n";
    {
        GridMap grid(64, 64);
        createWeightedMaze(grid, 88);
        auto vps = pickViewpoints(grid, 20, 404);

        int b = 8;
        HierarchyConfig hc;
        hc.tile_size_b = b;
        hc.group_w = 2;
        hc.group_h = 2;
        MinPlusOracle oracle(hc);
        if (!oracle.build(grid)) {
            std::cerr << "Build failed for Suite 3\n";
            return 1;
        }

        std::vector<dist_t> base_mat = oracle.queryAllPairs(vps, grid);

        // Select tile (2, 2) which spans rows [16, 23], cols [16, 23]
        int tx = 2;
        int ty = 2;
        std::vector<Point> tile_passable;
        for (int r = ty * b; r < (ty + 1) * b; ++r) {
            for (int c = tx * b; c < (tx + 1) * b; ++c) {
                if (grid.isPassable(r, c)) {
                    tile_passable.emplace_back(r, c);
                }
            }
        }
        if (tile_passable.size() < 10) {
            std::cerr << "Tile does not contain at least 10 passable cells\n";
            return 1;
        }

        // Pick 10 cells in this tile and assign new weights
        std::vector<std::pair<Point, dist_t>> batch_weights;
        for (size_t i = 0; i < 10; ++i) {
            dist_t new_w = static_cast<dist_t>(1 + i * 2);
            batch_weights.emplace_back(tile_passable[i], new_w);
        }

        oracle.onTileWeightsChanged(tx, ty, batch_weights, grid);

        if (!verifyExactness(oracle, grid, vps, "10-Cell Batch Tile Update (tx=2, ty=2)", &base_mat, false)) {
            all_ok = false;
        }
    }

    // -------------------------------------------------------------
    // Suite 4: Seam Boundary Weight Perturbations
    // -------------------------------------------------------------
    std::cout << "\n[Suite 4: Seam Boundary Weight Perturbations (Inter-Tile)]\n";
    {
        GridMap grid(64, 64);
        createWeightedMaze(grid, 99);
        auto vps = pickViewpoints(grid, 20, 505);

        HierarchyConfig hc;
        hc.tile_size_b = 8;
        hc.group_w = 2;
        hc.group_h = 2;
        MinPlusOracle oracle(hc);
        if (!oracle.build(grid)) {
            std::cerr << "Build failed for Suite 4\n";
            return 1;
        }
        std::vector<dist_t> base_mat = oracle.queryAllPairs(vps, grid);

        // Perturb cells directly on the seam boundary at column 7 and 8
        for (int r = 10; r < 25; ++r) {
            if (grid.isPassable(r, 7)) {
                oracle.onCellWeightChanged({r, 7}, 12, grid);
            }
            if (grid.isPassable(r, 8)) {
                oracle.onCellWeightChanged({r, 8}, 12, grid);
            }
        }

        if (!verifyExactness(oracle, grid, vps, "Seam Crossing Weight Mod", &base_mat, true)) {
            all_ok = false;
        }
    }

    // -------------------------------------------------------------
    // Suite 5: queryAllPairsIncremental on Weighted Grid with Fallback
    // -------------------------------------------------------------
    std::cout << "\n[Suite 5: queryAllPairsIncremental on Weighted Grid]\n";
    {
        GridMap grid(64, 64);
        createWeightedMaze(grid, 77);
        auto vps = pickViewpoints(grid, 20, 707);

        HierarchyConfig hc;
        hc.tile_size_b = 8;
        hc.group_w = 2;
        hc.group_h = 2;
        MinPlusOracle oracle(hc);
        if (!oracle.build(grid)) {
            std::cerr << "Build failed for Suite 5\n";
            return 1;
        }

        MinPlusOracle::HTripIncrementalState inc_state;
        MinPlusOracle::HTripQueryScratch scratch;
        std::vector<dist_t> inc_mat;

        // 1. Initial incremental query
        oracle.queryAllPairsIncremental(vps, grid, inc_mat, scratch, inc_state);
        std::vector<dist_t> dijk_mat1 = DijkstraOracle::queryAllPairs(grid, vps);
        bool ok1 = (inc_mat == dijk_mat1);
        std::cout << "  [Step 1: Initial Incremental Query] -> "
                  << (ok1 ? GREEN : RED) << (ok1 ? "PASSED" : "FAILED") << RESET << "\n";
        if (!ok1) all_ok = false;

        // 2. Weight changes
        Point p_test(15, 15);
        if (grid.isPassable(p_test)) {
            oracle.onCellWeightChanged(p_test, 2, grid);
        }
        oracle.queryAllPairsIncremental(vps, grid, inc_mat, scratch, inc_state);
        std::vector<dist_t> dijk_mat2 = DijkstraOracle::queryAllPairs(grid, vps);
        bool ok2 = (inc_mat == dijk_mat2);
        std::cout << "  [Step 2: Incremental Query after Weight Change] -> "
                  << (ok2 ? GREEN : RED) << (ok2 ? "PASSED" : "FAILED") << RESET << "\n";
        if (!ok2) all_ok = false;

        // 3. Structural blockage (triggers epoch bump / structural fallback)
        Point p_block(20, 20);
        grid.setObserved(p_block, CellState::Obstacle);
        oracle.onCellBlocked(p_block, grid);

        oracle.queryAllPairsIncremental(vps, grid, inc_mat, scratch, inc_state);
        std::vector<dist_t> dijk_mat3 = DijkstraOracle::queryAllPairs(grid, vps);
        bool ok3 = (inc_mat == dijk_mat3);
        std::cout << "  [Step 3: Incremental Query after Cell Blockage (Epoch Bump)] -> "
                  << (ok3 ? GREEN : RED) << (ok3 ? "PASSED" : "FAILED") << RESET << "\n";
        if (!ok3) all_ok = false;
    }

    std::cout << "\n=================================================================\n";
    if (all_ok) {
        std::cout << GREEN << "  ALL WEIGHTED TESTS PASSED PERFECTLY (100% BIT-EXACT VS DIJKSTRA)!" << RESET << "\n";
    } else {
        std::cout << RED << "  SOME WEIGHTED TESTS FAILED!" << RESET << "\n";
    }
    std::cout << "=================================================================\n";

    return all_ok ? 0 : 1;
}
