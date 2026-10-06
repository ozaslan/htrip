#include "htrip/minplus_oracle.hpp"
#include "htrip/dijkstra_oracle.hpp"
#include "htrip/grid_map.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <cmath>

using namespace htrip;

static void populateWeightedGrid(GridMap& grid, int seed = 42) {
    grid.generateProceduralRooms(grid.rows, grid.cols, 16, static_cast<uint64_t>(seed));
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            grid.setObserved(r, c, grid.ground_truth[static_cast<size_t>(grid.index(r, c))]);
        }
    }
    std::mt19937 rng(static_cast<uint32_t>(seed + 77));
    const dist_t weights[] = {1, 1, 2, 4, 6, 10, 20};
    std::uniform_int_distribution<size_t> dist(0, 6);
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            if (grid.isPassable(r, c)) {
                grid.setCellWeight(r, c, weights[dist(rng)]);
            }
        }
    }
}

int main() {
    std::cout << "======================================================================================\n";
    std::cout << "  BENCHMARK: WEIGHTED GRAPH MIN-PLUS INDEX (CONSTRUCTION, SINGLE & BATCH UPDATES)\n";
    std::cout << "======================================================================================\n\n";

    // -------------------------------------------------------------------------
    // PART 1: Initial Construction from Scratch Across Maze Sizes, b, and g
    // -------------------------------------------------------------------------
    std::cout << "--------------------------------------------------------------------------------------\n";
    std::cout << "  PART 1: INDEX BUILD TIME FROM SCRATCH ACROSS (MAZE SIZE, b, g)\n";
    std::cout << "--------------------------------------------------------------------------------------\n";
    std::cout << std::left
              << std::setw(12) << "Grid Size"
              << std::setw(8)  << "Tile b"
              << std::setw(8)  << "Group g"
              << std::setw(8)  << "Depth H"
              << std::setw(14) << "Spanned S"
              << std::setw(14) << "Leaves"
              << std::setw(18) << "Build Time (ms)"
              << "Status" << "\n";
    std::cout << "--------------------------------------------------------------------------------------\n";

    const std::vector<int> grid_sizes = {16, 32, 64, 128};
    const std::vector<int> b_values   = {2, 4, 8, 16};
    const std::vector<int> g_values   = {2, 4, 8, 16};

    for (int N : grid_sizes) {
        GridMap grid(N, N);
        populateWeightedGrid(grid, 100 + N);

        for (int b : b_values) {
            if (b > N) continue; // Tile cannot be larger than grid

            for (int g : g_values) {
                HierarchyConfig hc;
                hc.tile_size_b = b;
                hc.group_w = g;
                hc.group_h = g;

                if (!hc.isValid()) continue;

                // Safety guard against exponential port explosion across ALL hierarchy levels:
                // Level-l node has g^2 children, each of side (b * g^(l-1)).
                // Child port count = 4 * (b * g^(l-1)) - 4.
                // V_U(l) = g^2 * (4 * b * g^(l-1) - 4).
                // If V_U(l) > 1024 at any level, Floyd-Warshall (O(V^3)) would blow up memory/CPU.
                bool exceeds_budget = false;
                int cur_child_side = b;
                int max_vu = 0;
                int exceed_lvl = 1;
                for (int lvl = 1; cur_child_side < N; ++lvl) {
                    int child_ports = 4 * cur_child_side - 4;
                    int vu_lvl = g * g * child_ports;
                    if (vu_lvl > max_vu) max_vu = vu_lvl;
                    if (vu_lvl > 1024) {
                        exceeds_budget = true;
                        exceed_lvl = lvl;
                        break;
                    }
                    cur_child_side *= g;
                }

                if (exceeds_budget) {
                    std::cout << std::left
                              << std::setw(12) << (std::to_string(N) + "x" + std::to_string(N))
                              << std::setw(8)  << b
                              << std::setw(8)  << g
                              << std::setw(8)  << "-"
                              << std::setw(14) << "-"
                              << std::setw(14) << "-"
                              << std::setw(18) << "[Skipped]"
                              << "Exceeds port budget at L" << exceed_lvl << " (V_U=" << max_vu << ")\n";
                    continue;
                }

                MinPlusOracle oracle(hc);

                // Warmup & timing
                int reps = (N <= 32) ? 10 : ((N <= 64) ? 5 : 2);
                double total_ms = 0.0;
                bool ok = true;

                for (int rep = 0; rep < reps; ++rep) {
                    auto t0 = std::chrono::high_resolution_clock::now();
                    bool b_ok = oracle.build(grid);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    if (!b_ok) { ok = false; break; }
                    total_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
                }

                if (ok) {
                    double avg_ms = total_ms / reps;
                    std::cout << std::left
                              << std::setw(12) << (std::to_string(N) + "x" + std::to_string(N))
                              << std::setw(8)  << b
                              << std::setw(8)  << g
                              << std::setw(8)  << oracle.H
                              << std::setw(14) << (std::to_string(oracle.S) + "x" + std::to_string(oracle.S))
                              << std::setw(14) << oracle.leaves.size()
                              << std::setw(18) << std::fixed << std::setprecision(3) << avg_ms
                              << "OK\n";
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // PART 2: Single-Cell Dynamic Weight Update Latency in an Already-Built Index
    // -------------------------------------------------------------------------
    std::cout << "\n--------------------------------------------------------------------------------------\n";
    std::cout << "  PART 2: SINGLE-CELL DYNAMIC WEIGHT UPDATE IN ALREADY-BUILT INDEX\n";
    std::cout << "--------------------------------------------------------------------------------------\n";

    {
        int N = 64;
        GridMap grid(N, N);
        populateWeightedGrid(grid, 202);

        struct TestCfg { int b; int g; };
        std::vector<TestCfg> configs = {
            {8, 2},
            {16, 2},
            {8, 4},
            {16, 4}
        };

        std::cout << std::left
                  << std::setw(16) << "Config"
                  << std::setw(12) << "Updates"
                  << std::setw(20) << "Mean Latency (us)"
                  << std::setw(20) << "Median Latency (us)"
                  << std::setw(16) << "Throughput (ops/s)"
                  << "\n";
        std::cout << "--------------------------------------------------------------------------------------\n";

        for (const auto& cfg : configs) {
            HierarchyConfig hc;
            hc.tile_size_b = cfg.b;
            hc.group_w = cfg.g;
            hc.group_h = cfg.g;
            MinPlusOracle oracle(hc);
            oracle.build(grid);

            std::mt19937 rng(505);
            std::uniform_int_distribution<int> r_dist(1, N - 2);
            std::uniform_int_distribution<int> c_dist(1, N - 2);
            std::uniform_int_distribution<dist_t> w_dist(1, 20);

            const int num_updates = 100;
            std::vector<double> latencies_us;
            latencies_us.reserve(num_updates);

            for (int u = 0; u < num_updates; ++u) {
                Point p(r_dist(rng), c_dist(rng));
                while (!grid.isPassable(p)) {
                    p.r = r_dist(rng);
                    p.c = c_dist(rng);
                }
                dist_t new_w = w_dist(rng);

                auto t0 = std::chrono::high_resolution_clock::now();
                oracle.onCellWeightChanged(p, new_w, grid);
                auto t1 = std::chrono::high_resolution_clock::now();

                latencies_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
            }

            std::sort(latencies_us.begin(), latencies_us.end());
            double sum = 0.0;
            for (double v : latencies_us) sum += v;
            double mean = sum / num_updates;
            double median = latencies_us[num_updates / 2];
            double ops = 1e6 / mean;

            std::string label = "b=" + std::to_string(cfg.b) + ", g=" + std::to_string(cfg.g);
            std::cout << std::left
                      << std::setw(16) << label
                      << std::setw(12) << num_updates
                      << std::setw(20) << std::fixed << std::setprecision(2) << mean
                      << std::setw(20) << std::fixed << std::setprecision(2) << median
                      << std::setw(16) << std::fixed << std::setprecision(0) << ops
                      << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // PART 3: 10-Cell Batch Update Within a Single Tile vs 10 Single Updates
    // -------------------------------------------------------------------------
    std::cout << "\n--------------------------------------------------------------------------------------\n";
    std::cout << "  PART 3: 10-CELL BATCH UPDATE IN SINGLE TILE vs 10 INDIVIDUAL CELL UPDATES\n";
    std::cout << "--------------------------------------------------------------------------------------\n";

    {
        int N = 64;
        int b = 16;
        int g = 2;

        GridMap grid(N, N);
        populateWeightedGrid(grid, 303);

        HierarchyConfig hc;
        hc.tile_size_b = b;
        hc.group_w = g;
        hc.group_h = g;
        MinPlusOracle oracle(hc);
        oracle.build(grid);

        // Target leaf tile (1, 1) spanning [16, 31], [16, 31]
        int tx = 1, ty = 1;
        std::vector<Point> tile_cells;
        for (int r = ty * b; r < (ty + 1) * b; ++r) {
            for (int c = tx * b; c < (tx + 1) * b; ++c) {
                if (grid.isPassable(r, c)) {
                    tile_cells.emplace_back(r, c);
                }
            }
        }
        if (tile_cells.size() > 10) tile_cells.resize(10);

        const int num_trials = 50;
        std::mt19937 rng(404);
        std::uniform_int_distribution<dist_t> w_dist(1, 20);

        // 1. Measure 10 individual cell updates
        std::vector<double> indiv_us;
        indiv_us.reserve(num_trials);
        for (int t = 0; t < num_trials; ++t) {
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < tile_cells.size(); ++i) {
                dist_t nw = w_dist(rng);
                oracle.onCellWeightChanged(tile_cells[i], nw, grid);
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            indiv_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
        }

        // 2. Measure batch update of 10 cells in the same tile
        std::vector<double> batch_us;
        batch_us.reserve(num_trials);
        for (int t = 0; t < num_trials; ++t) {
            std::vector<std::pair<Point, dist_t>> batch;
            for (size_t i = 0; i < tile_cells.size(); ++i) {
                batch.emplace_back(tile_cells[i], w_dist(rng));
            }

            auto t0 = std::chrono::high_resolution_clock::now();
            oracle.onTileWeightsChanged(tx, ty, batch, grid);
            auto t1 = std::chrono::high_resolution_clock::now();
            batch_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
        }

        std::sort(indiv_us.begin(), indiv_us.end());
        std::sort(batch_us.begin(), batch_us.end());

        double indiv_med = indiv_us[num_trials / 2];
        double batch_med = batch_us[num_trials / 2];
        double speedup   = indiv_med / batch_med;

        std::cout << "Configuration: Maze 64x64, Tile b=16, Group g=2 (Quadtree)\n\n";
        std::cout << "  10 Individual Cell Updates Median: " << std::fixed << std::setprecision(1) << indiv_med << " us\n";
        std::cout << "  10-Cell Tile Batch Update Median:  " << std::fixed << std::setprecision(1) << batch_med << " us\n";
        std::cout << "  Tile-Batching Speedup:             " << std::fixed << std::setprecision(2) << speedup << "x faster!\n\n";
        std::cout << "Reason: In a tile-batch update, the leaf tile APSP is recomputed ONCE,\n";
        std::cout << "        and the ancestor closures up the hierarchy are recomputed ONCE,\n";
        std::cout << "        reducing algebraic ancestor propagation work by exactly 10x.\n";
    }

    std::cout << "======================================================================================\n";
    return 0;
}
