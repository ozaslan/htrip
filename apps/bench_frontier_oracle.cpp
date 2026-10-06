#include "htrip/grid_map.hpp"
#include "htrip/frontier.hpp"
#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <iomanip>
#include <random>

using namespace htrip;

struct BenchmarkRecord {
    std::string config_name;
    int k_frontiers = 0;
    int b = 16;
    int g = 2;
    int depth = 0;
    bool hybrid_fallback = false;
    double avg_doorway_us = 0.0;
    double avg_blockage_us = 0.0;
    double avg_query_us = 0.0;
    double avg_total_us = 0.0;
    double avg_bfs_us = 0.0;
    double doorway_speedup = 0.0;
    double e2e_speedup = 0.0;
    size_t memory_kb = 0;
    size_t effective_kb = 0;
    int discrepancies = 0;
};

int main(int argc, char** argv) {
    std::string map_path = "procedural";
    int tile_size = 16;
    int group_size = 2;
    int steps = 20;
    int hybrid_k = 0;
    bool full_query = false;
    std::string csv_path = "frontier_results.csv";
    std::vector<int> k_values = {10, 25, 50};

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--map" && i + 1 < argc) map_path = argv[++i];
        else if (arg == "--tile-size" && i + 1 < argc) tile_size = std::atoi(argv[++i]);
        else if (arg == "--group-size" && i + 1 < argc) group_size = std::atoi(argv[++i]);
        else if (arg == "--steps" && i + 1 < argc) steps = std::atoi(argv[++i]);
        else if (arg == "--hybrid-k" && i + 1 < argc) hybrid_k = std::atoi(argv[++i]);
        else if (arg == "--full-query") full_query = true;
        else if (arg == "--csv" && i + 1 < argc) csv_path = argv[++i];
    }

    std::cout << "=========================================================\n";
    std::cout << "  ICRA 2027: H-TRIP (Hierarchical Tropical) Benchmark\n";
    std::cout << "=========================================================\n";
    std::cout << "Map: " << map_path << ", Tile size b=" << tile_size
              << ", Group size g=" << group_size << ", Hybrid k0=" << hybrid_k << "\n\n";

    GridMap base_grid;
    if (map_path == "procedural") {
        base_grid.generateProceduralRooms(128, 128, 16, 777);
    } else {
        if (!base_grid.loadMovingAI(map_path)) {
            std::cerr << "Failed to load map: " << map_path << "\n";
            return 1;
        }
    }

    base_grid.initExplorationSnapshot({1, 1}, 0.50, 0.05, 42);

    HierarchyConfig cfg;
    cfg.tile_size_b = tile_size;
    cfg.group_w = group_size;
    cfg.group_h = group_size;
    // Hybrid dispatch is handled explicitly below (including skipping index
    // maintenance), so the oracle's own threshold must stay disabled here.

    HTripOracle base_oracle(cfg);
    base_oracle.build(base_grid);

    FrontierManager fm(2);
    auto all_clusters = fm.detectAndCluster(base_grid);
    std::cout << "Detected " << all_clusters.size() << " initial frontier clusters.\n\n";

    std::vector<BenchmarkRecord> records;

    for (int k : k_values) {
        if (static_cast<size_t>(k) > all_clusters.size()) {
            std::cout << "Skipping k=" << k << " (only " << all_clusters.size() << " frontiers available)\n";
            continue;
        }

        std::vector<Point> viewpoints;
        for (int i = 0; i < k; ++i) {
            viewpoints.push_back(all_clusters[static_cast<size_t>(i)].representative);
        }

        // Snapshot isolation: fresh copy of grid and oracle for every k
        GridMap grid = base_grid;
        HTripOracle oracle = base_oracle;

        BenchmarkRecord rec;
        rec.k_frontiers = k;
        rec.config_name = (group_size == 2 ? "H-TRIP (Quad)" : "H-TRIP (16-ary)");
        rec.b = tile_size;
        rec.g = group_size;
        rec.depth = oracle.depth();
        rec.hybrid_fallback = (hybrid_k > 0 && k <= hybrid_k);
        rec.memory_kb = oracle.totalMemoryBytes() / 1024;
        rec.effective_kb = oracle.effectiveMemoryBytes(grid) / 1024;

        double total_doorway_us = 0.0;
        double total_blockage_us = 0.0;
        double total_query_us = 0.0;
        double total_bfs_us = 0.0;
        int doorway_steps = 0;
        int blockage_steps = 0;

        std::mt19937_64 rng(12345);

        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState inc_state;
        std::vector<dist_t> prop_matrix;

        for (int s = 0; s < steps; ++s) {
            // Find an unexplored cell adjacent to explored cells to simulate doorway discovery
            std::vector<Point> doorway_cands;
            for (int r = 1; r < grid.rows - 1; ++r) {
                for (int c = 1; c < grid.cols - 1; ++c) {
                    Point p{r, c};
                    if (grid.getObserved(p) != CellState::Free && grid.ground_truth[grid.index(p)] == CellState::Free) {
                        doorway_cands.push_back(p);
                        if (doorway_cands.size() > 50) break;
                    }
                }
                if (doorway_cands.size() > 50) break;
            }

            if (!doorway_cands.empty()) {
                Point pt = doorway_cands[rng() % doorway_cands.size()];
                grid.openDoorway(pt);

                auto t0 = std::chrono::high_resolution_clock::now();
                if (!rec.hybrid_fallback) {
                    oracle.onCellOpened(pt, grid);
                }
                auto t1 = std::chrono::high_resolution_clock::now();
                total_doorway_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
                doorway_steps++;

                // Query and compare against BFS.
                // Hybrid fallback rows intentionally bypass the oracle (and its
                // index maintenance) so they measure a raw k x BFS pipeline.
                t0 = std::chrono::high_resolution_clock::now();
                if (rec.hybrid_fallback) {
                    prop_matrix = BfsOracle::queryAllPairs(grid, viewpoints);
                } else if (full_query) {
                    oracle.queryAllPairs(viewpoints, grid, prop_matrix, scratch);
                } else {
                    oracle.queryAllPairsIncremental(viewpoints, grid, prop_matrix, scratch, inc_state);
                }
                t1 = std::chrono::high_resolution_clock::now();
                total_query_us += std::chrono::duration<double, std::micro>(t1 - t0).count();

                t0 = std::chrono::high_resolution_clock::now();
                std::vector<dist_t> bfs_matrix = BfsOracle::queryAllPairs(grid, viewpoints);
                t1 = std::chrono::high_resolution_clock::now();
                total_bfs_us += std::chrono::duration<double, std::micro>(t1 - t0).count();

                for (size_t idx = 0; idx < prop_matrix.size(); ++idx) {
                    if (prop_matrix[idx] != bfs_matrix[idx]) {
                        rec.discrepancies++;
                        if (rec.discrepancies <= 8) {
                            const size_t ii = idx / viewpoints.size();
                            const size_t jj = idx % viewpoints.size();
                            const dist_t d_point = oracle.queryDistance(viewpoints[ii], viewpoints[jj], grid);
                            std::cerr << "  DIAG mismatch pair (" << viewpoints[ii].r << "," << viewpoints[ii].c
                                      << ")->(" << viewpoints[jj].r << "," << viewpoints[jj].c << "): prop="
                                      << prop_matrix[idx] << " bfs=" << bfs_matrix[idx]
                                      << " point=" << d_point
                                      << " lca=" << oracle.findLcaLevel(viewpoints[ii], viewpoints[jj]) << "\n";
                        }
                    }
                }
            }

            // Test dynamic blockage on an explored free corridor
            if (s % 4 == 0) {
                int r = static_cast<int>(rng() % static_cast<uint64_t>(grid.rows));
                int c = static_cast<int>(rng() % static_cast<uint64_t>(grid.cols));
                Point pt{r, c};
                if (grid.getObserved(pt) == CellState::Free) {
                    grid.blockCorridor(pt);

                    if (!rec.hybrid_fallback) {
                        auto t0 = std::chrono::high_resolution_clock::now();
                        oracle.onCellBlocked(pt, grid);
                        auto t1 = std::chrono::high_resolution_clock::now();
                        total_blockage_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
                        blockage_steps++;

                        // Re-open it to maintain test environment stability
                        grid.openDoorway(pt);
                        oracle.onCellOpened(pt, grid);
                    } else {
                        grid.openDoorway(pt);
                    }
                }
            }
        }

        rec.avg_doorway_us = doorway_steps > 0 ? total_doorway_us / doorway_steps : 0.0;
        rec.avg_blockage_us = blockage_steps > 0 ? total_blockage_us / blockage_steps : 0.0;
        rec.avg_query_us = doorway_steps > 0 ? total_query_us / doorway_steps : 0.0;
        rec.avg_total_us = rec.avg_doorway_us + rec.avg_query_us;
        rec.avg_bfs_us = doorway_steps > 0 ? total_bfs_us / doorway_steps : 0.0;

        rec.doorway_speedup = rec.avg_doorway_us > 0.0 ? (rec.avg_bfs_us / rec.avg_doorway_us) : 0.0;
        rec.e2e_speedup = rec.avg_total_us > 0.0 ? (rec.avg_bfs_us / rec.avg_total_us) : 0.0;

        records.push_back(rec);
    }

    // Print summary table
    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::left << std::setw(12) << "Config"
              << std::setw(5) << "k"
              << std::setw(7) << "Depth"
              << std::setw(8) << "Mode"
              << std::setw(14) << "Doorway (us)"
              << std::setw(14) << "Query (us)"
              << std::setw(14) << "Total (us)"
              << std::setw(14) << "k x BFS (us)"
              << std::setw(14) << "Update Spdup"
              << std::setw(12) << "E2E Spdup"
              << std::setw(8) << "Errors" << "\n";
    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------\n";

    for (const auto& r : records) {
        std::ostringstream d_ss, e_ss;
        d_ss << std::fixed << std::setprecision(1) << r.doorway_speedup << "x";
        e_ss << std::fixed << std::setprecision(1) << r.e2e_speedup << "x";
        std::cout << std::left << std::setw(12) << r.config_name
                  << std::setw(5) << r.k_frontiers
                  << std::setw(7) << r.depth
                  << std::setw(8) << (r.hybrid_fallback ? "BFS" : "H-TRIP")
                  << std::setw(14) << std::fixed << std::setprecision(1) << r.avg_doorway_us
                  << std::setw(14) << std::fixed << std::setprecision(1) << r.avg_query_us
                  << std::setw(14) << std::fixed << std::setprecision(1) << r.avg_total_us
                  << std::setw(14) << std::fixed << std::setprecision(1) << r.avg_bfs_us
                  << std::setw(14) << d_ss.str()
                  << std::setw(12) << e_ss.str()
                  << std::setw(8) << r.discrepancies << "\n";
    }
    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------\n";

    // Write CSV
    std::ofstream csv(csv_path);
    if (csv.is_open()) {
        csv << "Config,k,b,g,Depth,Mode,Doorway_us,Query_us,Total_Latency_us,BFS_us,Update_Speedup,End_to_End_Speedup,Memory_KB,Effective_KB,Errors\n";
        for (const auto& r : records) {
            csv << r.config_name << "," << r.k_frontiers << "," << r.b << "," << r.g << ","
                << r.depth << "," << (r.hybrid_fallback ? "BFS" : "H-TRIP") << ","
                << r.avg_doorway_us << "," << r.avg_query_us << ","
                << r.avg_total_us << "," << r.avg_bfs_us << "," << r.doorway_speedup << ","
                << r.e2e_speedup << "," << r.memory_kb << "," << r.effective_kb << ","
                << r.discrepancies << "\n";
        }
        std::cout << "\nResults written to " << csv_path << "\n";
    }

    return 0;
}
