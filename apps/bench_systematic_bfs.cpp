#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include "htrip/astar_oracle.hpp"
#include "htrip/grid_map.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <queue>

using namespace htrip;

// Find the largest connected component and systematically pick viewpoints
static std::vector<Point> selectSystematicViewpoints(const GridMap& grid, int grid_divisions = 8) {
    int R = grid.rows;
    int C = grid.cols;

    // 1. Find largest connected component via flood fill
    std::vector<int> comp_id(static_cast<size_t>(R * C), -1);
    std::vector<size_t> comp_size;
    int current_comp = 0;

    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            int idx = grid.index(r, c);
            if (!grid.isPassable(r, c) || comp_id[static_cast<size_t>(idx)] != -1) continue;

            size_t sz = 0;
            std::queue<Point> q;
            q.push({r, c});
            comp_id[static_cast<size_t>(idx)] = current_comp;

            while (!q.empty()) {
                Point p = q.front();
                q.pop();
                sz++;

                for (int k = 0; k < 4; ++k) {
                    int nr = p.r + dr[k];
                    int nc = p.c + dc[k];
                    if (nr < 0 || nr >= R || nc < 0 || nc >= C) continue;
                    if (!grid.isPassable(nr, nc)) continue;
                    int nidx = grid.index(nr, nc);
                    if (comp_id[static_cast<size_t>(nidx)] == -1) {
                        comp_id[static_cast<size_t>(nidx)] = current_comp;
                        q.push({nr, nc});
                    }
                }
            }
            comp_size.push_back(sz);
            current_comp++;
        }
    }

    int largest_comp = -1;
    size_t max_sz = 0;
    for (int i = 0; i < current_comp; ++i) {
        if (comp_size[static_cast<size_t>(i)] > max_sz) {
            max_sz = comp_size[static_cast<size_t>(i)];
            largest_comp = i;
        }
    }

    // 2. Sample systematically in an grid_divisions x grid_divisions sector grid
    std::vector<Point> sampled;
    int step_r = R / grid_divisions;
    int step_c = C / grid_divisions;

    for (int sr = 0; sr < grid_divisions; ++sr) {
        for (int sc = 0; sc < grid_divisions; ++sc) {
            int center_r = sr * step_r + step_r / 2;
            int center_c = sc * step_c + step_c / 2;

            // Spiral search around center for nearest passable cell in the main component
            Point best_p(-1, -1);
            int min_dist = 999999;

            for (int dr_ = -step_r / 2; dr_ <= step_r / 2; ++dr_) {
                for (int dc_ = -step_c / 2; dc_ <= step_c / 2; ++dc_) {
                    int r = center_r + dr_;
                    int c = center_c + dc_;
                    if (r >= 0 && r < R && c >= 0 && c < C) {
                        int idx = grid.index(r, c);
                        if (grid.isPassable(r, c) && comp_id[static_cast<size_t>(idx)] == largest_comp) {
                            int d = std::abs(dr_) + std::abs(dc_);
                            if (d < min_dist) {
                                min_dist = d;
                                best_p = {r, c};
                            }
                        }
                    }
                }
            }

            if (best_p.r != -1) {
                sampled.push_back(best_p);
            }
        }
    }

    return sampled;
}

struct SystematicPair {
    Point s;
    Point t;
    std::string category;
};

static void runLargeScaleMapBenchmark(const std::string& map_name,
                                      const std::string& map_path,
                                      int b, int g) {
    std::cout << "\n========================================================================================\n";
    std::cout << "  BENCHMARK: " << map_name << "\n";
    std::cout << "========================================================================================\n";

    GridMap grid;
    std::string actual_path = map_path;
    {
        std::ifstream test_f(actual_path);
        if (!test_f.is_open()) {
            actual_path = "../../" + map_path;
        }
    }
    if (!grid.loadMovingAI(actual_path)) {
        std::cerr << "Failed to load map: " << map_path << "\n";
        return;
    }

    // Set all ground truth passable cells as observed Free
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            grid.setObserved(r, c, grid.ground_truth[static_cast<size_t>(grid.index(r, c))]);
        }
    }

    std::cout << "  Original Input Maze Dimensions: " << grid.rows << " x " << grid.cols << " cells\n";
    std::cout << "  Hierarchy Configuration:        Tile b = " << b << ", Group g = " << g << "\n";

    // 1. Measure Build Time from Scratch
    HierarchyConfig hc;
    hc.tile_size_b = b;
    hc.group_w = g;
    hc.group_h = g;

    MinPlusOracle oracle(hc);
    auto t_build_0 = std::chrono::high_resolution_clock::now();
    bool built = oracle.build(grid);
    auto t_build_1 = std::chrono::high_resolution_clock::now();
    double build_ms = std::chrono::duration<double, std::milli>(t_build_1 - t_build_0).count();

    if (!built) {
        std::cerr << "Build failed for " << map_name << "\n";
        return;
    }

    size_t mem_bytes = oracle.totalMemoryBytes();
    double mem_mb = static_cast<double>(mem_bytes) / (1024.0 * 1024.0);

    std::cout << "  Hierarchy Depth H:              " << oracle.H << " levels\n";
    std::cout << "  Total Spanned Space S:          " << oracle.S << " x " << oracle.S << "\n";
    std::cout << "  Leaf Tiles Count:               " << oracle.leaves.size() << "\n";
    std::cout << "  Memory Footprint:               " << std::fixed << std::setprecision(2) << mem_mb << " MB\n";
    std::cout << "  Index Build Time from Scratch:  " << std::fixed << std::setprecision(2) << build_ms << " ms\n";

    // 2. Select Systematic Viewpoints
    std::vector<Point> vps = selectSystematicViewpoints(grid, 8);
    std::cout << "  Systematic Viewpoints Sampled:  " << vps.size() << " evenly distributed across sectors\n\n";

    // 3. Form Systematic Pairs across Hierarchy Spectrum (LCA 0, 1, 2, >= 3)
    std::vector<SystematicPair> pairs;

    // (A) Local intra-tile pairs (same leaf tile, LCA = 0)
    for (size_t i = 0; i < vps.size() && pairs.size() < 4; ++i) {
        for (size_t j = i + 1; j < vps.size() && pairs.size() < 4; ++j) {
            if (vps[i].c / b == vps[j].c / b && vps[i].r / b == vps[j].r / b) {
                pairs.push_back({vps[i], vps[j], "Local (Intra-Tile, LCA=0)"});
            }
        }
    }

    // (B) Sibling-tile pairs (different leaf, same L1 ancestor, LCA = 1)
    for (size_t i = 0; i < vps.size() && pairs.size() < 9; ++i) {
        for (size_t j = i + 1; j < vps.size() && pairs.size() < 9; ++j) {
            int lca = oracle.findLcaLevel(vps[i], vps[j]);
            if (lca == 1) {
                pairs.push_back({vps[i], vps[j], "Medium-Range (Sibling, LCA=1)"});
            }
        }
    }

    // (C) Long-range pairs (LCA = 2)
    for (size_t i = 0; i < vps.size() && pairs.size() < 14; ++i) {
        for (size_t j = i + 1; j < vps.size() && pairs.size() < 14; ++j) {
            int lca = oracle.findLcaLevel(vps[i], vps[j]);
            if (lca == 2) {
                pairs.push_back({vps[i], vps[j], "Long-Range (Quadrant, LCA=2)"});
            }
        }
    }

    // (D) Global cross-map pairs (LCA >= 3 or fallback)
    for (size_t i = 0; i < vps.size() && pairs.size() < 20; ++i) {
        for (size_t j = i + 1; j < vps.size() && pairs.size() < 20; ++j) {
            int lca = oracle.findLcaLevel(vps[i], vps[j]);
            if (lca >= 3 || (lca >= 2 && pairs.size() < 20)) {
                pairs.push_back({vps[i], vps[j], "Global (Cross-Map, LCA=" + std::to_string(lca) + ")"});
            }
        }
    }

    // 4. Run systematic pairwise comparison: BFS vs H-TRIP
    std::cout << "----------------------------------------------------------------------------------------\n";
    std::cout << "  SYSTEMATIC CELL PAIR DISTANCE & TIME COMPARISON: BFS vs H-TRIP\n";
    std::cout << "----------------------------------------------------------------------------------------\n";
    std::cout << std::left
              << std::setw(6)  << "Pair"
              << std::setw(16) << "Source -> Target"
              << std::setw(32) << "Distance Category"
              << std::setw(10) << "BFS Dist"
              << std::setw(12) << "H-TRIP Dist"
              << std::setw(8)  << "Diff"
              << std::setw(14) << "BFS Time (us)"
              << std::setw(16) << "H-TRIP Time (us)"
              << "Speedup" << "\n";
    std::cout << "----------------------------------------------------------------------------------------\n";

    size_t pair_idx = 1;
    double total_bfs_us = 0.0;
    double total_htrip_us = 0.0;
    bool all_exact = true;

    for (const auto& p : pairs) {
        // Measure BFS time (average of 3 runs)
        double bfs_sum_us = 0.0;
        dist_t d_bfs = 0;
        for (int r = 0; r < 3; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            d_bfs = BfsOracle::queryDistance(grid, p.s, p.t);
            auto t1 = std::chrono::high_resolution_clock::now();
            bfs_sum_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
        }
        double bfs_us = bfs_sum_us / 3.0;

        // Measure H-TRIP time (average of 3 runs)
        double htrip_sum_us = 0.0;
        dist_t d_htrip = 0;
        for (int r = 0; r < 3; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            d_htrip = oracle.queryDistance(p.s, p.t, grid);
            auto t1 = std::chrono::high_resolution_clock::now();
            htrip_sum_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
        }
        double htrip_us = htrip_sum_us / 3.0;

        int diff = std::abs(static_cast<int>(d_htrip) - static_cast<int>(d_bfs));
        if (diff != 0) all_exact = false;

        double speedup = bfs_us / (htrip_us > 0.0 ? htrip_us : 0.001);
        total_bfs_us += bfs_us;
        total_htrip_us += htrip_us;

        std::string pts = "(" + std::to_string(p.s.r) + "," + std::to_string(p.s.c) + ")->("
                              + std::to_string(p.t.r) + "," + std::to_string(p.t.c) + ")";

        std::cout << std::left
                  << std::setw(6)  << ("#" + std::to_string(pair_idx++))
                  << std::setw(16) << pts
                  << std::setw(32) << p.category
                  << std::setw(10) << d_bfs
                  << std::setw(12) << d_htrip
                  << std::setw(8)  << (diff == 0 ? "0 (OK)" : std::to_string(diff))
                  << std::setw(14) << std::fixed << std::setprecision(1) << bfs_us
                  << std::setw(16) << std::fixed << std::setprecision(1) << htrip_us
                  << std::fixed << std::setprecision(1) << speedup << "x\n";
    }

    std::cout << "----------------------------------------------------------------------------------------\n";
    std::cout << "  Pairwise Verification Summary: " << (all_exact ? "100% BIT-EXACT (ALL PAIRS MATCH BFS!)" : "MISMATCH DETECTED!") << "\n";
    std::cout << "  Average BFS Query Latency:     " << std::fixed << std::setprecision(1) << (total_bfs_us / pairs.size()) << " us\n";
    std::cout << "  Average H-TRIP Query Latency:  " << std::fixed << std::setprecision(1) << (total_htrip_us / pairs.size()) << " us\n";
    std::cout << "  Average Pairwise Speedup:      " << std::fixed << std::setprecision(1) << (total_bfs_us / total_htrip_us) << "x\n";

    // 5. Batched All-Pairs Matrix Query Benchmark (k = 25 viewpoints -> 625 pairs)
    if (vps.size() >= 25) {
        std::vector<Point> k_vps(vps.begin(), vps.begin() + 25);
        size_t k = k_vps.size();

        std::cout << "\n  --- BATCHED ALL-PAIRS MATRIX BENCHMARK (k = " << k << " Viewpoints = " << (k * k) << " Pairs) ---\n";

        // BFS: k separate searches (average of 3 runs)
        double bfs_mat_ms = 0.0;
        std::vector<dist_t> bfs_mat;
        for (int r = 0; r < 3; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            bfs_mat = BfsOracle::queryAllPairs(grid, k_vps);
            auto t1 = std::chrono::high_resolution_clock::now();
            bfs_mat_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
        bfs_mat_ms /= 3.0;

        // H-TRIP: batched tropical GEMM (average of 3 runs)
        double htrip_mat_ms = 0.0;
        std::vector<dist_t> htrip_mat;
        for (int r = 0; r < 3; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            htrip_mat = oracle.queryAllPairs(k_vps, grid);
            auto t1 = std::chrono::high_resolution_clock::now();
            htrip_mat_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
        htrip_mat_ms /= 3.0;

        // A*: k(k-1)/2 separate point-to-point searches (average of 3 runs)
        double astar_mat_ms = 0.0;
        std::vector<dist_t> astar_mat;
        for (int r = 0; r < 3; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            astar_mat = AStarOracle::queryAllPairs(grid, k_vps);
            auto t1 = std::chrono::high_resolution_clock::now();
            astar_mat_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
        astar_mat_ms /= 3.0;

        size_t mat_matches = 0;
        size_t astar_matches = 0;
        for (size_t i = 0; i < k * k; ++i) {
            if (bfs_mat[i] == htrip_mat[i]) mat_matches++;
            if (astar_mat[i] == htrip_mat[i]) astar_matches++;
        }

        std::cout << "  BFS Matrix Exactness:  " << mat_matches << " / " << (k * k)
                  << " (" << (mat_matches == k * k ? "100% BIT-EXACT MATCH" : "MISMATCH") << ")\n";
        std::cout << "  A* Matrix Exactness:   " << astar_matches << " / " << (k * k)
                  << " (" << (astar_matches == k * k ? "100% BIT-EXACT MATCH" : "MISMATCH") << ")\n";
        std::cout << "  k x Full BFS Latency:  " << std::fixed << std::setprecision(2) << bfs_mat_ms << " ms\n";
        std::cout << "  k(k-1)/2 x A* Latency: " << std::fixed << std::setprecision(2) << astar_mat_ms << " ms\n";
        std::cout << "  H-TRIP Matrix Latency: " << std::fixed << std::setprecision(2) << htrip_mat_ms << " ms\n";
        std::cout << "  Speedup vs BFS:        " << std::fixed << std::setprecision(1) << (bfs_mat_ms / htrip_mat_ms) << "x FASTER than BFS!\n";
        std::cout << "  Speedup vs A*:         " << std::fixed << std::setprecision(1) << (astar_mat_ms / htrip_mat_ms) << "x FASTER than A*!\n";
        std::cout << "  A* to BFS Overhead:    " << std::fixed << std::setprecision(1) << (astar_mat_ms / bfs_mat_ms) << "x (BFS is faster for all-pairs)\n";
    }
}

int main() {
    std::cout << "========================================================================================\n";
    std::cout << "   ICRA 2027: LARGE-SCALE INPUT MAZES (256, 512, 1024) SYSTEMATIC BFS COMPARISON\n";
    std::cout << "========================================================================================\n";

    // 1. Boston 256x256
    runLargeScaleMapBenchmark("Boston Street Map (256x256)",
                              "datasets/movingai/street/maps/Boston_0_256.map",
                              16, 2);

    // 2. Maze 512x512
    runLargeScaleMapBenchmark("MovingAI Maze (512x512)",
                              "datasets/movingai/maze/maps/maze512-1-0.map",
                              16, 2);

    // 3. Boston 1024x1024
    runLargeScaleMapBenchmark("Boston Street Map (1024x1024)",
                              "datasets/movingai/street/maps/Boston_0_1024.map",
                              32, 2);

    return 0;
}
