#include "htrip/bfs_oracle.hpp"
#include "htrip/campaign_support.hpp"
#include "htrip/minplus_oracle.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace htrip;
using namespace htrip::campaign;

namespace {

using Clock = std::chrono::steady_clock;

int leafIndexOf(Point p, int b, int num_x) {
    const int tx = p.c / b;
    const int ty = p.r / b;
    return ty * num_x + tx;
}

std::vector<Point> chooseWithheld(const std::string& pattern, const GridMap& grid, const ConnectedComponent& cc,
                                  const std::vector<Point>& protected_pts, int b, std::uint64_t seed, int want) {
    const int num_x = (grid.cols + b - 1) / b;
    const int num_y = (grid.rows + b - 1) / b;
    const int n_leaves = std::max(1, num_x * num_y);
    std::vector<uint8_t> prot(static_cast<size_t>(grid.rows * grid.cols), 0);
    for (Point p : protected_pts) {
        if (grid.inBounds(p)) prot[static_cast<size_t>(grid.index(p))] = 1;
    }
    std::vector<std::vector<Point>> by_leaf(static_cast<size_t>(n_leaves));
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const int idx = grid.index(r, c);
            if (cc.cell_comp[static_cast<size_t>(idx)] != cc.id) continue;
            if (prot[static_cast<size_t>(idx)]) continue;
            const int li = leafIndexOf({r, c}, b, num_x);
            if (li >= 0 && li < n_leaves) by_leaf[static_cast<size_t>(li)].push_back({r, c});
        }
    }
    std::mt19937_64 rng(seed);
    for (auto& v : by_leaf) std::shuffle(v.begin(), v.end(), rng);

    std::vector<Point> W;
    if (pattern == "same-leaf") {
        std::vector<int> ok;
        for (int i = 0; i < n_leaves; ++i) {
            if (static_cast<int>(by_leaf[static_cast<size_t>(i)].size()) >= want) ok.push_back(i);
        }
        if (ok.empty()) {
            int best = 0;
            for (int i = 1; i < n_leaves; ++i) {
                if (by_leaf[static_cast<size_t>(i)].size() > by_leaf[static_cast<size_t>(best)].size()) best = i;
            }
            ok.push_back(best);
        }
        std::uniform_int_distribution<size_t> dist(0, ok.size() - 1);
        const int li = ok[dist(rng)];
        const auto& cells = by_leaf[static_cast<size_t>(li)];
        const int take = std::min(want, static_cast<int>(cells.size()));
        W.assign(cells.begin(), cells.begin() + take);
    } else {
        std::vector<size_t> cursor(static_cast<size_t>(n_leaves), 0);
        bool progress = true;
        while (static_cast<int>(W.size()) < want && progress) {
            progress = false;
            for (int i = 0; i < n_leaves && static_cast<int>(W.size()) < want; ++i) {
                auto& cells = by_leaf[static_cast<size_t>(i)];
                if (cursor[static_cast<size_t>(i)] < cells.size()) {
                    W.push_back(cells[cursor[static_cast<size_t>(i)]++]);
                    progress = true;
                }
            }
        }
    }
    return W;
}

} // namespace

int main(int argc, char** argv) {
    std::string map_path;
    std::string map_id;
    int b = 16;
    int g = 2;
    std::vector<int> ms = {1, 4, 16, 64};
    std::vector<std::string> patterns = {"same-leaf", "multi-leaf"};
    int locations = 10;
    int warmup = 3;
    int repeats = 10;
    int rebuild_repeats = -1;
    int seed0 = 201;
    std::string csv = "campaigns/manuscript/E3/update_latency.csv";

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << name << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--map") map_path = need("--map");
        else if (arg == "--map-id") map_id = need("--map-id");
        else if (arg == "--b") b = std::stoi(need("--b"));
        else if (arg == "--g") g = std::stoi(need("--g"));
        else if (arg == "--m") ms = parseIntList(need("--m"), ms);
        else if (arg == "--locations") locations = std::stoi(need("--locations"));
        else if (arg == "--repeats") repeats = std::stoi(need("--repeats"));
        else if (arg == "--rebuild-repeats") rebuild_repeats = std::stoi(need("--rebuild-repeats"));
        else if (arg == "--warmup") warmup = std::stoi(need("--warmup"));
        else if (arg == "--seed0") seed0 = std::stoi(need("--seed0"));
        else if (arg == "--csv") csv = need("--csv");
        else if (arg == "--pattern") {
            patterns = {need("--pattern")};
        } else if (arg == "--help") {
            std::cout << "bench_update_latency --map PATH --map-id ID [--m 1,4,16,64] [--locations 10] [--repeats 10]\n"
                      << "  [--rebuild-repeats N]  (default: same as --repeats)\n";
            return 0;
        }
    }
    if (map_path.empty()) {
        std::cerr << "--map is required\n";
        return 2;
    }
    if (map_id.empty()) map_id = "map";
    if (rebuild_repeats < 0) rebuild_repeats = repeats;
    if (rebuild_repeats < 1) rebuild_repeats = 1;

    GridMap gt;
    if (!gt.loadMovingAI(resolveExistingPath(map_path))) {
        std::cerr << "failed to load " << map_path << "\n";
        return 1;
    }
    const auto cc = largestGtComponent(gt);
    auto check_pts = sampleSector(gt, cc, 8);
    if (check_pts.size() > 25) check_pts.resize(25);

    HierarchyConfig hc;
    hc.tile_size_b = b;
    hc.group_w = g;
    hc.group_h = g;
    hc.hybrid_k_threshold = 0;
    const std::string binary_sha = fileSha256("/proc/self/exe");

    for (const auto& pattern : patterns) {
        for (int m : ms) {
            const int want_w = std::max(m, 64);
            for (int loc = 0; loc < locations; ++loc) {
                const int seed = seed0 + loc;
                auto W = chooseWithheld(pattern, gt, cc, check_pts, b, static_cast<std::uint64_t>(seed), want_w);
                if (static_cast<int>(W.size()) < m) {
                    std::cerr << "not enough withheld cells for " << map_id << ' ' << pattern << " m=" << m << "\n";
                    return 1;
                }
                std::vector<Point> B(W.begin(), W.begin() + m);

                GridMap pre_grid = gt;
                markObservedFreeExcept(pre_grid, W);
                MinPlusOracle pre(hc);
                if (!pre.build(pre_grid)) {
                    std::cerr << "pre-opening build failed\n";
                    return 1;
                }

                GridMap post_grid = pre_grid;
                for (Point p : B) post_grid.setObserved(p, CellState::Free);

                if (warmup > 0) {
                    MinPlusOracle throwaway = pre;
                    throwaway.onBatchCellsOpened(B, post_grid);
                }

                std::vector<double> upd_ms;
                int touched = 0, dirty = 0, pivots = 0;
                MinPlusOracle live_last = pre;
                for (int r = 0; r < repeats; ++r) {
                    MinPlusOracle live = pre;
                    const auto t0 = Clock::now();
                    live.onBatchCellsOpened(B, post_grid);
                    const auto t1 = Clock::now();
                    upd_ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
                    touched = live.last_batch_stats.n_touched_leaves;
                    dirty = live.last_batch_stats.n_dirty_ancestors;
                    pivots = live.last_batch_stats.n_pivots;
                    live_last = std::move(live);
                }

                std::vector<double> rebuild_ms;
                for (int r = 0; r < rebuild_repeats; ++r) {
                    MinPlusOracle rebuilt(hc);
                    const auto t0 = Clock::now();
                    if (!rebuilt.build(post_grid)) return 1;
                    const auto t1 = Clock::now();
                    rebuild_ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
                }

                std::vector<double> kbfs_ms;
                for (int r = 0; r < repeats; ++r) {
                    const auto t0 = Clock::now();
                    (void)BfsOracle::queryAllPairs(post_grid, check_pts);
                    const auto t1 = Clock::now();
                    kbfs_ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
                }

                const auto htrip = live_last.queryAllPairs(check_pts, post_grid);
                const auto bfs = BfsOracle::queryAllPairs(post_grid, check_pts);
                const MatrixCheck chk = compareMatrices(htrip, bfs);
                if (chk.mismatches != 0) {
                    std::cerr << "post-update mismatch " << map_id << " m=" << m << " loc=" << loc << "\n";
                    return 1;
                }

                const TimingStats us = summarizeMs(upd_ms);
                const TimingStats rs = summarizeMs(rebuild_ms);
                const TimingStats ks = summarizeMs(kbfs_ms);
                std::ostringstream line;
                line << map_id << ',' << b << ',' << g << ',' << pattern << ',' << m << ',' << loc << ',' << seed << ','
                     << touched << ',' << dirty << ',' << pivots << ','
                     << us.median_ms << ',' << us.min_ms << ',' << us.max_ms << ','
                     << rs.median_ms << ',' << ks.median_ms << ','
                     << check_pts.size() << ',' << chk.entries << ',' << chk.mismatches << ','
                     << peakRssBytes() << ',' << binary_sha << ",query_only";
                csvAppendLine(csv, updateHeader(), line.str());
                std::cout << map_id << ' ' << pattern << " m=" << m << " loc=" << loc
                          << " upd_median_ms=" << us.median_ms << " rebuild_median_ms=" << rs.median_ms << "\n";
            }
        }
    }
    return 0;
}
