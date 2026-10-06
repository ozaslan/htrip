// Property-based differential test for H-TRIP.
//
// For every seed the test builds a procedural workspace, runs a random sequence
// of map mutations (doorway openings and block/reopen cycles) and continuously
// checks three invariants:
//   1. full H-TRIP matrix == BFS ground truth
//   2. incremental refresh == full recomputation (bit-exact)
//   3. queryDistance(s, t) == BFS(s, t) for random pairs
//
// Usage: test_property_based [--seeds N] [--size S] [--tile b] [--group g]
//                            [--steps T] [--viewpoints K]
#include "htrip/grid_map.hpp"
#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace htrip;

namespace {

struct Options {
    int seeds = 6;
    int size = 96;
    int tile = 16;
    int group = 2;
    int steps = 24;
    int viewpoints = 24;
};

Options parseArgs(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&](int fallback) { return (i + 1 < argc) ? std::atoi(argv[++i]) : fallback; };
        if (arg == "--seeds") opt.seeds = next(opt.seeds);
        else if (arg == "--size") opt.size = next(opt.size);
        else if (arg == "--tile") opt.tile = next(opt.tile);
        else if (arg == "--group") opt.group = next(opt.group);
        else if (arg == "--steps") opt.steps = next(opt.steps);
        else if (arg == "--viewpoints") opt.viewpoints = next(opt.viewpoints);
    }
    return opt;
}

std::vector<Point> pickViewpoints(const GridMap& grid, size_t count, std::mt19937_64& rng) {
    std::vector<Point> free_cells;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            if (grid.getObserved(r, c) == CellState::Free) free_cells.push_back({r, c});
        }
    }
    std::shuffle(free_cells.begin(), free_cells.end(), rng);
    if (free_cells.size() > count) free_cells.resize(count);
    return free_cells;
}

bool matricesEqual(const std::vector<dist_t>& a, const std::vector<dist_t>& b, size_t& first_bad) {
    if (a.size() != b.size()) {
        first_bad = 0;
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            first_bad = i;
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    const Options opt = parseArgs(argc, argv);
    std::cout << "H-TRIP property-based differential test\n"
              << "  seeds=" << opt.seeds << " size=" << opt.size
              << " tile=" << opt.tile << " group=" << opt.group
              << " steps=" << opt.steps << " viewpoints=" << opt.viewpoints << "\n";

    size_t total_verifications = 0;
    size_t total_openings = 0;

    for (int s = 0; s < opt.seeds; ++s) {
        std::mt19937_64 rng(1000u + static_cast<uint64_t>(s) * 7919u);

        GridMap grid;
        grid.generateProceduralRooms(opt.size, opt.size, 16, 100u + static_cast<uint64_t>(s));
        grid.initExplorationSnapshot({2, 2}, 0.60, 0.12, 200u + static_cast<uint64_t>(s));

        const std::vector<Point> vps = pickViewpoints(grid, static_cast<size_t>(opt.viewpoints), rng);

        HierarchyConfig cfg;
        cfg.tile_size_b = opt.tile;
        cfg.group_w = opt.group;
        cfg.group_h = opt.group;
        if (!cfg.isValid()) {
            std::cerr << "FAIL: invalid configuration\n";
            return 1;
        }
        MinPlusOracle oracle(cfg);
        if (!oracle.build(grid)) {
            std::cerr << "FAIL: oracle build failed (seed " << s << ")\n";
            return 1;
        }

        MinPlusOracle::HTripQueryScratch scratch;
        MinPlusOracle::HTripIncrementalState state;
        std::vector<dist_t> full_mat;
        std::vector<dist_t> inc_mat;

        oracle.queryAllPairs(vps, grid, full_mat, scratch);
        std::vector<dist_t> bfs_mat = BfsOracle::queryAllPairs(grid, vps);
        size_t bad = 0;
        if (!matricesEqual(full_mat, bfs_mat, bad)) {
            std::cerr << "FAIL: seed " << s << " initial full vs BFS mismatch at " << bad << "\n";
            return 1;
        }
        total_verifications++;

        for (int step = 0; step < opt.steps; ++step) {
            // Random operation: open an unexplored free-truth cell, or cycle a
            // free non-viewpoint cell through blocked -> reopened.
            const int op = static_cast<int>(rng() % 3);
            if (op == 0) {
                // Block a free cell (that is not a viewpoint) then reopen it.
                const int r = static_cast<int>(rng() % opt.size);
                const int c = static_cast<int>(rng() % opt.size);
                const Point p{r, c};
                const bool is_viewpoint = std::find(vps.begin(), vps.end(), p) != vps.end();
                if (grid.getObserved(p) == CellState::Free && !is_viewpoint) {
                    grid.blockCorridor(p);
                    oracle.onCellBlocked(p, grid);
                    grid.openDoorway(p);
                    oracle.onCellOpened(p, grid);
                }
            } else {
                // Open an unexplored but truly free cell (map growth).
                for (int attempt = 0; attempt < 200; ++attempt) {
                    const int r = static_cast<int>(rng() % opt.size);
                    const int c = static_cast<int>(rng() % opt.size);
                    const Point p{r, c};
                    if (grid.ground_truth[grid.index(p)] != CellState::Free) continue;
                    if (grid.getObserved(p) == CellState::Free) continue;
                    grid.openDoorway(p);
                    oracle.onCellOpened(p, grid);
                    total_openings++;
                    break;
                }
            }

            oracle.queryAllPairsIncremental(vps, grid, inc_mat, scratch, state);
            oracle.queryAllPairs(vps, grid, full_mat, scratch);
            if (!matricesEqual(inc_mat, full_mat, bad)) {
                std::cerr << "FAIL: seed " << s << " step " << step
                          << " incremental vs full mismatch at index " << bad
                          << " (inc=" << inc_mat[bad] << " full=" << full_mat[bad] << ")\n";
                return 1;
            }
            bfs_mat = BfsOracle::queryAllPairs(grid, vps);
            if (!matricesEqual(inc_mat, bfs_mat, bad)) {
                std::cerr << "FAIL: seed " << s << " step " << step
                          << " incremental vs BFS mismatch at index " << bad
                          << " (inc=" << inc_mat[bad] << " bfs=" << bfs_mat[bad] << ")\n";
                return 1;
            }
            total_verifications++;

            // Spot-check the single-pair API against BFS.
            if (!vps.empty()) {
                const size_t i = static_cast<size_t>(rng() % vps.size());
                const size_t j = static_cast<size_t>(rng() % vps.size());
                const dist_t single = oracle.queryDistance(vps[i], vps[j], grid);
                const dist_t truth = BfsOracle::queryDistance(grid, vps[i], vps[j]);
                if (single != truth) {
                    std::cerr << "FAIL: seed " << s << " step " << step
                              << " queryDistance mismatch: " << single << " vs " << truth << "\n";
                    return 1;
                }
            }
        }
        std::cout << "  seed " << s << ": OK\n";
    }

    std::cout << "ALL PROPERTY-BASED CHECKS PASSED (verifications=" << total_verifications
              << ", openings=" << total_openings << ")\n";
    return 0;
}
