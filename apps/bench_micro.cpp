// Google Benchmark micro/macro benchmarks for the H-TRIP oracle.
//
// Run:  ./build/dev/bench_micro
//       ./build/dev/bench_micro --benchmark_filter=Query
#include <benchmark/benchmark.h>

#include "htrip/bfs_oracle.hpp"
#include "htrip/grid_map.hpp"
#include "htrip/minplus_oracle.hpp"

#include <algorithm>
#include <memory>
#include <random>
#include <vector>

using namespace htrip;

namespace {

constexpr int kW = 128; // workspace side for the macro benchmarks

GridMap makeMap(uint64_t seed) {
    GridMap grid;
    grid.generateProceduralRooms(kW, kW, 16, seed);
    grid.initExplorationSnapshot({2, 2}, 0.60, 0.10, seed + 1);
    return grid;
}

std::vector<Point> pickViewpoints(const GridMap& grid, size_t k) {
    std::vector<Point> cells;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            if (grid.getObserved(r, c) == CellState::Free) cells.push_back({r, c});
        }
    }
    std::mt19937_64 rng(12345);
    std::shuffle(cells.begin(), cells.end(), rng);
    if (cells.size() > k) cells.resize(k);
    return cells;
}

std::vector<Point> pickCandidates(const GridMap& grid, uint64_t seed) {
    std::vector<Point> cells;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const Point p{r, c};
            if (grid.ground_truth[grid.index(p)] == CellState::Free &&
                grid.getObserved(p) != CellState::Free) {
                cells.push_back(p);
            }
        }
    }
    std::mt19937_64 rng(seed);
    std::shuffle(cells.begin(), cells.end(), rng);
    return cells;
}

// Index build cost for representative (b, g) settings.
void BM_Build(benchmark::State& state) {
    const int b = static_cast<int>(state.range(0));
    const int g = static_cast<int>(state.range(1));
    const GridMap grid = makeMap(42);
    HierarchyConfig cfg;
    cfg.tile_size_b = b;
    cfg.group_w = g;
    cfg.group_h = g;

    for (auto _ : state) {
        MinPlusOracle oracle(cfg);
        bool ok = oracle.build(grid);
        benchmark::DoNotOptimize(ok);
    }
    state.SetLabel("b=" + std::to_string(b) + ",g=" + std::to_string(g));
}

class QueryFixture : public benchmark::Fixture {
public:
    void SetUp(::benchmark::State& state) override {
        grid_ = std::make_shared<GridMap>(makeMap(42));
        vps_ = pickViewpoints(*grid_, static_cast<size_t>(state.range(0)));
        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 2;
        cfg.group_h = 2;
        oracle_ = std::make_shared<MinPlusOracle>(cfg);
        oracle_->build(*grid_);
    }

    std::shared_ptr<GridMap> grid_;
    std::shared_ptr<MinPlusOracle> oracle_;
    std::vector<Point> vps_;
    MinPlusOracle::HTripQueryScratch scratch_;
};

BENCHMARK_DEFINE_F(QueryFixture, Full)(benchmark::State& state) {
    std::vector<dist_t> matrix;
    for (auto _ : state) {
        oracle_->queryAllPairs(vps_, *grid_, matrix, scratch_);
        benchmark::DoNotOptimize(matrix.data());
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(vps_.size() * vps_.size()));
}

BENCHMARK_DEFINE_F(QueryFixture, IncrementalNoChange)(benchmark::State& state) {
    MinPlusOracle::HTripIncrementalState inc_state;
    std::vector<dist_t> matrix;
    oracle_->queryAllPairsIncremental(vps_, *grid_, matrix, scratch_, inc_state); // warm start
    for (auto _ : state) {
        oracle_->queryAllPairsIncremental(vps_, *grid_, matrix, scratch_, inc_state);
        benchmark::DoNotOptimize(matrix.data());
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(vps_.size() * vps_.size()));
}

BENCHMARK_DEFINE_F(QueryFixture, BfsBaseline)(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<dist_t> matrix = BfsOracle::queryAllPairs(*grid_, vps_);
        benchmark::DoNotOptimize(matrix.data());
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(vps_.size() * vps_.size()));
}

BENCHMARK_DEFINE_F(QueryFixture, QueryDistance)(benchmark::State& state) {
    std::mt19937_64 rng(7);
    const size_t k = vps_.size();
    for (auto _ : state) {
        const size_t i = static_cast<size_t>(rng() % k);
        const size_t j = static_cast<size_t>(rng() % k);
        dist_t d = oracle_->queryDistance(vps_[i], vps_[j], *grid_);
        benchmark::DoNotOptimize(d);
    }
}

// Per-scan workload: open one door and refresh the k x k matrix incrementally,
// compared against a from-scratch BFS refresh.
class CycleFixture : public benchmark::Fixture {
public:
    void SetUp(::benchmark::State& state) override {
        grid_ = std::make_shared<GridMap>(makeMap(42));
        vps_ = pickViewpoints(*grid_, static_cast<size_t>(state.range(0)));
        HierarchyConfig cfg;
        cfg.tile_size_b = 16;
        cfg.group_w = 2;
        cfg.group_h = 2;
        oracle_ = std::make_shared<MinPlusOracle>(cfg);
        oracle_->build(*grid_);
        candidates_ = pickCandidates(*grid_, 99);
        cursor_ = 0;
    }

    std::shared_ptr<GridMap> grid_;
    std::shared_ptr<MinPlusOracle> oracle_;
    std::vector<Point> vps_;
    std::vector<Point> candidates_;
    size_t cursor_ = 0;
    MinPlusOracle::HTripQueryScratch scratch_;
    MinPlusOracle::HTripIncrementalState inc_state_;
};

BENCHMARK_DEFINE_F(CycleFixture, OpenAndRefresh)(benchmark::State& state) {
    std::vector<dist_t> matrix;
    for (auto _ : state) {
        if (cursor_ >= candidates_.size()) break;
        const Point p = candidates_[cursor_++];
        grid_->openDoorway(p);
        oracle_->onCellOpened(p, *grid_);
        oracle_->queryAllPairsIncremental(vps_, *grid_, matrix, scratch_, inc_state_);
        benchmark::DoNotOptimize(matrix.data());
    }
    state.counters["updates"] = static_cast<double>(cursor_);
}

BENCHMARK_DEFINE_F(CycleFixture, OpenAndBfsRefresh)(benchmark::State& state) {
    std::vector<dist_t> matrix;
    for (auto _ : state) {
        if (cursor_ >= candidates_.size()) break;
        const Point p = candidates_[cursor_++];
        grid_->openDoorway(p);
        matrix = BfsOracle::queryAllPairs(*grid_, vps_);
        benchmark::DoNotOptimize(matrix.data());
    }
    state.counters["updates"] = static_cast<double>(cursor_);
}

} // namespace

BENCHMARK(BM_Build)
    ->Args({8, 2})
    ->Args({16, 2})
    ->Args({16, 4})
    ->Args({32, 2})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(QueryFixture, Full)
    ->Arg(10)->Arg(25)->Arg(50)->Arg(100)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(QueryFixture, IncrementalNoChange)
    ->Arg(10)->Arg(25)->Arg(50)->Arg(100)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(QueryFixture, BfsBaseline)
    ->Arg(10)->Arg(25)->Arg(50)->Arg(100)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(QueryFixture, QueryDistance)
    ->Arg(50)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_REGISTER_F(CycleFixture, OpenAndRefresh)
    ->Arg(25)->Arg(50)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(CycleFixture, OpenAndBfsRefresh)
    ->Arg(25)->Arg(50)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
