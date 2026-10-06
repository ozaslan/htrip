#include "htrip/astar_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include "htrip/campaign_support.hpp"
#include "htrip/minplus_oracle.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace htrip;
using namespace htrip::campaign;

namespace {

using Clock = std::chrono::steady_clock;

bool contains(const std::vector<int>& v, int x) {
    return std::find(v.begin(), v.end(), x) != v.end();
}

std::vector<double> timeMethod(const std::string& method, MinPlusOracle& oracle, const GridMap& grid,
                               const std::vector<Point>& pts, int warmup, int repeats) {
    auto once = [&]() {
        if (method == "HTRIP") {
            (void)oracle.queryAllPairs(pts, grid);
        } else if (method == "KBFS") {
            (void)BfsOracle::queryAllPairs(grid, pts);
        } else {
            (void)AStarOracle::queryAllPairs(grid, pts);
        }
    };
    for (int i = 0; i < warmup; ++i) once();
    std::vector<double> ms;
    ms.reserve(static_cast<size_t>(repeats));
    for (int i = 0; i < repeats; ++i) {
        const auto t0 = Clock::now();
        once();
        const auto t1 = Clock::now();
        ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return ms;
}

} // namespace

int main(int argc, char** argv) {
    std::string map_path;
    std::string map_id;
    int b = 16;
    int g = 2;
    std::string geometry = "sector";
    std::uint64_t seed = 101;
    std::vector<int> ks = {1, 2, 4, 8, 16, 25, 32, 40, 64};
    std::vector<int> astar_at = {4, 25};
    std::vector<std::string> methods = {"HTRIP", "KBFS", "ASTAR"};
    int warmup = 3;
    int repeats = 10;
    std::string csv = "campaigns/manuscript/E2/query_vs_k.csv";
    std::string repeats_csv = "campaigns/manuscript/E2/query_vs_k_repeats.csv";
    std::string rss_csv = "campaigns/manuscript/E4/rss.csv";
    std::string build_csv = "campaigns/manuscript/E7/full_build.csv";
    bool write_rss = true;
    bool write_build = true;
    bool build_only = false;
    int set_id = 0;
    int origin_r = -1;
    int origin_c = -1;

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
        else if (arg == "--geometry") geometry = need("--geometry");
        else if (arg == "--seed") seed = static_cast<std::uint64_t>(std::stoull(need("--seed")));
        else if (arg == "--set-id") set_id = std::stoi(need("--set-id"));
        else if (arg == "--k") ks = parseIntList(need("--k"), ks);
        else if (arg == "--astar-at") astar_at = parseIntList(need("--astar-at"), astar_at);
        else if (arg == "--warmup") warmup = std::stoi(need("--warmup"));
        else if (arg == "--repeats") repeats = std::stoi(need("--repeats"));
        else if (arg == "--csv") csv = need("--csv");
        else if (arg == "--repeats-csv") repeats_csv = need("--repeats-csv");
        else if (arg == "--rss-csv") rss_csv = need("--rss-csv");
        else if (arg == "--build-csv") build_csv = need("--build-csv");
        else if (arg == "--no-rss") write_rss = false;
        else if (arg == "--no-build") write_build = false;
        else if (arg == "--build-only") build_only = true;
        else if (arg == "--origin-r") origin_r = std::stoi(need("--origin-r"));
        else if (arg == "--origin-c") origin_c = std::stoi(need("--origin-c"));
        else if (arg == "--methods") {
            methods.clear();
            std::stringstream ss(need("--methods"));
            std::string tok;
            while (std::getline(ss, tok, ',')) if (!tok.empty()) methods.push_back(tok);
        } else if (arg == "--paper-e5") {
            ks = {25};
            geometry = "sector";
        } else if (arg == "--help") {
            std::cout << "bench_query_vs_k --map PATH --map-id ID [--b 16] [--g 2] [--geometry sector|spread|cluster]\n"
                      << "  [--k 1,2,4,8,16,25,32,40,64,128,256] [--set-id N] [--seed S]\n"
                      << "  [--warmup 3] [--repeats 10] [--csv PATH] [--build-csv PATH] [--build-only]\n";
            return 0;
        }
    }

    if (map_path.empty()) {
        std::cerr << "--map is required\n";
        return 2;
    }
    if (map_id.empty()) map_id = std::filesystem::path(map_path).stem().string();

    GridMap grid;
    if (!loadFullyRevealed(map_path, grid)) {
        std::cerr << "failed to load " << map_path << "\n";
        return 1;
    }
    const auto cc = largestGtComponent(grid);
    size_t n_need = 0;
    for (int k : ks) {
        if (k > 0) n_need = std::max(n_need, static_cast<size_t>(k));
    }
    if (!build_only && n_need == 0) {
        std::cerr << "no positive k\n";
        return 1;
    }
    if (n_need == 0) n_need = 1;
    if (origin_r < 0 || origin_c < 0) {
        origin_r = 0;
        origin_c = 0;
        if (geometry == "sector" && set_id > 0) {
            const int divisions = 8;
            const int step_r = std::max(1, grid.rows / divisions);
            const int step_c = std::max(1, grid.cols / divisions);
            const int ir = set_id % 4;
            const int ic = (set_id / 4) % 4;
            origin_r = ir * std::max(1, step_r / 4);
            origin_c = ic * std::max(1, step_c / 4);
        }
    }
    auto pool = sampleGeometry(geometry, grid, cc, seed, n_need, origin_r, origin_c);
    if (pool.size() < n_need) {
        std::cerr << "sampled " << pool.size() << " points, need " << n_need
                  << " (largest CC size " << cc.size << ")\n";
        return 1;
    }
    std::cout << map_id << " sampled " << pool.size() << " points (CC " << cc.size
              << ") origin=(" << origin_r << ',' << origin_c << ") set_id=" << set_id << "\n";

    HierarchyConfig hc;
    hc.tile_size_b = b;
    hc.group_w = g;
    hc.group_h = g;
    hc.hybrid_k_threshold = 0;
    MinPlusOracle oracle(hc);
    const auto build_t0 = Clock::now();
    if (!oracle.build(grid)) {
        std::cerr << "oracle build failed\n";
        return 1;
    }
    const auto build_t1 = Clock::now();
    const double build_ms = std::chrono::duration<double, std::milli>(build_t1 - build_t0).count();
    const std::size_t index_bytes = oracle.totalMemoryBytes();
    const std::string binary_sha = fileSha256("/proc/self/exe");
    std::cout << map_id << " full_build_ms=" << build_ms << " index_bytes=" << index_bytes << "\n";
    if (write_build) {
        std::ostringstream bl;
        bl << map_id << ',' << b << ',' << g << ',' << set_id << ',' << seed << ','
           << build_ms << ',' << index_bytes << ',' << peakRssBytes() << ','
           << binary_sha << ",query_only";
        csvAppendLine(build_csv, fullBuildHeader(), bl.str());
    }
    if (build_only) return 0;

    for (int k : ks) {
        if (k <= 0) continue;
        const size_t kk = std::min(static_cast<size_t>(k), pool.size());
        std::vector<Point> pts(pool.begin(), pool.begin() + static_cast<std::ptrdiff_t>(kk));
        const auto htrip_mat = oracle.queryAllPairs(pts, grid);
        const auto bfs_mat = BfsOracle::queryAllPairs(grid, pts);
        const MatrixCheck chk = compareMatrices(htrip_mat, bfs_mat);
        if (chk.mismatches != 0) {
            std::cerr << "exactness mismatch at k=" << k << " mismatches=" << chk.mismatches << "\n";
            return 1;
        }
        for (const auto& method : methods) {
            if (method == "ASTAR" && !contains(astar_at, k)) continue;
            auto samples = timeMethod(method, oracle, grid, pts, warmup, repeats);
            const TimingStats st = summarizeMs(samples);
            std::ostringstream line;
            line << map_id << ',' << b << ',' << g << ',' << geometry << ',' << set_id << ',' << seed << ','
                 << k << ',' << method << ',' << warmup << ',' << repeats << ','
                 << st.median_ms << ',' << st.mean_ms << ',' << st.min_ms << ',' << st.max_ms << ','
                 << (kk * kk) << ',' << index_bytes << ',' << peakRssBytes() << ','
                 << chk.max_finite << ',' << chk.mismatches << ',' << binary_sha << ",query_only";
            csvAppendLine(csv, queryVsKHeader(), line.str());
            for (int r = 0; r < static_cast<int>(samples.size()); ++r) {
                std::ostringstream rl;
                rl << map_id << ',' << b << ',' << g << ',' << geometry << ',' << set_id << ',' << seed << ','
                   << k << ',' << method << ',' << r << ',' << samples[static_cast<size_t>(r)] << ",query_only";
                csvAppendLine(repeats_csv, queryVsKRepeatsHeader(), rl.str());
            }
            std::cout << map_id << " k=" << k << ' ' << method << " median_ms=" << st.median_ms << "\n";
        }
        if (write_rss && k == 25 && geometry == "sector") {
            std::ostringstream rss;
            rss << map_id << ',' << b << ',' << g << ",query_only,1,25,"
                << index_bytes << ',' << peakRssBytes() << ",single_oracle," << binary_sha;
            csvAppendLine(rss_csv, rssHeader(), rss.str());
        }
    }
    return 0;
}
