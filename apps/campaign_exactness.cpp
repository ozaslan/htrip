#include "htrip/campaign_support.hpp"
#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace htrip;
using namespace htrip::campaign;

namespace {

struct SuiteAcc {
    std::string suite_id;
    std::string family;
    std::string config;
    int b = 0;
    int g = 0;
    int seed = 0;
    int n_openings = 0;
    int n_blockages = 0;
    std::int64_t n_entries = 0;
    std::int64_t mismatches = 0;
    int max_abs_diff = 0;
    dist_t max_finite = 0;
    std::string reference = "BFS";
};

void flushSuite(std::vector<SuiteAcc>& rows, SuiteAcc& cur, bool& open) {
    if (!open) return;
    rows.push_back(cur);
    cur = SuiteAcc{};
    open = false;
}

void parseExactnessLog(const std::string& text, const std::string& family, std::vector<SuiteAcc>& rows,
                       int default_seed) {
    std::istringstream in(text);
    std::string line;
    SuiteAcc cur;
    bool open = false;
    const std::regex suite_re(R"(\[Suite[ ]*([^:\]]+))");
    const std::regex query_re(R"(Queries:[ ]*([0-9]+),[ ]*Matches:[ ]*([0-9]+))");
    const std::regex maxdiff_re(R"(MaxDiff:[ ]*([0-9]+))");
    const std::regex prop_re(R"(verifications=([0-9]+),[ ]*openings=([0-9]+))");
    while (std::getline(in, line)) {
        std::smatch m;
        if (std::regex_search(line, m, suite_re)) {
            flushSuite(rows, cur, open);
            cur.suite_id = m[1].str();
            while (!cur.suite_id.empty() && cur.suite_id.back() == ' ') cur.suite_id.pop_back();
            cur.family = family;
            cur.seed = default_seed;
            cur.reference = (family == "weighted_named") ? "Dijkstra" : "BFS";
            open = true;
            continue;
        }
        if (std::regex_search(line, m, query_re)) {
            if (!open) {
                cur.family = family;
                cur.suite_id = "ungrouped";
                cur.seed = default_seed;
                cur.reference = (family == "weighted_named") ? "Dijkstra" : "BFS";
                open = true;
            }
            const std::int64_t q = std::stoll(m[1].str());
            const std::int64_t matches = std::stoll(m[2].str());
            cur.n_entries += q;
            cur.mismatches += (q - matches);
            std::smatch d;
            if (std::regex_search(line, d, maxdiff_re)) {
                cur.max_abs_diff = std::max(cur.max_abs_diff, std::stoi(d[1].str()));
            }
            continue;
        }
        if (std::regex_search(line, m, prop_re)) {
            flushSuite(rows, cur, open);
            SuiteAcc p;
            p.suite_id = "property";
            p.family = "property";
            p.config = "seeds/size/steps from CLI";
            p.seed = default_seed;
            p.n_openings = std::stoi(m[2].str());
            const std::int64_t ver = std::stoll(m[1].str());
            p.n_entries = ver; // each verification is one matrix compare; entry count filled by --viewpoints below
            p.reference = "BFS";
            rows.push_back(p);
        }
    }
    flushSuite(rows, cur, open);
}

std::string runCapture(const std::string& cmd, const std::string& log_path, int& rc) {
    std::filesystem::path lp(log_path);
    if (lp.has_parent_path()) std::filesystem::create_directories(lp.parent_path());
    const std::string wrapped = cmd + " > \"" + log_path + "\" 2>&1";
    rc = std::system(wrapped.c_str());
    std::ifstream in(log_path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void writeRows(const std::string& csv, const std::vector<SuiteAcc>& rows,
               const std::string& binary_sha, const std::string& git_sha) {
    for (const auto& r : rows) {
        std::ostringstream line;
        line << r.suite_id << ',' << r.family << ',' << r.config << ','
             << r.b << ',' << r.g << ',' << r.seed << ','
             << r.n_openings << ',' << r.n_blockages << ','
             << r.n_entries << ',' << r.reference << ','
             << r.mismatches << ',' << r.max_abs_diff << ','
             << r.max_finite << ',' << binary_sha << ',' << git_sha;
        csvAppendLine(csv, exactnessHeader(), line.str());
    }
}

} // namespace

int main(int argc, char** argv) {
    std::string csv = "campaigns/manuscript/E1/exactness_summary.csv";
    std::string bin_dir;
    std::string stdout_dir = "campaigns/manuscript/E1/stdout";
    std::string maps_json = "campaigns/manuscript/maps.json";
    bool run_tests = true;
    bool map_static = true;
    int seed = 12345;
    int prop_viewpoints = 24;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << name << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--csv") csv = need("--csv");
        else if (arg == "--bin-dir") bin_dir = need("--bin-dir");
        else if (arg == "--stdout-dir") stdout_dir = need("--stdout-dir");
        else if (arg == "--maps-json") maps_json = need("--maps-json");
        else if (arg == "--seed") seed = std::stoi(need("--seed"));
        else if (arg == "--no-tests") run_tests = false;
        else if (arg == "--no-map-static") map_static = false;
        else if (arg == "--help") {
            std::cout << "campaign_exactness [--csv PATH] [--bin-dir DIR] [--no-tests] [--no-map-static]\n";
            return 0;
        }
    }

    const std::string repo = findRepoRoot();
    if (bin_dir.empty()) bin_dir = (std::filesystem::path(repo) / "build/dev").string();
    const std::string git_sha = gitSha(repo);
    const std::string exact_bin = (std::filesystem::path(bin_dir) / "test_exactness").string();
    const std::string binary_sha = fileSha256(exact_bin);

    std::filesystem::create_directories(std::filesystem::path(csv).parent_path());
    std::filesystem::create_directories(stdout_dir);

    if (run_tests) {
        struct Job {
            const char* name;
            std::string cmd;
            const char* family;
        };
        const std::string pbin = (std::filesystem::path(bin_dir) / "test_property_based").string();
        const std::string wbin = (std::filesystem::path(bin_dir) / "test_weighted").string();
        const std::string rbin = (std::filesystem::path(bin_dir) / "test_reproducers").string();
        const Job jobs[] = {
            {"test_exactness", exact_bin + " " + std::to_string(seed), "unweighted_named"},
            {"test_weighted", wbin, "weighted_named"},
            {"test_property_based",
             pbin + " --seeds 10 --size 96 --tile 16 --group 2 --steps 24 --viewpoints " +
                 std::to_string(prop_viewpoints),
             "property"},
            {"test_reproducers", rbin, "reproducer"},
        };
        for (const auto& job : jobs) {
            if (!std::filesystem::exists(std::filesystem::path(bin_dir) / job.name)) {
                std::cerr << "skip missing " << job.name << "\n";
                continue;
            }
            int rc = 0;
            const std::string log = (std::filesystem::path(stdout_dir) / (std::string(job.name) + ".log")).string();
            std::cout << "running " << job.name << " ...\n" << std::flush;
            const std::string text = runCapture(job.cmd, log, rc);
            std::vector<SuiteAcc> rows;
            parseExactnessLog(text, job.family, rows, seed);
            if (std::string(job.family) == "property") {
                for (auto& r : rows) {
                    if (r.family == "property") {
                        r.config = "seeds=10;size=96;b=16;g=2;steps=24;k=24";
                        r.b = 16;
                        r.g = 2;
                        r.n_entries *= static_cast<std::int64_t>(prop_viewpoints) * prop_viewpoints;
                    }
                }
            }
            if (std::string(job.family) == "reproducer" && rows.empty()) {
                SuiteAcc r;
                r.suite_id = "reproducers";
                r.family = "reproducer";
                r.config = "bug_repro_suite";
                r.reference = "BFS";
                r.mismatches = (rc == 0 ? 0 : 1);
                rows.push_back(r);
            }
            if (rc != 0) {
                std::cerr << job.name << " failed, rc=" << rc << " log=" << log << "\n";
                if (rows.empty()) {
                    SuiteAcc fail;
                    fail.suite_id = job.name;
                    fail.family = job.family;
                    fail.mismatches = 1;
                    fail.reference = "BFS";
                    rows.push_back(fail);
                }
                for (auto& r : rows) {
                    if (r.mismatches == 0) r.mismatches = 1;
                }
            }
            writeRows(csv, rows, fileSha256((std::filesystem::path(bin_dir) / job.name).string()), git_sha);
        }
    }

    if (map_static) {
        std::ostringstream maps;
        maps << "[\n";
        bool first_map = true;
        for (const auto& pm : kPaperMaps) {
            GridMap grid;
            if (!loadFullyRevealed(pm.path, grid)) {
                std::cerr << "failed to load " << pm.path << "\n";
                return 1;
            }
            const auto cc = largestGtComponent(grid);
            auto vps = sampleSector(grid, cc, 8);
            if (vps.size() > 25) vps.resize(25);
            HierarchyConfig hc;
            hc.tile_size_b = pm.b;
            hc.group_w = pm.g;
            hc.group_h = pm.g;
            hc.hybrid_k_threshold = 0;
            MinPlusOracle oracle(hc);
            if (!oracle.build(grid)) {
                std::cerr << "build failed for " << pm.id << "\n";
                return 1;
            }
            const auto htrip = oracle.queryAllPairs(vps, grid);
            const auto bfs = BfsOracle::queryAllPairs(grid, vps);
            const MatrixCheck chk = compareMatrices(htrip, bfs);
            SuiteAcc r;
            r.suite_id = std::string(pm.id) + "_k25";
            r.family = "map_static_k25";
            r.config = pm.path;
            r.b = pm.b;
            r.g = pm.g;
            r.seed = 0;
            r.n_entries = static_cast<std::int64_t>(chk.entries);
            r.mismatches = static_cast<std::int64_t>(chk.mismatches);
            r.max_abs_diff = chk.max_abs_diff;
            r.max_finite = chk.max_finite;
            r.reference = "BFS";
            writeRows(csv, {r}, binary_sha, git_sha);

            if (!first_map) maps << ",\n";
            first_map = false;
            maps << "  {\"id\": \"" << pm.id << "\", \"path\": \"" << jsonEscape(resolveExistingPath(pm.path))
                 << "\", \"sha256\": \"" << fileSha256(resolveExistingPath(pm.path))
                 << "\", \"rows\": " << grid.rows << ", \"cols\": " << grid.cols
                 << ", \"n_passable\": " << countGtPassable(grid)
                 << ", \"n_largest_cc\": " << cc.size
                 << ", \"max_finite_geodesic_on_sample\": " << chk.max_finite
                 << ", \"b\": " << pm.b << ", \"g\": " << pm.g << "}";
            std::cout << pm.id << " k=25 entries=" << chk.entries << " mismatches=" << chk.mismatches << "\n";
        }
        maps << "\n]\n";
        std::filesystem::create_directories(std::filesystem::path(maps_json).parent_path());
        std::ofstream mj(maps_json);
        mj << maps.str();
    }

    std::cout << "wrote " << csv << "\n";
    return 0;
}
