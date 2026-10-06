#include "htrip/campaign_support.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace htrip::campaign;

int main(int argc, char** argv) {
    std::string out_dir = "campaigns/manuscript/E0";
    std::string bin_dir;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << name << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--out") out_dir = need("--out");
        else if (arg == "--bin-dir") bin_dir = need("--bin-dir");
        else if (arg == "--help") {
            std::cout << "campaign_e0 --out DIR [--bin-dir DIR]\n";
            return 0;
        }
    }

    const std::string repo = findRepoRoot();
    if (bin_dir.empty()) {
        bin_dir = (std::filesystem::path(repo) / "build/dev").string();
    }

    std::filesystem::create_directories(out_dir);

    const std::string hostname = shellCapture("hostname");
    std::string cpu = shellCapture("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //'");
    if (cpu.empty()) cpu = "unknown";
    const std::string physical = shellCapture("lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l | tr -d '[:space:]'");
    const std::string logical = shellCapture("nproc");
    const std::string ram_kb = shellCapture("awk '/MemTotal/ {print $2}' /proc/meminfo");
    double ram_gb = 0.0;
    if (!ram_kb.empty()) ram_gb = std::stod(ram_kb) / (1024.0 * 1024.0);
    const std::string os = shellCapture("uname -a");
    std::string distro = shellCapture("grep PRETTY_NAME /etc/os-release 2>/dev/null | cut -d= -f2- | tr -d '\"'");
    const std::string governor = shellCapture("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null");
    const std::string smt = shellCapture("cat /sys/devices/system/cpu/smt/active 2>/dev/null");

    {
        std::ofstream js(std::filesystem::path(out_dir) / "machine.json");
        js << "{\n"
           << "  \"hostname\": \"" << jsonEscape(hostname) << "\",\n"
           << "  \"cpu_model\": \"" << jsonEscape(cpu) << "\",\n"
           << "  \"physical_cores\": " << (physical.empty() ? "null" : physical) << ",\n"
           << "  \"logical_cpus\": " << (logical.empty() ? "null" : logical) << ",\n"
           << "  \"installed_ram_gb\": " << ram_gb << ",\n"
           << "  \"os\": \"" << jsonEscape(os) << "\",\n"
           << "  \"distro\": \"" << jsonEscape(distro) << "\",\n"
           << "  \"governor\": \"" << jsonEscape(governor.empty() ? "unknown" : governor) << "\",\n"
           << "  \"turbo\": \"unknown\",\n"
           << "  \"pinned\": \"none\",\n"
           << "  \"smt\": \"" << jsonEscape(smt.empty() ? "unknown" : smt) << "\"\n"
           << "}\n";
    }

    const std::vector<std::string> bins = {
        "test_exactness", "test_weighted", "test_property_based", "test_reproducers",
        "bench_query_vs_k", "bench_update_latency", "campaign_exactness", "campaign_e0",
        "bench_systematic_bfs", "run_multi_robot_sim"
    };

    std::ofstream bj(std::filesystem::path(out_dir) / "build.json");
    bj << "{\n"
       << "  \"git_sha\": \"" << jsonEscape(gitSha(repo)) << "\",\n"
       << "  \"git_dirty\": " << (gitDirty(repo) ? "true" : "false") << ",\n"
       << "  \"compiler\": \"" << jsonEscape(shellCapture("c++ --version | head -n 1")) << "\",\n"
       << "  \"flags\": \"-O3 -march=native\",\n"
       << "  \"avx2\": true,\n"
       << "  \"bin_dir\": \"" << jsonEscape(bin_dir) << "\",\n"
       << "  \"executables\": {\n";
    bool first = true;
    for (const auto& name : bins) {
        const std::string path = (std::filesystem::path(bin_dir) / name).string();
        if (!std::filesystem::exists(path)) continue;
        if (!first) bj << ",\n";
        first = false;
        bj << "    \"" << name << "\": {\"path\": \"" << jsonEscape(path)
           << "\", \"sha256\": \"" << fileSha256(path) << "\"}";
    }
    bj << "\n  }\n}\n";

    std::cout << "wrote " << out_dir << "/machine.json and build.json\n";
    return 0;
}
