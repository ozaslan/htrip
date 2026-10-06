#include "htrip/trial_manifest.hpp"

#include <array>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unistd.h>

#ifdef __AVX2__
constexpr bool kAvx2Built = true;
#else
constexpr bool kAvx2Built = false;
#endif

namespace htrip {

namespace {

std::string readPopen(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (pipe == nullptr) return {};
    std::string out;
    char buf[512];
    while (fgets(buf, sizeof(buf), pipe) != nullptr) {
        out.append(buf);
    }
    pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) {
        out.pop_back();
    }
    return out;
}

void writeBool(std::ostream& out, const char* key, bool v, bool last = false) {
    out << "    \"" << key << "\": " << (v ? "true" : "false");
    out << (last ? "\n" : ",\n");
}

void writeStr(std::ostream& out, const char* key, const std::string& v, bool last = false) {
    out << "    \"" << key << "\": \"" << jsonEscape(v) << "\"";
    out << (last ? "\n" : ",\n");
}

}  // namespace

std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (unsigned char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) {
                    char hex[8];
                    std::snprintf(hex, sizeof(hex), "\\u%04x", c);
                    o += hex;
                } else {
                    o += static_cast<char>(c);
                }
        }
    }
    return o;
}

std::string captureGitCommit() {
    return readPopen("git rev-parse HEAD 2>/dev/null");
}

bool captureGitDirty() {
    const std::string porcelain = readPopen("git status --porcelain 2>/dev/null");
    return !porcelain.empty();
}

std::string captureHostname() {
    std::array<char, 256> buf{};
    if (gethostname(buf.data(), buf.size() - 1) == 0) {
        return std::string(buf.data());
    }
    return {};
}

std::string captureCpuModel() {
    return readPopen("sh -c \"grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'\"");
}

std::string sha256File(const std::string& path) {
    if (path.empty()) return {};
    std::string cmd = "sha256sum -- \"" + path + "\" 2>/dev/null | awk '{print $1}'";
    return readPopen(cmd.c_str());
}

std::string joinCommandLine(int argc, char** argv) {
    std::ostringstream oss;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) oss << ' ';
        oss << argv[i];
    }
    return oss.str();
}

VideoQuality parseVideoQuality(const std::string& name) {
    if (name == "low" || name == "0") return VideoQuality::Low;
    if (name == "high" || name == "2") return VideoQuality::High;
    return VideoQuality::Medium;
}

void fillRuntimeEnvironment(TrialManifest& m) {
    m.git_commit = captureGitCommit();
    m.git_dirty = captureGitDirty();
    m.hostname = captureHostname();
    m.cpu_model = captureCpuModel();
    m.avx2 = kAvx2Built;
}

void fillFromSimulator(TrialManifest& m, const MultiRobotSimulator& sim) {
    m.seed = sim.getSeed();
    m.num_robots = static_cast<int>(sim.robots.size());
    m.max_frontiers = sim.max_frontiers_eval;
    m.scale = sim.scale_factor;
    m.cell_size_cm = sim.cell_size_m * 100.0;
    m.tile_size_b = sim.oracle.config.tile_size_b;
    m.group_size_g = sim.oracle.config.group_w;
    m.lidar_range_m = sim.robot_config.lidar.range_m;
    m.lidar_span_deg = sim.robot_config.lidar.span_deg;
    m.lidar_num_rays = sim.robot_config.lidar.num_rays;
    m.v_max = sim.robot_config.v_max;
    m.omega_max = sim.robot_config.omega_max;
    m.radius_m = sim.robot_config.radius_m;
    m.replan_period_ticks = sim.replan_period_ticks;
    m.persistence_window_ticks = sim.persistence_window_ticks;
    m.enable_async = sim.enable_async_updates;
    m.enable_auto_respawn = sim.enable_auto_respawn;
    m.auto_respawn_ticks = sim.auto_respawn_ticks;
    m.enable_live_frontier_projection = sim.enable_live_frontier_projection;
    m.enable_bfs = sim.enable_bfs_baseline;
    m.bfs_compare_interval = sim.bfs_compare_interval;
    m.total_ticks = sim.total_ticks;
    m.final_coverage = sim.explorationFraction();
    m.grid_rows = sim.observed.rows;
    m.grid_cols = sim.observed.cols;
    m.world_w_m = sim.world_w_m;
    m.world_h_m = sim.world_h_m;
    m.tree_depth = sim.oracle.depth() + 1;
    m.total_leaves = static_cast<int>(sim.oracle.leaves.size());
    m.bfs_compare_count = sim.bfs_compare_samples.size();
    double bfs_sum = 0.0;
    double htrip_sum = 0.0;
    double speedup_sum = 0.0;
    for (const auto& sample : sim.bfs_compare_samples) {
        bfs_sum += sample.bfs_query_us;
        htrip_sum += sample.htrip_query_us;
        speedup_sum += sample.speedup;
    }
    if (!sim.bfs_compare_samples.empty()) {
        const double n = static_cast<double>(sim.bfs_compare_samples.size());
        m.bfs_query_mean_us = bfs_sum / n;
        m.paired_htrip_mean_us = htrip_sum / n;
        m.mean_speedup_vs_bfs = speedup_sum / n;
    }
}

bool writeTrialManifest(const std::string& path, const TrialManifest& m) {
    std::ofstream out(path);
    if (!out.is_open()) return false;

    out << std::boolalpha;
    out << "{\n";
    out << "  \"status\": \"" << jsonEscape(m.status) << "\",\n";
    out << "  \"stop_reason\": \"" << jsonEscape(m.stop_reason) << "\",\n";
    out << "  \"command_line\": \"" << jsonEscape(m.command_line) << "\",\n";
    out << "  \"log_dir\": \"" << jsonEscape(m.log_dir) << "\",\n";
    out << "  \"environment\": {\n";
    writeStr(out, "git_commit", m.git_commit);
    writeBool(out, "git_dirty", m.git_dirty);
    writeStr(out, "hostname", m.hostname);
    writeStr(out, "cpu_model", m.cpu_model);
    writeBool(out, "avx2", m.avx2, true);
    out << "  },\n";
    out << "  \"map\": {\n";
    writeStr(out, "path", m.map_path);
    writeStr(out, "id", m.map_id);
    writeStr(out, "sha256", m.map_sha256);
    out << "    \"grid_rows\": " << m.grid_rows << ",\n";
    out << "    \"grid_cols\": " << m.grid_cols << ",\n";
    out << "    \"world_w_m\": " << m.world_w_m << ",\n";
    out << "    \"world_h_m\": " << m.world_h_m << ",\n";
    out << "    \"scale_factor\": " << m.scale << ",\n";
    out << "    \"cell_size_cm\": " << m.cell_size_cm << "\n";
    out << "  },\n";
    out << "  \"htrip\": {\n";
    out << "    \"tile_size_b\": " << m.tile_size_b << ",\n";
    out << "    \"group_size_g\": " << m.group_size_g << ",\n";
    out << "    \"tree_depth\": " << m.tree_depth << ",\n";
    out << "    \"total_leaves\": " << m.total_leaves << "\n";
    out << "  },\n";
    out << "  \"fleet\": {\n";
    out << "    \"num_robots\": " << m.num_robots << ",\n";
    out << "    \"max_frontiers\": " << m.max_frontiers << ",\n";
    out << "    \"seed\": " << m.seed << ",\n";
    writeBool(out, "auto_clock_seed_on_reset", m.auto_clock_seed_on_reset, true);
    out << "  },\n";
    out << "  \"robot\": {\n";
    out << "    \"v_max\": " << m.v_max << ",\n";
    out << "    \"omega_max\": " << m.omega_max << ",\n";
    out << "    \"radius_m\": " << m.radius_m << ",\n";
    out << "    \"lidar_range_m\": " << m.lidar_range_m << ",\n";
    out << "    \"lidar_span_deg\": " << m.lidar_span_deg << ",\n";
    out << "    \"lidar_num_rays\": " << m.lidar_num_rays << "\n";
    out << "  },\n";
    out << "  \"control\": {\n";
    out << "    \"dt\": " << m.dt << ",\n";
    out << "    \"sim_step_multiplier\": " << m.sim_step_multiplier << ",\n";
    out << "    \"max_steps\": " << m.max_steps << ",\n";
    out << "    \"stop_coverage\": " << m.stop_coverage << ",\n";
    out << "    \"replan_period_ticks\": " << m.replan_period_ticks << ",\n";
    out << "    \"persistence_window_ticks\": " << m.persistence_window_ticks << ",\n";
    writeBool(out, "enable_async", m.enable_async);
    writeBool(out, "enable_auto_respawn", m.enable_auto_respawn);
    out << "    \"auto_respawn_ticks\": " << m.auto_respawn_ticks << ",\n";
    writeBool(out, "enable_live_frontier_projection", m.enable_live_frontier_projection);
    writeBool(out, "enable_fog_of_war", m.enable_fog_of_war);
    writeBool(out, "enable_bfs", m.enable_bfs);
    writeBool(out, "enable_astar", m.enable_astar);
    out << "    \"bfs_compare_interval\": " << m.bfs_compare_interval << "\n";
    out << "  },\n";
    out << "  \"capture\": {\n";
    writeStr(out, "video_quality", m.video_quality);
    out << "    \"video_fps\": " << m.video_fps << ",\n";
    out << "    \"video_resolution_px\": " << m.video_resolution_px << ",\n";
    writeBool(out, "video_delete_ppm", m.video_delete_ppm);
    out << "    \"stills\": " << m.stills << ",\n";
    writeStr(out, "video_path", m.video_path);
    out << "    \"still_paths\": [";
    for (size_t i = 0; i < m.still_paths.size(); ++i) {
        if (i > 0) out << ", ";
        out << "\"" << jsonEscape(m.still_paths[i]) << "\"";
    }
    out << "]\n";
    out << "  },\n";
    out << "  \"outcome\": {\n";
    out << "    \"total_ticks\": " << m.total_ticks << ",\n";
    out << "    \"final_coverage\": " << m.final_coverage << ",\n";
    out << "    \"wall_time_s\": " << m.wall_time_s << ",\n";
    out << "    \"bfs_compare_count\": " << m.bfs_compare_count << ",\n";
    out << "    \"paired_htrip_mean_us\": " << m.paired_htrip_mean_us << ",\n";
    out << "    \"bfs_query_mean_us\": " << m.bfs_query_mean_us << ",\n";
    out << "    \"mean_speedup_vs_bfs\": " << m.mean_speedup_vs_bfs << "\n";
    out << "  }\n";
    out << "}\n";
    return true;
}

}  // namespace htrip
