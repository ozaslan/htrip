#include "htrip/robot_simulator.hpp"
#include "htrip/arena_recorder.hpp"
#include "htrip/trial_manifest.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void printUsage(const char* argv0) {
    std::cout <<
        "Usage: " << argv0 << " [options]\n"
        "\n"
        "H-TRIP multi-robot exploration trial (headless). Writes experiment.json,\n"
        "mission telemetry, optional Medium-quality arena video, and coverage stills.\n"
        "\n"
        "Required for a mapped trial:\n"
        "  --map PATH                 MovingAI .map file\n"
        "\n"
        "Hierarchy / fleet:\n"
        "  --b N                      leaf tile size (8, 16, 32). Forbidden: 64\n"
        "  --g N                      grouping factor (2 or 4). Forbidden: g>=16\n"
        "  --robots N                 fleet size 1..32 (default 4)\n"
        "  --frontiers K              max frontier targets (default 1024, slider max)\n"
        "  --seed N                   RNG seed (required for a repeatable run)\n"
        "\n"
        "Stop conditions:\n"
        "  --steps N                  max ticks (default 4000)\n"
        "  --stop-coverage F          stop at this explored fraction (default 0.90)\n"
        "  --dt S                     sim step seconds (default 0.05)\n"
        "  --sim_step N               ticks advanced per printed step (default 1)\n"
        "\n"
        "Map / robot (GUI defaults):\n"
        "  --scale N                  1, 2 or 4 (default 1)\n"
        "  --cell_size_cm F           unscaled cell size (default 5)\n"
        "  --range F                  lidar metres (default 8)\n"
        "  --span F                   lidar FOV degrees (default 360)\n"
        "  --rays N                   lidar beams (default 150)\n"
        "  --v-max F                  linear speed m/s (default 0.60)\n"
        "  --omega-max F              turn rate rad/s (default 1.80)\n"
        "  --radius F                 chassis radius m (default 0.18)\n"
        "  --replan-period N          m-TSP cadence ticks (default 20)\n"
        "  --persistence N            goal persistence ticks (default 50)\n"
        "  --async 0|1                background H-TRIP updates (default 1)\n"
        "  --auto-respawn 0|1         stuck-robot respawn (default 1)\n"
        "  --no-fog                   reveal entire floorplan (no fog of war; full-graph BFS/A*)\n"
        "  --fog 0|1                  fog-of-war toggle (default 1)\n"
        "  --bfs 0|1                  concurrent k x BFS alongside H-TRIP (default 0)\n"
        "  --astar 0|1                paired A* baseline comparison (default 0, slow)\n"
        "\n"
        "Capture:\n"
        "  --record                   write arena video\n"
        "  --keep-ppm                 preserve all raw uncompressed PPM frames (no cleanup)\n"
        "  --video-quality low|medium|high   (default medium)\n"
        "  --video-fps N              (default 20)\n"
        "  --stills N                 coverage stills (default 5: start/25/50/75/end)\n"
        "  --log-dir PATH             trial directory (telemetry + video + experiment.json)\n"
        "  --output CSV               optional extra per-step CSV\n"
        "  --help\n";
}

std::string mapIdFromPath(const std::string& path) {
    std::filesystem::path p(path);
    return p.stem().string();
}

bool parseBoolFlag(const std::string& v) {
    return !(v == "0" || v == "false" || v == "False" || v == "no");
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(true);
    std::cout << std::unitbuf;
    std::cerr << std::unitbuf;
    int scale = 1;
    double cell_size_cm = 5.0;
    int num_robots = 4;
    double range_m = 8.0;
    double span_deg = 360.0;
    int num_rays = 150;
    int max_steps = 4000;
    int sim_step_multiplier = 1;
    int max_frontiers = 1024;
    int tile_b = 16;
    int group_g = 2;
    uint32_t seed = 0;
    double dt = 0.05;
    double stop_coverage = 0.90;
    double v_max = 0.60;
    double omega_max = 1.80;
    double radius_m = 0.18;
    int replan_period = 20;
    int persistence = 50;
    bool enable_async = true;
    bool enable_auto_respawn = true;
    int auto_respawn_ticks = 100;
    bool enable_fog = true;
    bool enable_bfs = false;
    bool enable_astar = false;
    bool do_record = false;
    bool keep_ppm = false;
    std::string video_quality_name = "medium";
    int video_fps = 20;
    int stills = 5;
    std::string csv_out;
    std::string log_dir;
    std::string map_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << name << " requires a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--scale") scale = std::stoi(need("--scale"));
        else if (arg == "--cell_size_cm") cell_size_cm = std::stod(need("--cell_size_cm"));
        else if (arg == "--robots") num_robots = std::stoi(need("--robots"));
        else if (arg == "--frontiers") max_frontiers = std::clamp(std::stoi(need("--frontiers")), 8, 1024);
        else if (arg == "--range") range_m = std::stod(need("--range"));
        else if (arg == "--span") span_deg = std::stod(need("--span"));
        else if (arg == "--rays") num_rays = std::stoi(need("--rays"));
        else if (arg == "--steps") max_steps = std::stoi(need("--steps"));
        else if (arg == "--sim_step") sim_step_multiplier = std::clamp(std::stoi(need("--sim_step")), 1, 8);
        else if (arg == "--output") csv_out = need("--output");
        else if (arg == "--log-dir") log_dir = need("--log-dir");
        else if (arg == "--map") map_path = need("--map");
        else if (arg == "--b" || arg == "--tile") tile_b = std::stoi(need(arg.c_str()));
        else if (arg == "--g" || arg == "--group") group_g = std::stoi(need(arg.c_str()));
        else if (arg == "--seed") seed = static_cast<uint32_t>(std::stoul(need("--seed")));
        else if (arg == "--dt") dt = std::stod(need("--dt"));
        else if (arg == "--stop-coverage") stop_coverage = std::stod(need("--stop-coverage"));
        else if (arg == "--v-max") v_max = std::stod(need("--v-max"));
        else if (arg == "--omega-max") omega_max = std::stod(need("--omega-max"));
        else if (arg == "--radius") radius_m = std::stod(need("--radius"));
        else if (arg == "--replan-period") replan_period = std::stoi(need("--replan-period"));
        else if (arg == "--persistence") persistence = std::stoi(need("--persistence"));
        else if (arg == "--async") enable_async = parseBoolFlag(need("--async"));
        else if (arg == "--auto-respawn") enable_auto_respawn = parseBoolFlag(need("--auto-respawn"));
        else if (arg == "--respawn-ticks") auto_respawn_ticks = std::stoi(need("--respawn-ticks"));
        else if (arg == "--no-fog") enable_fog = false;
        else if (arg == "--fog") enable_fog = parseBoolFlag(need("--fog"));
        else if (arg == "--bfs") enable_bfs = parseBoolFlag(need("--bfs"));
        else if (arg == "--astar") enable_astar = parseBoolFlag(need("--astar"));
        else if (arg == "--record") do_record = true;
        else if (arg == "--keep-ppm" || arg == "--dump-ppm") keep_ppm = true;
        else if (arg == "--video-quality") video_quality_name = need("--video-quality");
        else if (arg == "--video-fps") video_fps = std::clamp(std::stoi(need("--video-fps")), 10, 60);
        else if (arg == "--stills") stills = std::max(0, std::stoi(need("--stills")));
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return 2;
        }
    }

    if (group_g >= 16) {
        std::cerr << "Error: g>=16 is forbidden (too slow for this campaign).\n";
        return 2;
    }
    if (tile_b != 8 && tile_b != 16 && tile_b != 32 && tile_b != 64) {
        std::cerr << "Error: --b must be 8, 16, 32, or 64.\n";
        return 2;
    }
    if (group_g != 2 && group_g != 4) {
        std::cerr << "Error: --g must be 2 or 4.\n";
        return 2;
    }
    num_robots = std::clamp(num_robots, 1, 32);
    if (seed == 0) {
        seed = 2026091301u;
        std::cerr << "Warning: --seed not set; using default " << seed << "\n";
    }
    if (log_dir.empty()) {
        log_dir = "campaigns/closed_loop/adhoc";
    }

    std::error_code ec;
    std::filesystem::create_directories(log_dir, ec);

    htrip::TrialManifest manifest;
    manifest.command_line = htrip::joinCommandLine(argc, argv);
    manifest.map_path = map_path;
    manifest.map_id = mapIdFromPath(map_path);
    manifest.map_sha256 = htrip::sha256File(map_path);
    manifest.scale = scale;
    manifest.cell_size_cm = cell_size_cm;
    manifest.tile_size_b = tile_b;
    manifest.group_size_g = group_g;
    manifest.num_robots = num_robots;
    manifest.max_frontiers = max_frontiers;
    manifest.seed = seed;
    manifest.max_steps = max_steps;
    manifest.stop_coverage = stop_coverage;
    manifest.dt = dt;
    manifest.sim_step_multiplier = sim_step_multiplier;
    manifest.lidar_range_m = range_m;
    manifest.lidar_span_deg = span_deg;
    manifest.lidar_num_rays = num_rays;
    manifest.v_max = v_max;
    manifest.omega_max = omega_max;
    manifest.radius_m = radius_m;
    manifest.replan_period_ticks = replan_period;
    manifest.persistence_window_ticks = persistence;
    manifest.enable_async = enable_async;
    manifest.enable_auto_respawn = enable_auto_respawn;
    manifest.auto_respawn_ticks = auto_respawn_ticks;
    manifest.enable_live_frontier_projection = true;
    manifest.enable_fog_of_war = enable_fog;
    manifest.enable_bfs = enable_bfs;
    manifest.enable_astar = enable_astar;
    manifest.bfs_compare_interval = enable_bfs ? 1 : 0;
    manifest.auto_clock_seed_on_reset = false;
    manifest.video_quality = video_quality_name;
    manifest.video_fps = video_fps;
    manifest.stills = stills;
    manifest.log_dir = log_dir;
    htrip::fillRuntimeEnvironment(manifest);
    manifest.status = "running";
    htrip::writeTrialManifest((std::filesystem::path(log_dir) / "experiment.json").string(), manifest);

    std::cout << "========================================================================================\n";
    std::cout << "  H-TRIP MULTI-ROBOT CONTINUOUS SE(2) EXPLORATION TRIAL\n";
    std::cout << "========================================================================================\n";
    if (!map_path.empty()) {
        std::cout << "  Map:                   " << map_path << "\n";
        std::cout << "  Map SHA256:            " << manifest.map_sha256 << "\n";
    }
    std::cout << "  Resolution Multiplier: " << scale << "x\n";
    std::cout << "  Cell Size:             " << cell_size_cm << " cm\n";
    std::cout << "  H-TRIP b,g:            " << tile_b << ", " << group_g << "\n";
    std::cout << "  Active Robots:         " << num_robots << "\n";
    std::cout << "  Frontier Targets (K):  " << max_frontiers << "\n";
    std::cout << "  Seed:                  " << seed << "\n";
    std::cout << "  LiDAR:                 " << range_m << " m / " << span_deg << " deg / " << num_rays << " rays\n";
    std::cout << "  Record:                " << (do_record ? (keep_ppm ? "YES (all PPM frames preserved)" : "YES (MP4 stream)") : "NO") << "\n";
    std::cout << "  BFS / A* Baselines:    " << (enable_bfs ? "ENABLED (paired speedup on every tick)" : "DISABLED") << "\n";
    std::cout << "  Log dir:               " << log_dir << "\n";
    std::cout << "----------------------------------------------------------------------------------------\n\n";

    htrip::GridMap base_map(64, 64, htrip::CellState::Free);
    if (!map_path.empty()) {
        if (!base_map.loadMovingAI(map_path)) {
            std::cerr << "Error: Failed to load map from " << map_path << "\n";
            manifest.status = "failed";
            manifest.stop_reason = "map_load_failed";
            htrip::writeTrialManifest((std::filesystem::path(log_dir) / "experiment.json").string(), manifest);
            return 1;
        }
    } else {
        int R0 = 64;
        int C0 = 64;
        for (int r = 0; r < R0; ++r) {
            base_map.setObserved(r, 0, htrip::CellState::Obstacle);
            base_map.setObserved(r, C0 - 1, htrip::CellState::Obstacle);
        }
        for (int c = 0; c < C0; ++c) {
            base_map.setObserved(0, c, htrip::CellState::Obstacle);
            base_map.setObserved(R0 - 1, c, htrip::CellState::Obstacle);
        }
        for (int r = 0; r < R0; ++r) {
            if (r != 16 && r != 48) base_map.setObserved(r, 32, htrip::CellState::Obstacle);
        }
        for (int c = 0; c < C0; ++c) {
            if (c != 16 && c != 48) base_map.setObserved(32, c, htrip::CellState::Obstacle);
        }
    }

    htrip::RobotConfig rcfg;
    rcfg.radius_m = radius_m;
    rcfg.v_max = v_max;
    rcfg.omega_max = omega_max;
    rcfg.lidar.range_m = range_m;
    rcfg.lidar.span_deg = span_deg;
    rcfg.lidar.num_rays = num_rays;

    htrip::MultiRobotSimulator sim;
    sim.enable_fog_of_war = enable_fog;
    sim.enable_bfs_baseline = enable_bfs;
    sim.enable_astar_baseline = enable_astar;
    sim.bfs_compare_interval = enable_bfs ? 1 : 0;
    if (enable_bfs) {
        sim.baseline_mode = htrip::BaselineExecutionMode::Synchronous;
    }
    sim.max_frontiers_eval = max_frontiers;
    sim.init(base_map, cell_size_cm / 100.0, scale, num_robots, rcfg, tile_b, group_g, seed, true, enable_fog);
    sim.replan_period_ticks = replan_period;
    sim.persistence_window_ticks = persistence;
    sim.enable_async_updates = enable_async;
    sim.enable_auto_respawn = enable_auto_respawn;
    sim.auto_respawn_ticks = auto_respawn_ticks;
    sim.enable_live_frontier_projection = true;

    std::cout << "  Initialized World: " << sim.ground_truth.rows << "x" << sim.ground_truth.cols
              << " cells (" << std::fixed << std::setprecision(2) << sim.world_w_m << " m x "
              << sim.world_h_m << " m physical space)\n";
    std::cout << "  Total Free Floor:  " << sim.total_free_cells_gt << " cells\n";
    std::cout << "  Map State:         " << (enable_fog ? "Fog of War (Exploration)" : "Fully Revealed (Patrol/Inspection)") << "\n\n";

    std::cout << std::left
              << std::setw(8)  << "Tick"
              << std::setw(12) << (enable_fog ? "Explored" : "Targets")
              << std::setw(12) << (enable_fog ? "Frontiers" : "Active VPs")
              << std::setw(16) << "H-TRIP Upd(us)"
              << std::setw(16) << "H-TRIP Qry(us)"
              << "\n";
    std::cout << std::string(64, '-') << "\n";

    std::ofstream csv;
    if (!csv_out.empty()) {
        csv.open(csv_out);
        csv << "tick,explored_fraction,frontiers_k,htrip_update_us,htrip_query_us\n";
    }

    sim.startMissionLogging(log_dir, "mission", false);

    htrip::ArenaRecorder recorder;
    recorder.setQuality(htrip::parseVideoQuality(video_quality_name));
    recorder.fps = video_fps;
    recorder.delete_ppm_after_compile = !keep_ppm;
    manifest.video_quality = recorder.qualityTag();
    manifest.video_resolution_px = recorder.resolution_px;
    if (do_record) {
        recorder.stream_mp4 = !keep_ppm;
        recorder.startRecordingAt(log_dir);
    }

    const std::filesystem::path still_dir = std::filesystem::path(log_dir) / "stills";
    std::filesystem::create_directories(still_dir, ec);
    auto save_named_still = [&](const std::string& stem) {
        if (stills <= 0) return;
        const std::string path = (still_dir / (stem + ".ppm")).string();
        if (recorder.saveStill(sim, path)) {
            manifest.still_paths.push_back(path);
        }
    };

    bool still_25 = false;
    bool still_50 = false;
    bool still_75 = false;
    save_named_still("still_000_start");

    const auto wall0 = std::chrono::steady_clock::now();
    std::string stop_reason = "max_ticks";

    for (int step = 1; step <= max_steps; ++step) {
        for (int s = 0; s < sim_step_multiplier; ++s) {
            sim.step(dt);
        }

        if (do_record) {
            recorder.captureTick(sim);
        }

        const double cov = sim.explorationFraction();
        if (stills >= 5) {
            if (!sim.enable_fog_of_war) {
                if (!still_25 && step >= max_steps * 0.25) {
                    save_named_still("still_025_progress");
                    still_25 = true;
                }
                if (!still_50 && step >= max_steps * 0.50) {
                    save_named_still("still_050_progress");
                    still_50 = true;
                }
                if (!still_75 && step >= max_steps * 0.75) {
                    save_named_still("still_075_progress");
                    still_75 = true;
                }
            } else {
                if (!still_25 && cov >= 0.25) {
                    save_named_still("still_025_coverage");
                    still_25 = true;
                }
                if (!still_50 && cov >= 0.50) {
                    save_named_still("still_050_coverage");
                    still_50 = true;
                }
                if (!still_75 && cov >= 0.75) {
                    save_named_still("still_075_coverage");
                    still_75 = true;
                }
            }
        }

        const bool print_row = (step == 1) || (step % 20 == 0) || (sim.enable_fog_of_war && cov >= stop_coverage);
        if (print_row) {
            std::string cov_col = sim.enable_fog_of_war ? (std::to_string(static_cast<int>(cov * 100.0)) + "%")
                                                        : ("Vis: " + std::to_string(sim.targets_visited_count));
            std::cout << std::left
                      << std::setw(8)  << sim.total_ticks
                      << std::setw(12) << cov_col
                      << std::setw(12) << sim.frontiers.size()
                      << std::setw(16) << std::fixed << std::setprecision(1) << sim.last_htrip_update_us.load()
                      << std::setw(16) << std::fixed << std::setprecision(1) << sim.last_htrip_query_us.load()
                      << "\n";
        }

        if (csv.is_open()) {
            csv << sim.total_ticks << ","
                << cov << ","
                << sim.frontiers.size() << ","
                << sim.last_htrip_update_us.load() << ","
                << sim.last_htrip_query_us.load() << "\n";
        }

        if (sim.enable_fog_of_war && cov >= stop_coverage) {
            stop_reason = "coverage";
            break;
        }
    }

    save_named_still("still_end");

    if (do_record) {
        recorder.stopRecording();
        if (!recorder.compileVideoBlocking()) {
            std::cerr << "Warning: video compile failed: " << recorder.status_msg << "\n";
        }
        manifest.video_path = recorder.last_video_path;
    }

    if (csv.is_open()) {
        csv.close();
    }

    if (sim.isMissionLogging()) {
        if (sim.enable_bfs_baseline) {
            std::cout << "Waiting for in-flight BFS so paired H-TRIP comparisons can be written.\n";
        }
        sim.stopMissionLogging();
    }

    const auto wall1 = std::chrono::steady_clock::now();
    htrip::fillFromSimulator(manifest, sim);
    manifest.stop_reason = stop_reason;
    manifest.status = "complete";
    manifest.wall_time_s = std::chrono::duration<double>(wall1 - wall0).count();
    manifest.seed = seed;
    htrip::writeTrialManifest((std::filesystem::path(log_dir) / "experiment.json").string(), manifest);

    std::cout << std::string(76, '=') << "\n";
    std::cout << "  Stop reason:           " << stop_reason << "\n";
    if (!sim.enable_fog_of_war) {
        std::cout << "  Targets Visited:       " << sim.targets_visited_count << "\n";
        std::cout << "  Map State:             Fully Revealed (No Fog of War)\n";
    } else {
        std::cout << "  Final coverage:        " << std::fixed << std::setprecision(2)
                  << (sim.explorationFraction() * 100.0) << "%\n";
    }
    std::cout << "  Ticks:                 " << sim.total_ticks << "\n";
    std::cout << "  BFS compares:          " << manifest.bfs_compare_count << "\n";
    if (manifest.bfs_compare_count > 0) {
        std::cout << "  Paired H-TRIP mean:    " << manifest.paired_htrip_mean_us << " us\n";
        std::cout << "  Paired BFS mean:       " << manifest.bfs_query_mean_us << " us\n";
        std::cout << "  Mean speedup vs BFS:   " << manifest.mean_speedup_vs_bfs << "\n";
    }
    std::cout << "  Wall time:             " << manifest.wall_time_s << " s\n";
    std::cout << "  experiment.json:       " << (std::filesystem::path(log_dir) / "experiment.json") << "\n";
    if (!manifest.video_path.empty()) {
        std::cout << "  Video:                 " << manifest.video_path << "\n";
    }
    std::cout << "========================================================================================\n";
    return 0;
}
