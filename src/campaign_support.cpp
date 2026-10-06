#include "htrip/campaign_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <sys/resource.h>
#include <system_error>

namespace htrip {
namespace campaign {

namespace {

constexpr int kDr[4] = {-1, 1, 0, 0};
constexpr int kDc[4] = {0, 0, -1, 1};

bool gtPassable(const GridMap& grid, int r, int c) {
    if (!grid.inBounds(r, c)) return false;
    return grid.ground_truth[static_cast<size_t>(grid.index(r, c))] == CellState::Free;
}

} // namespace

std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size() + 8);
    for (char ch : s) {
        switch (ch) {
        case '\\': o += "\\\\"; break;
        case '"': o += "\\\""; break;
        case '\n': o += "\\n"; break;
        case '\r': o += "\\r"; break;
        case '\t': o += "\\t"; break;
        default: o += ch; break;
        }
    }
    return o;
}

std::string shellCapture(const std::string& cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) return {};
    std::string out;
    char buf[512];
    while (fgets(buf, sizeof(buf), pipe) != nullptr) out += buf;
    pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out;
}

std::string findRepoRoot() {
    std::filesystem::path p = std::filesystem::current_path();
    for (int i = 0; i < 8; ++i) {
        if (std::filesystem::exists(p / ".git") && std::filesystem::exists(p / "include" / "htrip")) {
            return p.string();
        }
        if (!p.has_parent_path() || p == p.parent_path()) break;
        p = p.parent_path();
    }
    return std::filesystem::current_path().string();
}

std::string resolveExistingPath(const std::string& path) {
    if (std::filesystem::exists(path)) return path;
    const std::string root = findRepoRoot();
    std::filesystem::path rel = std::filesystem::path(root) / path;
    if (std::filesystem::exists(rel)) return rel.string();
    std::filesystem::path up = std::filesystem::current_path() / ".." / ".." / path;
    if (std::filesystem::exists(up)) return std::filesystem::weakly_canonical(up).string();
    return path;
}

std::string fileSha256(const std::string& path) {
    if (path.empty() || !std::filesystem::exists(path)) return "unknown";
    std::string out = shellCapture("sha256sum -b \"" + path + "\" 2>/dev/null | awk '{print $1}'");
    return out.empty() ? "unknown" : out;
}

std::string gitSha(const std::string& repo) {
    std::string out = shellCapture("git -C \"" + repo + "\" rev-parse HEAD 2>/dev/null");
    return out.empty() ? "unknown" : out;
}

bool gitDirty(const std::string& repo) {
    std::string out = shellCapture("git -C \"" + repo + "\" status --porcelain 2>/dev/null");
    return !out.empty();
}

std::int64_t peakRssBytes() {
    rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
#if defined(__linux__)
    return static_cast<std::int64_t>(ru.ru_maxrss) * 1024;
#else
    return static_cast<std::int64_t>(ru.ru_maxrss);
#endif
}

double percentileSorted(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) return 0.0;
    const double idx = (p / 100.0) * static_cast<double>(sorted.size() - 1);
    const size_t i0 = static_cast<size_t>(std::floor(idx));
    const size_t i1 = std::min(i0 + 1, sorted.size() - 1);
    const double frac = idx - static_cast<double>(i0);
    return sorted[i0] + frac * (sorted[i1] - sorted[i0]);
}

TimingStats summarizeMs(std::vector<double> samples) {
    TimingStats s;
    if (samples.empty()) return s;
    std::sort(samples.begin(), samples.end());
    s.min_ms = samples.front();
    s.max_ms = samples.back();
    s.mean_ms = std::accumulate(samples.begin(), samples.end(), 0.0) / static_cast<double>(samples.size());
    s.median_ms = percentileSorted(samples, 50.0);
    s.q1_ms = percentileSorted(samples, 25.0);
    s.q3_ms = percentileSorted(samples, 75.0);
    return s;
}

std::vector<int> parseIntList(const std::string& csv, const std::vector<int>& fallback) {
    if (csv.empty()) return fallback;
    std::vector<int> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        out.push_back(std::atoi(tok.c_str()));
    }
    return out.empty() ? fallback : out;
}

std::string joinInts(const std::vector<int>& v) {
    std::string o;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) o += ',';
        o += std::to_string(v[i]);
    }
    return o;
}

bool loadFullyRevealed(const std::string& path, GridMap& grid) {
    const std::string resolved = resolveExistingPath(path);
    if (!grid.loadMovingAI(resolved)) return false;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            grid.setObserved(r, c, grid.ground_truth[static_cast<size_t>(grid.index(r, c))]);
        }
    }
    return true;
}

size_t countGtPassable(const GridMap& grid) {
    size_t n = 0;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            if (gtPassable(grid, r, c)) ++n;
        }
    }
    return n;
}

ConnectedComponent largestGtComponent(const GridMap& grid) {
    ConnectedComponent cc;
    const int R = grid.rows;
    const int C = grid.cols;
    cc.cell_comp.assign(static_cast<size_t>(R * C), -1);
    int current = 0;
    std::vector<size_t> sizes;
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            const int idx = grid.index(r, c);
            if (!gtPassable(grid, r, c) || cc.cell_comp[static_cast<size_t>(idx)] != -1) continue;
            size_t sz = 0;
            std::queue<Point> q;
            q.push({r, c});
            cc.cell_comp[static_cast<size_t>(idx)] = current;
            while (!q.empty()) {
                const Point p = q.front();
                q.pop();
                ++sz;
                for (int k = 0; k < 4; ++k) {
                    const int nr = p.r + kDr[k];
                    const int nc = p.c + kDc[k];
                    if (!gtPassable(grid, nr, nc)) continue;
                    const int nidx = grid.index(nr, nc);
                    if (cc.cell_comp[static_cast<size_t>(nidx)] == -1) {
                        cc.cell_comp[static_cast<size_t>(nidx)] = current;
                        q.push({nr, nc});
                    }
                }
            }
            sizes.push_back(sz);
            ++current;
        }
    }
    size_t max_sz = 0;
    int best = -1;
    for (int i = 0; i < current; ++i) {
        if (sizes[static_cast<size_t>(i)] > max_sz) {
            max_sz = sizes[static_cast<size_t>(i)];
            best = i;
        }
    }
    cc.id = best;
    cc.size = max_sz;
    return cc;
}

void markObservedFreeExcept(GridMap& grid, const std::vector<Point>& withheld) {
    std::vector<uint8_t> hide(static_cast<size_t>(grid.rows * grid.cols), 0);
    for (Point p : withheld) {
        if (grid.inBounds(p)) hide[static_cast<size_t>(grid.index(p))] = 1;
    }
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const size_t idx = static_cast<size_t>(grid.index(r, c));
            if (grid.ground_truth[idx] == CellState::Free && hide[idx] == 0) {
                grid.setObserved(r, c, CellState::Free);
            } else if (grid.ground_truth[idx] == CellState::Obstacle) {
                grid.setObserved(r, c, CellState::Obstacle);
            } else {
                grid.setObserved(r, c, CellState::Unknown);
            }
        }
    }
}

std::vector<Point> sampleSector(const GridMap& grid, const ConnectedComponent& cc, int divisions,
                                int origin_r, int origin_c) {
    const int R = grid.rows;
    const int C = grid.cols;
    std::vector<Point> sampled;
    const int step_r = std::max(1, R / divisions);
    const int step_c = std::max(1, C / divisions);
    const int half_r = std::max(1, step_r / 2);
    const int half_c = std::max(1, step_c / 2);
    for (int sr = 0; sr < divisions; ++sr) {
        for (int sc = 0; sc < divisions; ++sc) {
            const int center_r = sr * step_r + step_r / 2 + origin_r;
            const int center_c = sc * step_c + step_c / 2 + origin_c;
            Point best_p(-1, -1);
            int min_dist = 999999;
            for (int dr = -half_r; dr <= half_r; ++dr) {
                for (int dc = -half_c; dc <= half_c; ++dc) {
                    const int r = center_r + dr;
                    const int c = center_c + dc;
                    if (r < 0 || r >= R || c < 0 || c >= C) continue;
                    const int idx = grid.index(r, c);
                    if (gtPassable(grid, r, c) && cc.cell_comp[static_cast<size_t>(idx)] == cc.id) {
                        const int d = std::abs(dr) + std::abs(dc);
                        if (d < min_dist) {
                            min_dist = d;
                            best_p = {r, c};
                        }
                    }
                }
            }
            if (best_p.r != -1) sampled.push_back(best_p);
        }
    }
    return sampled;
}

std::vector<Point> sampleSpread(const GridMap& grid, const ConnectedComponent& cc, std::uint64_t seed, size_t n) {
    std::vector<Point> cells;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const int idx = grid.index(r, c);
            if (cc.cell_comp[static_cast<size_t>(idx)] == cc.id) cells.push_back({r, c});
        }
    }
    std::vector<Point> out;
    if (cells.empty() || n == 0) return out;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> dist(0, cells.size() - 1);
    out.push_back(cells[dist(rng)]);
    while (out.size() < n && out.size() < cells.size()) {
        int best_d = -1;
        Point best = out.back();
        for (const Point& cand : cells) {
            int min_d = 1 << 30;
            for (const Point& chosen : out) {
                min_d = std::min(min_d, cand.manhattan(chosen));
            }
            if (min_d > best_d) {
                best_d = min_d;
                best = cand;
            }
        }
        bool dup = false;
        for (const Point& p : out) {
            if (p == best) {
                dup = true;
                break;
            }
        }
        if (dup) break;
        out.push_back(best);
    }
    return out;
}

std::vector<Point> sampleCluster(const GridMap& grid, const ConnectedComponent& cc, std::uint64_t seed, size_t n) {
    std::vector<Point> cells;
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const int idx = grid.index(r, c);
            if (cc.cell_comp[static_cast<size_t>(idx)] == cc.id) cells.push_back({r, c});
        }
    }
    if (cells.empty() || n == 0) return {};
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> dist(0, cells.size() - 1);
    const Point c0 = cells[dist(rng)];
    std::sort(cells.begin(), cells.end(), [&](const Point& a, const Point& b) {
        const int ca = std::max(std::abs(a.r - c0.r), std::abs(a.c - c0.c));
        const int cb = std::max(std::abs(b.r - c0.r), std::abs(b.c - c0.c));
        if (ca != cb) return ca < cb;
        return a.manhattan(c0) < b.manhattan(c0);
    });
    if (cells.size() > n) cells.resize(n);
    return cells;
}

void padFarthest(const GridMap& grid, const ConnectedComponent& cc, std::vector<Point>& out, size_t n) {
    (void)grid;
    if (out.size() >= n) return;
    std::vector<Point> cells;
    cells.reserve(cc.size);
    for (int r = 0; r < grid.rows; ++r) {
        for (int c = 0; c < grid.cols; ++c) {
            const int idx = grid.index(r, c);
            if (cc.cell_comp[static_cast<size_t>(idx)] == cc.id) cells.push_back({r, c});
        }
    }
    while (out.size() < n && out.size() < cells.size()) {
        int best_d = -1;
        Point best = out.empty() ? cells.front() : out.back();
        for (const Point& cand : cells) {
            int min_d = 1 << 30;
            for (const Point& chosen : out) {
                min_d = std::min(min_d, cand.manhattan(chosen));
            }
            if (min_d > best_d) {
                best_d = min_d;
                best = cand;
            }
        }
        bool dup = false;
        for (const Point& p : out) {
            if (p == best) {
                dup = true;
                break;
            }
        }
        if (dup) break;
        out.push_back(best);
    }
}

std::vector<Point> sampleGeometry(const std::string& name, const GridMap& grid, const ConnectedComponent& cc,
                                  std::uint64_t seed, size_t n, int origin_r, int origin_c) {
    if (name == "spread") return sampleSpread(grid, cc, seed, n);
    if (name == "cluster") return sampleCluster(grid, cc, seed, n);
    int divisions = 8;
    while (static_cast<size_t>(divisions) * static_cast<size_t>(divisions) < n && divisions < 64) {
        divisions *= 2;
    }
    auto pts = sampleSector(grid, cc, divisions, origin_r, origin_c);
    std::vector<Point> uniq;
    uniq.reserve(pts.size());
    for (const Point& p : pts) {
        if (std::find(uniq.begin(), uniq.end(), p) == uniq.end()) uniq.push_back(p);
    }
    pts = std::move(uniq);
    if (pts.size() < n) padFarthest(grid, cc, pts, n);
    if (pts.size() > n) pts.resize(n);
    return pts;
}

MatrixCheck compareMatrices(const std::vector<dist_t>& a, const std::vector<dist_t>& b) {
    MatrixCheck chk;
    const size_t n = std::min(a.size(), b.size());
    chk.entries = std::max(a.size(), b.size());
    if (a.size() != b.size()) chk.mismatches = chk.entries;
    for (size_t i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            ++chk.mismatches;
            chk.max_abs_diff = std::max(chk.max_abs_diff, std::abs(static_cast<int>(a[i]) - static_cast<int>(b[i])));
        }
        if (a[i] < INF) chk.max_finite = std::max(chk.max_finite, a[i]);
        if (b[i] < INF) chk.max_finite = std::max(chk.max_finite, b[i]);
    }
    return chk;
}

bool csvNeedsHeader(const std::string& path) {
    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) return true;
    return std::filesystem::file_size(path, ec) == 0;
}

void csvAppendLine(const std::string& path, const std::string& header, const std::string& line) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
    const bool header_needed = csvNeedsHeader(path);
    std::ofstream out(path, std::ios::app);
    if (header_needed) out << header << '\n';
    out << line << '\n';
}

std::string exactnessHeader() {
    return "suite_id,family,maps_or_config,b,g,seed,n_openings,n_blockages,n_entries,reference,"
           "mismatches,max_abs_diff,max_finite_dist,binary_sha,git_sha";
}

std::string queryVsKHeader() {
    return "map_id,b,g,geometry,set_id,seed,k,method,warmup,repeats,"
           "median_ms,mean_ms,min_ms,max_ms,n_pairs,index_bytes,peak_rss_bytes,"
           "max_finite_dist,mismatches,binary_sha,regime";
}

std::string queryVsKRepeatsHeader() {
    return "map_id,b,g,geometry,set_id,seed,k,method,repeat_idx,ms,regime";
}

std::string updateHeader() {
    return "map_id,b,g,pattern,m,batch_id,seed,n_touched_leaves,n_dirty_ancestors,n_pivots,"
           "htrip_upd_median_ms,htrip_upd_min_ms,htrip_upd_max_ms,rebuild_median_ms,kbfs_dt_median_ms,"
           "k_check,entries_checked,mismatches,peak_rss_bytes,binary_sha,regime";
}

std::string rssHeader() {
    return "map_id,b,g,regime,replica_count,k,index_bytes,peak_rss_bytes,notes,binary_sha";
}

std::string fullBuildHeader() {
    return "map_id,b,g,set_id,seed,build_ms,index_bytes,peak_rss_bytes,binary_sha,regime";
}

} // namespace campaign
} // namespace htrip
