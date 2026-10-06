#include "htrip/frontier.hpp"
#include <queue>
#include <algorithm>
#include <limits>

namespace htrip {

// ============================================================================
// Frontier Cell Classification with Clearance & Depth Filtering
// ============================================================================

bool FrontierManager::isFrontierCell(const GridMap& grid, Point p) noexcept {
    // A frontier cell must be in currently observed traversable free space
    if (!grid.isPassable(p)) return false;

    // Stage 1: Check 4-connected cardinal neighbors for at least one unknown cell
    static const int dr4[4] = {-1, 1, 0, 0};
    static const int dc4[4] = {0, 0, -1, 1};

    int num_unknown_4 = 0;
    for (int i = 0; i < 4; ++i) {
        const int nr = p.r + dr4[i];
        const int nc = p.c + dc4[i];
        if (grid.inBounds(nr, nc) && grid.getObserved(nr, nc) == CellState::Unknown) {
            num_unknown_4++;
        }
    }
    if (num_unknown_4 == 0) return false;

    // Stage 2: 5x5 neighborhood composition analysis around candidate cell
    int num_unknown_5x5 = 0;
    int num_obstacle_5x5 = 0;
    for (int dr = -2; dr <= 2; ++dr) {
        for (int dc = -2; dc <= 2; ++dc) {
            if (dr == 0 && dc == 0) continue;
            const int nr = p.r + dr;
            const int nc = p.c + dc;
            if (!grid.inBounds(nr, nc)) {
                num_obstacle_5x5++;
                continue;
            }
            const CellState cs = grid.getObserved(nr, nc);
            if (cs == CellState::Unknown) num_unknown_5x5++;
            else if (cs == CellState::Obstacle) num_obstacle_5x5++;
        }
    }

    // A real explorable frontier must open up into a meaningful void of unknown space.
    // If the 5x5 patch contains fewer than 4 unknown cells, this is an isolated wall
    // shadow or corner crack already surrounded by explored ground.
    if (num_unknown_5x5 < 4) {
        return false;
    }

    // Stage 3: Clearance-aware wall proximity suppression
    // If cell is hugging an obstacle wall (clearance <= 1) and obstacles occupy significant space,
    // require substantial unknown space (>= 6 cells) to prevent wall-skirting pseudo-frontiers.
    const uint8_t clr = grid.getClearance(p);
    if (clr <= 1 && (num_obstacle_5x5 >= 6 || num_unknown_5x5 < 6)) {
        return false;
    }

    // Stage 4: Unknown region depth verification
    // Require that the adjacent unknown area has depth (at least 2 unknown cells in the local patch)
    // rather than being a 1-pixel gap/corner inside a wall.
    int valid_unknown_regions = 0;
    for (int i = 0; i < 4; ++i) {
        const int nr = p.r + dr4[i];
        const int nc = p.c + dc4[i];
        if (!grid.inBounds(nr, nc) || grid.getObserved(nr, nc) != CellState::Unknown) continue;

        int unk_count = 0;
        for (int j = 0; j < 4; ++j) {
            const int nnr = nr + dr4[j];
            const int nnc = nc + dc4[j];
            if (grid.inBounds(nnr, nnc) && grid.getObserved(nnr, nnc) == CellState::Unknown) {
                unk_count++;
            }
        }
        if (unk_count >= 2) {
            valid_unknown_regions++;
        }
    }

    return valid_unknown_regions > 0;
}

// ============================================================================
// Clearance-Aware Medoid Viewpoint Selection
// ============================================================================

Point FrontierManager::computeMedoid(const GridMap& grid, const std::vector<Point>& cluster_cells) noexcept {
    if (cluster_cells.empty()) return {0, 0};

    // Helper: searches a 7x7 neighborhood to locate a safe stand-off viewpoint with clearance >= 2
    auto findSafeViewpoint = [&](Point base) noexcept -> Point {
        if (grid.getClearance(base) >= 2) return base;
        Point best_safe = base;
        uint8_t best_clr = grid.getClearance(base);
        double min_dist = 1e9;

        for (int dr = -3; dr <= 3; ++dr) {
            for (int dc = -3; dc <= 3; ++dc) {
                const Point cand{base.r + dr, base.c + dc};
                if (!grid.inBounds(cand.r, cand.c) || !grid.isPassable(cand)) continue;
                const uint8_t clr = grid.getClearance(cand);
                const double d = std::hypot(dr, dc);
                if (clr > best_clr || (clr == best_clr && d < min_dist)) {
                    best_clr = clr;
                    min_dist = d;
                    best_safe = cand;
                }
            }
        }
        return best_safe;
    };

    if (cluster_cells.size() == 1) {
        return findSafeViewpoint(cluster_cells[0]);
    }

    // Pick best medoid heavily penalizing wall proximity
    Point best_p = cluster_cells[0];
    int64_t min_cost = std::numeric_limits<int64_t>::max();

    // Small-to-moderate clusters: evaluate exact L1 Manhattan medoid penalized by wall proximity
    if (cluster_cells.size() <= 128) {
        for (const auto& candidate : cluster_cells) {
            if (!grid.isPassable(candidate)) continue;
            int64_t dist_sum = 0;
            for (const auto& other : cluster_cells) {
                dist_sum += candidate.manhattan(other);
            }
            const uint8_t clr = grid.getClearance(candidate);
            const int64_t clr_penalty = (clr <= 1) ? 1000 : ((clr == 2) ? 100 : 0);
            const int64_t cost = dist_sum + clr_penalty;
            if (cost < min_cost) {
                min_cost = cost;
                best_p = candidate;
            }
        }
        return findSafeViewpoint(best_p);
    }

    // Large clusters: pick cell closest to Euclidean centroid with clearance bonus
    double avg_r = 0.0, avg_c = 0.0;
    for (const auto& pt : cluster_cells) {
        avg_r += pt.r;
        avg_c += pt.c;
    }
    avg_r /= static_cast<double>(cluster_cells.size());
    avg_c /= static_cast<double>(cluster_cells.size());

    double min_cost_sq = std::numeric_limits<double>::max();
    for (const auto& pt : cluster_cells) {
        if (!grid.isPassable(pt)) continue;
        const double dr = pt.r - avg_r;
        const double dc = pt.c - avg_c;
        const double d2 = dr * dr + dc * dc;
        const uint8_t clr = grid.getClearance(pt);
        const double clr_penalty = (clr <= 1) ? 200.0 : ((clr == 2) ? 20.0 : 0.0);
        const double cost = d2 + clr_penalty;
        if (cost < min_cost_sq) {
            min_cost_sq = cost;
            best_p = pt;
        }
    }
    return findSafeViewpoint(best_p);
}

// ============================================================================
// 8-Connected Component Extraction & Size Partitioning
// ============================================================================

std::vector<FrontierCluster> FrontierManager::detectAndCluster(const GridMap& grid) const {
    std::vector<FrontierCluster> clusters;
    const int R = grid.rows;
    const int C = grid.cols;
    if (R <= 0 || C <= 0) return clusters;

    // Identify all frontier cells across the grid
    std::vector<bool> is_frontier(static_cast<size_t>(R * C), false);
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            const Point p{r, c};
            if (isFrontierCell(grid, p)) {
                is_frontier[static_cast<size_t>(grid.index(p))] = true;
            }
        }
    }

    std::vector<bool> visited(static_cast<size_t>(R * C), false);
    int cluster_id = 0;

    // 8-connected clustering offsets for contiguous frontier fronts
    static const int dr8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    static const int dc8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            const int start_idx = grid.index(r, c);
            if (!is_frontier[static_cast<size_t>(start_idx)] || visited[static_cast<size_t>(start_idx)]) continue;

            FrontierCluster cluster;
            cluster.id = cluster_id++;

            std::queue<Point> q;
            q.push({r, c});
            visited[static_cast<size_t>(start_idx)] = true;

            while (!q.empty()) {
                const Point curr = q.front();
                q.pop();
                cluster.cells.push_back(curr);

                for (int i = 0; i < 8; ++i) {
                    const int nr = curr.r + dr8[i];
                    const int nc = curr.c + dc8[i];
                    if (!grid.inBounds(nr, nc)) continue;
                    const int nidx = grid.index(nr, nc);
                    if (is_frontier[static_cast<size_t>(nidx)] && !visited[static_cast<size_t>(nidx)]) {
                        visited[static_cast<size_t>(nidx)] = true;
                        q.push({nr, nc});
                    }
                }
            }

            if (cluster.size() >= min_cluster_size) {
                // Verification: ensure cluster has at least one cell with clearance >= 2,
                // guaranteeing a physical robot can approach safely
                uint8_t max_clr = 0;
                for (const auto& cell : cluster.cells) {
                    max_clr = std::max(max_clr, grid.getClearance(cell));
                }
                if (max_clr <= 1) {
                    continue; // Discard wall-trapped pseudo-frontier cluster
                }

                // Partition oversized frontier components into localized viewpoint clusters
                if (max_cluster_size > 0 && cluster.size() > max_cluster_size) {
                    for (size_t start = 0; start < cluster.cells.size(); start += static_cast<size_t>(max_cluster_size)) {
                        const size_t end = std::min(start + static_cast<size_t>(max_cluster_size), cluster.cells.size());
                        if (end - start >= static_cast<size_t>(min_cluster_size)) {
                            FrontierCluster sub;
                            sub.id = cluster_id++;
                            sub.cells.assign(cluster.cells.begin() + static_cast<ptrdiff_t>(start),
                                             cluster.cells.begin() + static_cast<ptrdiff_t>(end));
                            sub.representative = computeMedoid(grid, sub.cells);
                            clusters.push_back(std::move(sub));
                        }
                    }
                } else {
                    cluster.representative = computeMedoid(grid, cluster.cells);
                    clusters.push_back(std::move(cluster));
                }
            }
        }
    }

    return clusters;
}

} // namespace htrip
