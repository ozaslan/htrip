#include "htrip/tour_planner.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>

namespace htrip {

// ============================================================================
// Multi-Robot Fleet Tour Planning (m-TSP)
// ============================================================================

void MTSPPlanner::planFleetTours(
    const std::vector<Point>& robot_positions,
    const std::vector<Point>& frontier_points,
    const std::vector<dist_t>& dist_matrix,
    std::vector<RobotTour>& out_tours,
    size_t max_frontiers_per_robot) {

    const size_t N = robot_positions.size();
    const size_t K = frontier_points.size();
    const size_t M = N + K;

    out_tours.resize(N);
    for (size_t i = 0; i < N; ++i) {
        out_tours[i].robot_id = static_cast<int>(i);
        out_tours[i].clear();
    }

    if (N == 0 || K == 0 || dist_matrix.size() < M * M) {
        return;
    }

    // 1. Capacity-Constrained Partitioning:
    // Compute dynamic capacity cap per robot: cap = min(max_cap, max(1, (K + N - 1) / N + 4)).
    // The +4 slack allows geographically natural clustering without fragmenting nearby frontiers.
    const size_t dynamic_cap = std::max(size_t{1}, (K + N - 1) / N + 4);
    const size_t cap = std::min(max_frontiers_per_robot, dynamic_cap);

    std::vector<std::vector<int>> assigned(N);
    for (size_t i = 0; i < N; ++i) {
        assigned[i].reserve(std::max(cap, size_t{1}));
    }

    // Assign each frontier to its closest reachable robot that has not exceeded capacity
    for (size_t f = 0; f < K; ++f) {
        dist_t best_dist = INF;
        int best_r = -1;

        for (size_t r = 0; r < N; ++r) {
            if (assigned[r].size() >= cap) continue;
            const dist_t d = dist_matrix[r * M + (N + f)];
            if (d < best_dist) {
                best_dist = d;
                best_r = static_cast<int>(r);
            }
        }

        // Fallback: if all robots within capacity are unreachable or full, relax cap
        if (best_r == -1) {
            for (size_t r = 0; r < N; ++r) {
                const dist_t d = dist_matrix[r * M + (N + f)];
                if (d < best_dist) {
                    best_dist = d;
                    best_r = static_cast<int>(r);
                }
            }
        }

        if (best_r >= 0 && best_dist < INF) {
            assigned[static_cast<size_t>(best_r)].push_back(static_cast<int>(f));
        }
    }

    // 1b. Coverage Fill for Idle Robots:
    // Nearest-robot assignment can leave trailing or peripheral robots empty if every
    // frontier is slightly closer to another teammate. To maintain 100% fleet productivity,
    // idle robots take the least-loaded reachable frontier (shared target).
    // If topological geodesic distance is INF due to an incomplete map snapshot,
    // fallback to Euclidean distance to steer the idle robot toward the unexplored territory.
    std::vector<int> load(K, 0);
    for (size_t r = 0; r < N; ++r) {
        for (int f : assigned[r]) {
            load[static_cast<size_t>(f)]++;
        }
    }

    for (size_t r = 0; r < N; ++r) {
        if (!assigned[r].empty()) continue;

        int best_f = -1;
        dist_t best_d = INF;
        int best_load = std::numeric_limits<int>::max();

        for (size_t f = 0; f < K; ++f) {
            const dist_t d = dist_matrix[r * M + (N + f)];
            if (d >= INF) continue;
            if (load[f] < best_load || (load[f] == best_load && d < best_d)) {
                best_load = load[f];
                best_d = d;
                best_f = static_cast<int>(f);
            }
        }

        // Euclidean fallback if topological connectivity is not yet established
        if (best_f < 0) {
            int64_t best_e2 = std::numeric_limits<int64_t>::max();
            best_load = std::numeric_limits<int>::max();
            for (size_t f = 0; f < K; ++f) {
                const int64_t dr = static_cast<int64_t>(robot_positions[r].r) -
                                   static_cast<int64_t>(frontier_points[f].r);
                const int64_t dc = static_cast<int64_t>(robot_positions[r].c) -
                                   static_cast<int64_t>(frontier_points[f].c);
                const int64_t e2 = dr * dr + dc * dc;
                if (load[f] < best_load || (load[f] == best_load && e2 < best_e2)) {
                    best_load = load[f];
                    best_e2 = e2;
                    best_f = static_cast<int>(f);
                }
            }
        }

        if (best_f >= 0) {
            assigned[r].push_back(best_f);
            load[static_cast<size_t>(best_f)]++;
        }
    }

    // 2. Compute Open TSP Tour for each robot
    for (size_t r = 0; r < N; ++r) {
        if (assigned[r].empty()) continue;

        const std::vector<int> ordered_f = solveOpenTSP(
            static_cast<int>(r), assigned[r], N, M, dist_matrix);

        out_tours[r].waypoints.reserve(ordered_f.size());
        uint32_t total_cost = 0;
        int prev_idx = static_cast<int>(r);

        for (int f : ordered_f) {
            out_tours[r].waypoints.push_back(frontier_points[static_cast<size_t>(f)]);
            const int curr_idx = static_cast<int>(N + static_cast<size_t>(f));
            const dist_t leg = dist_matrix[static_cast<size_t>(prev_idx) * M + static_cast<size_t>(curr_idx)];
            if (leg < INF) total_cost += leg;
            prev_idx = curr_idx;
        }

        out_tours[r].total_tour_cost = total_cost;
    }
}

// ============================================================================
// Open TSP Solver: Nearest Neighbor + 2-Opt Local Search
// ============================================================================

std::vector<int> MTSPPlanner::solveOpenTSP(
    int robot_idx,
    const std::vector<int>& assigned_frontier_indices,
    size_t num_robots,
    size_t total_points,
    const std::vector<dist_t>& dist_matrix) {

    const size_t m = assigned_frontier_indices.size();
    if (m == 0) return {};
    if (m == 1) return {assigned_frontier_indices.front()};

    const size_t M = total_points;
    auto get_dist = [&](int u, int v) -> dist_t {
        if (u < 0 || v < 0 || static_cast<size_t>(u) >= M || static_cast<size_t>(v) >= M) return INF;
        return dist_matrix[static_cast<size_t>(u) * M + static_cast<size_t>(v)];
    };

    // Stage 1: Greedy Nearest Neighbor tour construction from robot pose
    std::vector<int> tour;
    tour.reserve(m);
    std::vector<bool> visited(m, false);
    int curr_matrix_idx = robot_idx;

    for (size_t step = 0; step < m; ++step) {
        dist_t best_d = INF;
        int best_k = -1;

        for (size_t k = 0; k < m; ++k) {
            if (visited[k]) continue;
            const int f_idx = assigned_frontier_indices[k];
            const int f_matrix_idx = static_cast<int>(num_robots + static_cast<size_t>(f_idx));
            const dist_t d = get_dist(curr_matrix_idx, f_matrix_idx);
            if (d < best_d) {
                best_d = d;
                best_k = static_cast<int>(k);
            }
        }

        if (best_k == -1) {
            // Pick first remaining unvisited target
            for (size_t k = 0; k < m; ++k) {
                if (!visited[k]) {
                    best_k = static_cast<int>(k);
                    break;
                }
            }
        }

        if (best_k >= 0) {
            visited[static_cast<size_t>(best_k)] = true;
            tour.push_back(assigned_frontier_indices[static_cast<size_t>(best_k)]);
            curr_matrix_idx = static_cast<int>(num_robots + static_cast<size_t>(tour.back()));
        }
    }

    if (tour.size() <= 2) return tour;

    // Helper: evaluate total open tour travel cost
    auto eval_tour_cost = [&](const std::vector<int>& t) -> int64_t {
        if (t.empty()) return 0;
        int64_t cost = 0;
        int prev = robot_idx;
        for (int f : t) {
            const int cur = static_cast<int>(num_robots + static_cast<size_t>(f));
            const dist_t d = get_dist(prev, cur);
            cost += (d < INF ? d : 1000000);
            prev = cur;
        }
        return cost;
    };

    // Stage 2: 2-Opt local search improvement for open tour
    int64_t best_cost = eval_tour_cost(tour);
    bool improved = true;
    int iterations = 0;

    while (improved && iterations < 20) {
        improved = false;
        iterations++;

        for (size_t i = 0; i < tour.size(); ++i) {
            for (size_t j = i + 1; j < tour.size(); ++j) {
                std::reverse(tour.begin() + static_cast<ptrdiff_t>(i),
                             tour.begin() + static_cast<ptrdiff_t>(j + 1));
                const int64_t new_cost = eval_tour_cost(tour);

                if (new_cost < best_cost) {
                    best_cost = new_cost;
                    improved = true;
                } else {
                    // Revert 2-opt segment reversal
                    std::reverse(tour.begin() + static_cast<ptrdiff_t>(i),
                                 tour.begin() + static_cast<ptrdiff_t>(j + 1));
                }
            }
        }
    }

    return tour;
}

} // namespace htrip
