#include "htrip/fuel_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include <algorithm>
#include <cmath>

namespace htrip {

// ============================================================================
// Initial FIS Cache Construction
// ============================================================================

void FuelCacheOracle::build(const GridMap& grid, const std::vector<Point>& viewpoints) {
    viewpoints_ = viewpoints;
    const size_t k = viewpoints_.size();
    cost_matrix_.assign(k * k, INF);

    const auto t0 = std::chrono::high_resolution_clock::now();

    // Compute dense all-pairs matrix initially via k single-source BFS sweeps
    for (size_t i = 0; i < k; ++i) {
        cost_matrix_[i * k + i] = 0;
        const auto row = BfsOracle::queryOneToAll(grid, viewpoints_[i], viewpoints_);
        for (size_t j = 0; j < k; ++j) {
            cost_matrix_[i * k + j] = row[j];
        }
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    last_update_us_ = std::chrono::duration<double, std::micro>(t1 - t0).count();
    last_bfs_count_ = k;
}

// ============================================================================
// Event-Driven Sensor Updates & Staleness Generation
// ============================================================================

void FuelCacheOracle::onCellOpened(
    Point cell_opened, const GridMap& grid, const std::vector<Point>& viewpoints) {
    onSensorUpdate(cell_opened, sensor_radius_, grid, viewpoints);
}

void FuelCacheOracle::onSensorUpdate(
    Point sensor_pos, int box_radius, const GridMap& grid, const std::vector<Point>& viewpoints) {

    const auto t0 = std::chrono::high_resolution_clock::now();

    // Preserve previous viewpoint set and cost matrix for cache adoption
    const auto old_vps = std::move(viewpoints_);
    const auto old_matrix = std::move(cost_matrix_);
    const size_t old_k = old_vps.size();

    viewpoints_ = viewpoints;
    const size_t k = viewpoints_.size();
    cost_matrix_.assign(k * k, INF);
    for (size_t i = 0; i < k; ++i) cost_matrix_[i * k + i] = 0;

    // Map new viewpoint indices to old viewpoint indices to preserve surviving pairwise entries
    std::vector<int> new_to_old(k, -1);
    if (old_k > 0 && old_matrix.size() == old_k * old_k) {
        for (size_t ni = 0; ni < k; ++ni) {
            for (size_t oi = 0; oi < old_k; ++oi) {
                if (viewpoints_[ni] == old_vps[oi]) {
                    new_to_old[ni] = static_cast<int>(oi);
                    break;
                }
            }
        }
        // Copy over previously cached distance entries between surviving viewpoints
        for (size_t ni = 0; ni < k; ++ni) {
            const int oi = new_to_old[ni];
            if (oi < 0) continue;
            for (size_t nj = 0; nj < k; ++nj) {
                const int oj = new_to_old[nj];
                if (oj < 0) continue;
                cost_matrix_[ni * k + nj] = old_matrix[static_cast<size_t>(oi) * old_k + static_cast<size_t>(oj)];
            }
        }
    }

    // Identify which viewpoints fall within the sensor bounding box (AABB)
    // or are newly created viewpoints. In FUEL, only these are marked as modified.
    std::vector<bool> is_modified(k, false);
    size_t modified_count = 0;

    for (size_t i = 0; i < k; ++i) {
        if (new_to_old[i] < 0) {
            is_modified[i] = true;
            modified_count++;
            continue;
        }
        const int dr = std::abs(viewpoints_[i].r - sensor_pos.r);
        const int dc = std::abs(viewpoints_[i].c - sensor_pos.c);
        if (dr <= box_radius && dc <= box_radius) {
            is_modified[i] = true;
            modified_count++;
        }
    }

    last_bfs_count_ = 0;

    // ------------------------------------------------------------------------
    // FUEL Caching Policy (FrontierFinder::updateFrontierCostMatrix):
    // Single-source BFS/A* searches are executed ONLY for modified viewpoints.
    //
    // CRITICAL OBSERVATION / FAILURE MODE:
    // Viewpoint pairs where BOTH viewpoints lie outside the sensor bounding box
    // are completely bypassed and retain their previously cached distance entries.
    // If a newly opened doorway or corridor creates a dramatic shortcut between
    // two distant clusters outside the sensor box, FUEL preserves the old high-cost
    // path! This produces an overestimation error (staleness) in the cost matrix,
    // which misleads downstream TSP tour planning.
    // ------------------------------------------------------------------------
    for (size_t i = 0; i < k; ++i) {
        if (!is_modified[i]) continue;

        // Run single-source BFS from this modified viewpoint to all viewpoints
        const auto row = BfsOracle::queryOneToAll(grid, viewpoints_[i], viewpoints_);
        last_bfs_count_++;

        for (size_t j = 0; j < k; ++j) {
            cost_matrix_[i * k + j] = row[j];
            cost_matrix_[j * k + i] = row[j];
        }
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    last_update_us_ = std::chrono::duration<double, std::micro>(t1 - t0).count();
}

dist_t FuelCacheOracle::queryDistance(size_t i, size_t j) const noexcept {
    const size_t k = viewpoints_.size();
    if (i >= k || j >= k) return INF;
    return cost_matrix_[i * k + j];
}

// ============================================================================
// Staleness Quantification against Ground Truth
// ============================================================================

FuelStalenessReport FuelCacheOracle::evaluateAgainstGroundTruth(
    const std::vector<dist_t>& gt_matrix) const {

    FuelStalenessReport rep;
    const size_t k = viewpoints_.size();
    if (gt_matrix.size() != k * k || cost_matrix_.size() != k * k) return rep;

    rep.total_pairs = 0;
    rep.stale_pairs = 0;
    rep.max_absolute_error = 0;
    rep.max_relative_error = 0.0;
    double sum_rel_err = 0.0;
    rep.bfs_searches_executed = last_bfs_count_;
    rep.update_time_us = last_update_us_;

    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            if (i == j) continue;
            rep.total_pairs++;

            const dist_t d_cache = cost_matrix_[i * k + j];
            const dist_t d_gt = gt_matrix[i * k + j];

            // If the cached distance strictly exceeds ground truth, it is STALE!
            if (d_cache > d_gt && d_gt < INF) {
                rep.stale_pairs++;
                const dist_t abs_err = d_cache - d_gt;
                if (abs_err > rep.max_absolute_error) {
                    rep.max_absolute_error = abs_err;
                }
                const double rel_err = static_cast<double>(abs_err) / static_cast<double>(d_gt);
                if (rel_err > rep.max_relative_error) {
                    rep.max_relative_error = rel_err;
                }
                sum_rel_err += rel_err;
            }
        }
    }

    if (rep.total_pairs > 0) {
        rep.stale_pair_fraction = static_cast<double>(rep.stale_pairs) / static_cast<double>(rep.total_pairs);
        rep.mean_relative_error = sum_rel_err / static_cast<double>(rep.total_pairs);
    }
    return rep;
}

} // namespace htrip
