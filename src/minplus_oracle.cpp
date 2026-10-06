#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include "htrip/dijkstra_oracle.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#ifndef __AVX2__
#error "AVX2 instruction set is required. Compile with -mavx2 or -march=native."
#endif
#include <immintrin.h>

namespace htrip {

int MinPlusOracle::tileIndexAtLevel(int lvl, int tx, int ty) const noexcept {
    // Determine the spatial side length of tiles at hierarchy level lvl: side_lvl = b * g^lvl
    int side_lvl = config.tile_size_b;
    for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
    int num_tiles_x = S / side_lvl;
    if (num_tiles_x <= 0) num_tiles_x = 1;
    return ty * num_tiles_x + tx;
}

int MinPlusOracle::findLcaLevel(Point p1, Point p2) const noexcept {
    const int b = config.tile_size_b;
    int t1_x = p1.c / b;
    int t1_y = p1.r / b;
    int t2_x = p2.c / b;
    int t2_y = p2.r / b;

    // Level 0: both points reside within the exact same base leaf tile
    if (t1_x == t2_x && t1_y == t2_y) return 0;

    // Ascend hierarchy levels: at each step, divide by branching factor g
    // The lowest level where parent coordinates coincide is the Lowest Common Ancestor (LCA).
    for (int lvl = 1; lvl <= H; ++lvl) {
        t1_x /= config.group_w;
        t1_y /= config.group_h;
        t2_x /= config.group_w;
        t2_y /= config.group_h;
        if (t1_x == t2_x && t1_y == t2_y) {
            return lvl;
        }
    }
    return H + 1;
}

bool MinPlusOracle::build(const GridMap& grid) {
    if (!config.isValid()) return false;
    const int b = config.tile_size_b;
    const int max_dim = std::max(grid.rows, grid.cols);
    if (max_dim <= 0) return false;

    // Determine bounding world dimension S as a power-of-g multiple of leaf size b:
    // S = b * g^H >= max(rows, cols)
    S = b;
    while (S < max_dim) {
        S *= config.group_w;
    }

    // Hierarchy depth H: number of recursive internal aggregation levels above the leaves
    H = 0;
    int temp = S / b;
    while (temp > 1) {
        temp /= config.group_w;
        H++;
    }

    num_leaf_tiles_x = S / b;
    num_leaf_tiles_y = S / b;

    // 1. Precompute base leaf tiles in parallel grid layout.
    // Each leaf tile computes its local APSP metric D (dimension b^2 x b^2) and projects
    // boundary distances onto perimeter ports to form summary S (dimension (4b - 4)^2).
    leaves.clear();
    leaves.reserve(static_cast<size_t>(num_leaf_tiles_x * num_leaf_tiles_y));

    for (int ty = 0; ty < num_leaf_tiles_y; ++ty) {
        for (int tx = 0; tx < num_leaf_tiles_x; ++tx) {
            LeafTile leaf(tx, ty, b);
            leaf.recomputeFromGrid(grid);
            leaves.push_back(std::move(leaf));
        }
    }

    // 2. Build internal aggregation levels bottom-up (levels 1 through H).
    // Level lvl aggregates g x g child regions from level lvl - 1.
    levels.clear();
    levels.resize(static_cast<size_t>(H));

    for (int lvl = 1; lvl <= H; ++lvl) {
        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
        int num_tiles_x = S / side_lvl;
        int num_tiles_y = S / side_lvl;
        if (num_tiles_x <= 0) num_tiles_x = 1;
        if (num_tiles_y <= 0) num_tiles_y = 1;

        levels[static_cast<size_t>(lvl - 1)].reserve(static_cast<size_t>(num_tiles_x * num_tiles_y));

        for (int ty = 0; ty < num_tiles_y; ++ty) {
            for (int tx = 0; tx < num_tiles_x; ++tx) {
                InternalNode node(lvl, tx, ty, b, config.group_w, config.group_h);

                // Gather child summaries using zero-copy spans
                const int num_children = config.group_w * config.group_h;
                std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));

                int child_side = side_lvl / config.group_w;
                for (int k = 0; k < num_children; ++k) {
                    int cy = k / config.group_w;
                    int cx = k % config.group_w;
                    int child_tx = tx * config.group_w + cx;
                    int child_ty = ty * config.group_h + cy;

                    if (lvl == 1) {
                        int l_idx = child_ty * num_leaf_tiles_x + child_tx;
                        child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(l_idx)].S;
                    } else {
                        int prev_tiles_x = S / child_side;
                        int p_idx = child_ty * prev_tiles_x + child_tx;
                        child_summaries[static_cast<size_t>(k)] =
                            levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(p_idx)].S_U;
                    }
                }

                // Compose child summaries, insert touching seams, and execute AVX2 min-plus closure
                node.recomputeFromChildren(child_summaries, grid);
                levels[static_cast<size_t>(lvl - 1)].push_back(std::move(node));
            }
        }
    }

    // Structural reset: increment epoch to invalidate any previous incremental cache
    epoch++;
    return true;
}

void MinPlusOracle::onCellOpened(Point p, const GridMap& grid) {
    const int b = config.tile_size_b;
    const int tx = p.c / b;
    const int ty = p.r / b;
    if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) return;

    // 1. Update the containing leaf tile via Conway rank-1 tropical relaxation (Zero BFS)
    LeafTile& leaf = leaves[static_cast<size_t>(leafIndex(tx, ty))];
    leaf.version = ++update_counter;
    auto decreasedPorts = leaf.onCellOpened(p, grid);

    // 2. Identify newly opened 4-connected seam edges incident to p that cross tile boundaries
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};
    std::vector<SeamEdge> seamEdges;
    seamEdges.reserve(4);

    for (int i = 0; i < 4; ++i) {
        int nr = p.r + dr[i];
        int nc = p.c + dc[i];
        if (grid.isPassable(nr, nc)) {
            int ntx = nc / b;
            int nty = nr / b;
            if (ntx != tx || nty != ty) {
                Point np{nr, nc};
                dist_t w = grid.edgeWeight(p, np);
                if (w < INF) {
                    int lca_lvl = findLcaLevel(p, np);
                    if (lca_lvl <= H) {
                        seamEdges.push_back({p, np, lca_lvl, w});
                    }
                }
            }
        }
    }

    // 3. Propagate decreased boundary summaries and seam edges upward along ancestor branch
    propagateUpward(tx, ty, decreasedPorts, seamEdges, grid);
}

void MinPlusOracle::propagateUpward(int leaf_tx, int leaf_ty,
                                    const std::vector<PortDecrease>& leafDecreases,
                                    const std::vector<SeamEdge>& seamEdges,
                                    const GridMap& grid) {
    // If decrease volume is large, delegate to batched upward propagation
    if (leafDecreases.size() > 16 || seamEdges.size() > 4) {
        propagateBatchUpward({{leaf_tx, leaf_ty}}, grid);
        return;
    }
    (void)grid;
    int cur_tx = leaf_tx;
    int cur_ty = leaf_ty;
    std::vector<PortDecrease> curDecreases = leafDecreases;

    // Ascend ancestor chain from level 1 to root level H
    for (int lvl = 1; lvl <= H; ++lvl) {
        int parent_tx = cur_tx / config.group_w;
        int parent_ty = cur_ty / config.group_h;
        int child_quad = (cur_ty % config.group_h) * config.group_w + (cur_tx % config.group_w);

        int p_idx = tileIndexAtLevel(lvl, parent_tx, parent_ty);
        InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];
        int childPortBase = parent.childOffset[static_cast<size_t>(child_quad)];

        // 1. Incorporate decreased child summary edges via rank-1 updates
        for (const auto& dec : curDecreases) {
            int u = childPortBase + dec.p_a;
            int v = childPortBase + dec.p_b;
            parent.updateRank1(u, v, dec.w);
        }

        // 2. Incorporate newly opened inter-child seam edges enclosed at this level
        for (const auto& se : seamEdges) {
            if (se.lca_lvl == lvl) {
                parent.insertSeamEdge(se.p1, se.p2, se.w);
            }
        }

        // 3. Extract new exterior boundary summary decreases
        auto nextDecreases = parent.refreshBoundarySummary();
        curDecreases = std::move(nextDecreases);

        // Keep contiguous child-to-exterior lift cache consistent with D_U
        const bool changed_here = parent.child_to_ext_dirty;
        if (changed_here) {
            parent.refreshChildToExt();
            parent.version = ++update_counter;
        }

        cur_tx = parent_tx;
        cur_ty = parent_ty;

        // Check if any higher-level seam events are still pending
        bool pending_higher_seam = false;
        for (const auto& se : seamEdges) {
            if (se.lca_lvl > lvl) {
                pending_higher_seam = true;
                break;
            }
        }

        // Stabilization: early-exit iff exterior summary stabilized AND no pending higher seams
        if (curDecreases.empty() && !pending_higher_seam) {
            break;
        }
    }
}

void MinPlusOracle::onCellBlocked(Point p, const GridMap& grid) {
    const int b = config.tile_size_b;
    const int tx = p.c / b;
    const int ty = p.r / b;
    if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) return;

    // 1. Non-monotone edge deletion: recompute intra-tile metric D from the grid via local BFS
    LeafTile& leaf = leaves[static_cast<size_t>(leafIndex(tx, ty))];
    leaf.onCellBlocked(p, grid);
    leaf.version = ++update_counter;

    // 2. Bottom-up full rebuild of ancestor internal nodes along the single tree branch to root
    int cur_tx = tx;
    int cur_ty = ty;

    for (int lvl = 1; lvl <= H; ++lvl) {
        int parent_tx = cur_tx / config.group_w;
        int parent_ty = cur_ty / config.group_h;

        int p_idx = tileIndexAtLevel(lvl, parent_tx, parent_ty);
        InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];

        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
        int child_side = side_lvl / config.group_w;
        const int num_children = config.group_w * config.group_h;
        std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));

        for (int k = 0; k < num_children; ++k) {
            int cy = k / config.group_w;
            int cx = k % config.group_w;
            int child_tx = parent_tx * config.group_w + cx;
            int child_ty = parent_ty * config.group_h + cy;

            if (lvl == 1) {
                int l_idx = child_ty * num_leaf_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(l_idx)].S;
            } else {
                int prev_tiles_x = S / child_side;
                int prev_idx = child_ty * prev_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] =
                    levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(prev_idx)].S_U;
            }
        }

        // Recompute parent closure from child summaries + seam edges
        parent.recomputeFromChildren(child_summaries, grid);
        parent.version = ++update_counter;

        cur_tx = parent_tx;
        cur_ty = parent_ty;
    }

    // Structural reset: blockages invalidate the incremental-query cache, forcing full recomputation
    epoch++;
}

void MinPlusOracle::onCellWeightChanged(Point p, dist_t new_weight, GridMap& grid) {
    if (!grid.inBounds(p)) return;
    if (grid.getCellWeight(p) == new_weight) return; // Short-circuit if weight is unchanged
    grid.setCellWeight(p, new_weight);

    const int b = config.tile_size_b;
    const int tx = p.c / b;
    const int ty = p.r / b;
    if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) return;

    // 1. Recompute containing leaf tile
    int leaf_idx = leafIndex(tx, ty);
    LeafTile& leaf = leaves[static_cast<size_t>(leaf_idx)];
    leaf.recomputeFromGrid(grid);
    leaf.version = ++update_counter;

    // 2. Rebuild ancestor node closures bottom-up to the root
    int cur_tx = tx;
    int cur_ty = ty;

    for (int lvl = 1; lvl <= H; ++lvl) {
        int parent_tx = cur_tx / config.group_w;
        int parent_ty = cur_ty / config.group_h;

        int p_idx = tileIndexAtLevel(lvl, parent_tx, parent_ty);
        InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];

        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
        int child_side = side_lvl / config.group_w;
        const int num_children = config.group_w * config.group_h;
        std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));

        for (int k = 0; k < num_children; ++k) {
            int cy = k / config.group_w;
            int cx = k % config.group_w;
            int child_tx = parent_tx * config.group_w + cx;
            int child_ty = parent_ty * config.group_h + cy;

            if (lvl == 1) {
                int l_idx = child_ty * num_leaf_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(l_idx)].S;
            } else {
                int prev_tiles_x = S / child_side;
                int prev_idx = child_ty * prev_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] =
                    levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(prev_idx)].S_U;
            }
        }

        parent.recomputeFromChildren(child_summaries, grid);
        parent.version = ++update_counter;

        cur_tx = parent_tx;
        cur_ty = parent_ty;
    }

    epoch++;
}

void MinPlusOracle::onTileWeightsChanged(int tx, int ty,
                                         const std::vector<std::pair<Point, dist_t>>& cell_weights,
                                         GridMap& grid) {
    if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) return;

    // Apply all cell weight updates to the grid costmap
    for (const auto& [p, w] : cell_weights) {
        if (grid.inBounds(p)) {
            grid.setCellWeight(p, w);
        }
    }

    // Recompute the containing leaf tile ONCE for the entire batch
    int leaf_idx = leafIndex(tx, ty);
    LeafTile& leaf = leaves[static_cast<size_t>(leaf_idx)];
    leaf.recomputeFromGrid(grid);
    leaf.version = ++update_counter;

    int cur_tx = tx;
    int cur_ty = ty;
    const int b = config.tile_size_b;

    // Propagate closures up the ancestor branch ONCE
    for (int lvl = 1; lvl <= H; ++lvl) {
        int parent_tx = cur_tx / config.group_w;
        int parent_ty = cur_ty / config.group_h;

        int p_idx = tileIndexAtLevel(lvl, parent_tx, parent_ty);
        InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];

        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
        int child_side = side_lvl / config.group_w;
        const int num_children = config.group_w * config.group_h;
        std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));

        for (int k = 0; k < num_children; ++k) {
            int cy = k / config.group_w;
            int cx = k % config.group_w;
            int child_tx = parent_tx * config.group_w + cx;
            int child_ty = parent_ty * config.group_h + cy;

            if (lvl == 1) {
                int l_idx = child_ty * num_leaf_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(l_idx)].S;
            } else {
                int prev_tiles_x = S / child_side;
                int prev_idx = child_ty * prev_tiles_x + child_tx;
                child_summaries[static_cast<size_t>(k)] =
                    levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(prev_idx)].S_U;
            }
        }

        parent.recomputeFromChildren(child_summaries, grid);
        parent.version = ++update_counter;

        cur_tx = parent_tx;
        cur_ty = parent_ty;
    }

    epoch++;
}

void MinPlusOracle::propagateBatchUpward(
    std::vector<std::pair<int, int>> dirty_tiles,
    const GridMap& grid) {
    if (dirty_tiles.empty()) return;

    const int b = config.tile_size_b;
    const int gw = config.group_w;
    const int gh = config.group_h;
    const int num_children = gw * gh;

    std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));

    // Ascend hierarchy levels: at each level, rebuild only unique dirty parent tiles
    for (int lvl = 1; lvl <= H; ++lvl) {
        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= gw;
        int child_side = side_lvl / gw;

        // Deduplicate dirty parent coordinates
        std::vector<std::pair<int, int>> dirty_parents;
        dirty_parents.reserve(dirty_tiles.size());

        for (const auto& [cx, cy] : dirty_tiles) {
            int px = cx / gw;
            int py = cy / gh;
            if (std::find(dirty_parents.begin(), dirty_parents.end(), std::make_pair(px, py)) == dirty_parents.end()) {
                dirty_parents.push_back({px, py});
            }
        }

        if (dirty_parents.empty()) break;

        std::vector<std::pair<int, int>> next_dirty_tiles;
        next_dirty_tiles.reserve(dirty_parents.size());

        for (const auto& [ptx, pty] : dirty_parents) {
            int p_idx = tileIndexAtLevel(lvl, ptx, pty);
            if (p_idx < 0 || p_idx >= static_cast<int>(levels[static_cast<size_t>(lvl - 1)].size())) continue;
            InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];

            for (int k = 0; k < num_children; ++k) {
                int cx = ptx * gw + (k % gw);
                int cy = pty * gh + (k / gw);
                if (lvl == 1) {
                    child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(leafIndex(cx, cy))].S;
                } else {
                    int prev_tiles_x = S / child_side;
                    int prev_idx = cy * prev_tiles_x + cx;
                    child_summaries[static_cast<size_t>(k)] = levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(prev_idx)].S_U;
                }
            }

            std::vector<dist_t> old_S_U = parent.S_U;
            parent.recomputeFromChildren(child_summaries, grid);
            parent.version = ++update_counter;

            // Detect if exterior summary changed; only propagating if outer ports changed
            bool su_changed = false;
            if (parent.S_U.size() != old_S_U.size()) {
                su_changed = true;
            } else {
                for (size_t i = 0; i < old_S_U.size(); ++i) {
                    if (old_S_U[i] != parent.S_U[i]) {
                        su_changed = true;
                        break;
                    }
                }
            }

            if (su_changed) {
                next_dirty_tiles.push_back({ptx, pty});
            }
        }

        dirty_tiles = std::move(next_dirty_tiles);
    }

    epoch++;
}

void MinPlusOracle::onBatchLeavesChanged(const std::vector<int>& leaf_indices, const GridMap& grid) {
    if (leaf_indices.empty()) return;

    std::vector<std::pair<int, int>> dirty_tiles;
    dirty_tiles.reserve(leaf_indices.size());

    // 1. Recompute each dirty leaf tile and collect distinct 2D tile coordinates
    for (int l_idx : leaf_indices) {
        if (l_idx < 0 || l_idx >= static_cast<int>(leaves.size())) continue;
        LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        leaf.recomputeFromGrid(grid);
        leaf.version = ++update_counter;

        int tx = l_idx % num_leaf_tiles_x;
        int ty = l_idx / num_leaf_tiles_x;
        dirty_tiles.push_back({tx, ty});
    }

    // 2. Propagate closures upward across the tree
    propagateBatchUpward(std::move(dirty_tiles), grid);
}


void MinPlusOracle::onBatchCellsOpened(const std::vector<Point>& newly_opened, const GridMap& grid) {
    last_batch_stats = {};
    if (newly_opened.empty()) return;

    const int b = config.tile_size_b;
    const int gw = config.group_w;
    const int gh = config.group_h;

    // 1. Group newly opened cells by enclosing base leaf tile
    std::vector<std::vector<Point>> cells_by_leaf(leaves.size());
    std::vector<int> dirty_leaf_list;
    dirty_leaf_list.reserve(newly_opened.size());

    for (Point p : newly_opened) {
        int tx = p.c / b;
        int ty = p.r / b;
        if (tx >= 0 && tx < num_leaf_tiles_x && ty >= 0 && ty < num_leaf_tiles_y) {
            int l_idx = leafIndex(tx, ty);
            if (cells_by_leaf[static_cast<size_t>(l_idx)].empty()) {
                dirty_leaf_list.push_back(l_idx);
            }
            cells_by_leaf[static_cast<size_t>(l_idx)].push_back(p);
        }
    }

    if (dirty_leaf_list.empty()) return;

    // 2. Identify all newly opened inter-tile seam edges crossing leaf boundaries.
    // When a newly opened cell p lies on a leaf perimeter and touches a passable neighbor np
    // in an adjacent leaf tile, this forms a cross-boundary shortcut.
    // The lowest level containing both tiles is lca_lvl = findLcaLevel(p, np).
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};
    std::vector<std::vector<std::vector<SeamEdge>>> seams_by_lvl(static_cast<size_t>(H + 1));
    for (int l = 1; l <= H; ++l) {
        seams_by_lvl[static_cast<size_t>(l)].resize(levels[static_cast<size_t>(l - 1)].size());
    }

    for (Point p : newly_opened) {
        int tx = p.c / b;
        int ty = p.r / b;
        for (int i = 0; i < 4; ++i) {
            int nr = p.r + dr[i];
            int nc = p.c + dc[i];
            if (grid.isPassable(nr, nc)) {
                int ntx = nc / b;
                int nty = nr / b;
                if (ntx != tx || nty != ty) {
                    Point np{nr, nc};
                    dist_t w = grid.edgeWeight(p, np);
                    if (w < INF) {
                        int lca_lvl = findLcaLevel(p, np);
                        if (lca_lvl >= 1 && lca_lvl <= H) {
                            int px = tx;
                            int py = ty;
                            for (int l = 1; l <= lca_lvl; ++l) {
                                px /= gw;
                                py /= gh;
                            }
                            int p_idx = tileIndexAtLevel(lca_lvl, px, py);
                            if (p_idx >= 0 && p_idx < static_cast<int>(seams_by_lvl[static_cast<size_t>(lca_lvl)].size())) {
                                seams_by_lvl[static_cast<size_t>(lca_lvl)][static_cast<size_t>(p_idx)].push_back({p, np, lca_lvl, w});
                            }
                        }
                    }
                }
            }
        }
    }

    // 3. Process dirty leaf tiles algebraically via Conway rank-1 updates (Zero BFS / Zero Dijkstra).
    // Theorem 1 (Monotone Exactness):
    // Because edge additions are monotonically non-increasing (w_new < w_old), the Kleene closure
    // of the perturbed leaf is exactly maintained by relaxing over dirty boundary pivots K.
    if (levels.empty() || H == 0) {
        int pivots = 0;
        for (int l_idx : dirty_leaf_list) {
            LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
            leaf.onBatchCellsOpened(cells_by_leaf[static_cast<size_t>(l_idx)], grid);
            leaf.version = ++update_counter;
            pivots += leaf.last_pivot_count;
        }
        last_batch_stats.n_touched_leaves = static_cast<int>(dirty_leaf_list.size());
        last_batch_stats.n_pivots = pivots;
        return;
    }

    std::vector<std::vector<ChildDecrease>> cur_level_decreases(levels[0].size());
    std::vector<int> cur_active_parents;
    int pivots = 0;

    for (int l_idx : dirty_leaf_list) {
        LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        leaf.version = ++update_counter;
        auto leaf_decs = leaf.onBatchCellsOpened(cells_by_leaf[static_cast<size_t>(l_idx)], grid);
        pivots += leaf.last_pivot_count;
        if (!leaf_decs.empty()) {
            int tx = l_idx % num_leaf_tiles_x;
            int ty = l_idx / num_leaf_tiles_x;
            int ptx = tx / gw;
            int pty = ty / gh;
            int quad = (ty % gh) * gw + (tx % gw);
            int p_idx = tileIndexAtLevel(1, ptx, pty);

            if (p_idx >= 0 && p_idx < static_cast<int>(cur_level_decreases.size())) {
                if (cur_level_decreases[static_cast<size_t>(p_idx)].empty()) {
                    cur_active_parents.push_back(p_idx);
                }
                auto& dec_list = cur_level_decreases[static_cast<size_t>(p_idx)];
                for (const auto& ld : leaf_decs) {
                    dec_list.push_back({quad, ld.p_a, ld.p_b, ld.w});
                }
            }
        }
    }

    int dirty_anc = 0;

    // 4. Hierarchical Conway Upward Propagation across internal levels (Zero Floyd-Warshall).
    // Propagates decreased port distances and enclosed seam edges bottom-up from level 1 to root H.
    for (int lvl = 1; lvl <= H; ++lvl) {
        const auto& level_seams = seams_by_lvl[static_cast<size_t>(lvl)];

        // Enqueue any parents that enclose newly formed seam edges at this level
        for (size_t p_idx = 0; p_idx < level_seams.size(); ++p_idx) {
            if (!level_seams[p_idx].empty()) {
                if (cur_level_decreases[p_idx].empty()) {
                    cur_active_parents.push_back(static_cast<int>(p_idx));
                }
            }
        }

        if (cur_active_parents.empty()) {
            // Check if any higher level has pending seam edges
            bool has_higher_seams = false;
            for (int hl = lvl + 1; hl <= H; ++hl) {
                for (const auto& s_list : seams_by_lvl[static_cast<size_t>(hl)]) {
                    if (!s_list.empty()) {
                        has_higher_seams = true;
                        break;
                    }
                }
                if (has_higher_seams) break;
            }
            if (!has_higher_seams) break;
        }

        dirty_anc += static_cast<int>(cur_active_parents.size());

        int next_lvl_size = (lvl < H) ? static_cast<int>(levels[static_cast<size_t>(lvl)].size()) : 0;
        std::vector<std::vector<ChildDecrease>> next_level_decreases(static_cast<size_t>(next_lvl_size));
        std::vector<int> next_active_parents;

        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= gw;
        int num_tiles_x_lvl = S / side_lvl;
        if (num_tiles_x_lvl <= 0) num_tiles_x_lvl = 1;

        // Apply Conway perturbations to each active parent node at level lvl
        for (int p_idx : cur_active_parents) {
            if (p_idx < 0 || p_idx >= static_cast<int>(levels[static_cast<size_t>(lvl - 1)].size())) continue;
            InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];

            const auto& child_decs = cur_level_decreases[static_cast<size_t>(p_idx)];
            const auto& seam_edges = level_seams[static_cast<size_t>(p_idx)];

            bool is_root = (lvl == H);
            auto ext_decs = parent.applyConwayPerturbations(child_decs, seam_edges, is_root);

            // If exterior boundary decreased and not at root, route decreases to grandparent level
            if (!ext_decs.empty() && lvl < H) {
                int ptx = p_idx % num_tiles_x_lvl;
                int pty = p_idx / num_tiles_x_lvl;
                int gptx = ptx / gw;
                int gpty = pty / gh;
                int quad = (pty % gh) * gw + (ptx % gw);
                int gp_idx = tileIndexAtLevel(lvl + 1, gptx, gpty);

                if (gp_idx >= 0 && gp_idx < next_lvl_size) {
                    if (next_level_decreases[static_cast<size_t>(gp_idx)].empty()) {
                        next_active_parents.push_back(gp_idx);
                    }
                    auto& next_list = next_level_decreases[static_cast<size_t>(gp_idx)];
                    for (const auto& ed : ext_decs) {
                        next_list.push_back({quad, ed.p_a, ed.p_b, ed.w});
                    }
                }
            }
        }

        cur_level_decreases = std::move(next_level_decreases);
        cur_active_parents = std::move(next_active_parents);
    }

    last_batch_stats.n_touched_leaves = static_cast<int>(dirty_leaf_list.size());
    last_batch_stats.n_dirty_ancestors = dirty_anc;
    last_batch_stats.n_pivots = pivots;

    epoch++;
}

void MinPlusOracle::liftPoint(Point p, LiftedViewpoint& lifted) const {
    const int b = config.tile_size_b;
    const int p_tx = p.c / b;
    const int p_ty = p.r / b;

    if (lifted.X.size() != static_cast<size_t>(H + 1)) {
        lifted.X.resize(static_cast<size_t>(H + 1));
        lifted.active.resize(static_cast<size_t>(H + 1));
        lifted.parent_node_idx.resize(static_cast<size_t>(H + 1));
        lifted.child_quad.resize(static_cast<size_t>(H + 1));
        lifted.child_offset.resize(static_cast<size_t>(H + 1));
    }
    std::fill_n(lifted.parent_node_idx.data(), H + 1, -1);
    std::fill_n(lifted.child_quad.data(), H + 1, -1);
    std::fill_n(lifted.child_offset.data(), H + 1, 0);

    // Level 0: extract distances from point p to all perimeter ports of its base leaf tile from D
    const LeafTile& leaf = leaves[static_cast<size_t>(leafIndex(p_tx, p_ty))];
    int u = leaf.localIdx(p.r, p.c);
    int nP = static_cast<int>(leaf.ports.size());
    if (lifted.X[0].size() != static_cast<size_t>(nP)) {
        lifted.X[0].resize(static_cast<size_t>(nP));
    }
    lifted.active[0].clear();

    for (int i = 0; i < nP; ++i) {
        int pCell = leaf.localIdx(leaf.ports[static_cast<size_t>(i)].r,
                                  leaf.ports[static_cast<size_t>(i)].c);
        dist_t d = leaf.D[static_cast<size_t>(u * (b * b) + pCell)];
        lifted.X[0][static_cast<size_t>(i)] = d;
        if (d < INF) {
            lifted.active[0].push_back(i);
        }
    }
    std::sort(lifted.active[0].begin(), lifted.active[0].end(), [&](int a, int b_idx) {
        return lifted.X[0][static_cast<size_t>(a)] < lifted.X[0][static_cast<size_t>(b_idx)];
    });

    int cur_tx = p_tx, cur_ty = p_ty;
    for (int lvl = 1; lvl <= H; ++lvl) {
        int parent_tx = cur_tx / config.group_w;
        int parent_ty = cur_ty / config.group_h;
        int child_quad = (cur_ty % config.group_h) * config.group_w + (cur_tx % config.group_w);

        int p_idx = tileIndexAtLevel(lvl, parent_tx, parent_ty);
        const InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];
        int childOffset = parent.childOffset[static_cast<size_t>(child_quad)];
        int nExt = static_cast<int>(parent.extPorts.size());

        lifted.parent_node_idx[static_cast<size_t>(lvl)] = p_idx;
        lifted.child_quad[static_cast<size_t>(lvl)] = child_quad;
        lifted.child_offset[static_cast<size_t>(lvl)] = childOffset;

        if (lifted.X[static_cast<size_t>(lvl)].size() != static_cast<size_t>(nExt)) {
            lifted.X[static_cast<size_t>(lvl)].resize(static_cast<size_t>(nExt));
        }
        std::fill_n(lifted.X[static_cast<size_t>(lvl)].data(), nExt, INF);
        lifted.active[static_cast<size_t>(lvl)].clear();

        const auto& prev_active = lifted.active[static_cast<size_t>(lvl - 1)];
        const auto& submat = parent.child_to_ext[static_cast<size_t>(child_quad)];

        for (int cp : prev_active) {
            dist_t d_child = lifted.X[static_cast<size_t>(lvl - 1)][static_cast<size_t>(cp)];
            const dist_t* row = &submat[static_cast<size_t>(cp * nExt)];

            for (int e = 0; e < nExt; ++e) {
                dist_t d_ext = row[e];
                if (d_ext < INF) {
                    dist_t cand = add_dist(d_child, d_ext);
                    if (cand < lifted.X[static_cast<size_t>(lvl)][static_cast<size_t>(e)]) {
                        lifted.X[static_cast<size_t>(lvl)][static_cast<size_t>(e)] = cand;
                    }
                }
            }
        }

        for (int e = 0; e < nExt; ++e) {
            if (lifted.X[static_cast<size_t>(lvl)][static_cast<size_t>(e)] < INF) {
                lifted.active[static_cast<size_t>(lvl)].push_back(e);
            }
        }
        std::sort(lifted.active[static_cast<size_t>(lvl)].begin(),
                  lifted.active[static_cast<size_t>(lvl)].end(),
                  [&](int a, int b_idx) {
                      return lifted.X[static_cast<size_t>(lvl)][static_cast<size_t>(a)]
                           < lifted.X[static_cast<size_t>(lvl)][static_cast<size_t>(b_idx)];
                  });

        cur_tx = parent_tx;
        cur_ty = parent_ty;
    }
}

// ============================================================================
// Single-Pair Query Evaluation: Lowest Common Ancestor (LCA) Meet
// ============================================================================

dist_t MinPlusOracle::queryDistance(Point s, Point t, const GridMap& grid) const {
    return queryDistance(s, t, grid, s_lift_scratch, t_lift_scratch);
}

dist_t MinPlusOracle::queryDistance(Point s, Point t, const GridMap& grid,
                                    LiftedViewpoint& X_s, LiftedViewpoint& Y_t) const {
    // 1. Boundary & passability verification
    if (!grid.isPassable(s) || !grid.isPassable(t)) return INF;
    if (s == t) return 0;

    const int b = config.tile_size_b;
    const int s_tx = s.c / b;
    const int s_ty = s.r / b;
    const int t_tx = t.c / b;
    const int t_ty = t.r / b;

    if (s_tx < 0 || s_tx >= num_leaf_tiles_x || s_ty < 0 || s_ty >= num_leaf_tiles_y) return INF;
    if (t_tx < 0 || t_tx >= num_leaf_tiles_x || t_ty < 0 || t_ty >= num_leaf_tiles_y) return INF;

    dist_t bestGlobal = INF;
    int lca_lvl = 1;

    // 2. Intra-leaf Fast Path:
    // When both endpoints lie within the same base leaf tile, evaluate the exact
    // intra-tile shortest path D(u, v) in O(1) table lookup time. This serves
    // as an immediate upper bound for any path crossing external tile seams.
    if (s_tx == t_tx && s_ty == t_ty) {
        const LeafTile& leaf = leaves[static_cast<size_t>(leafIndex(s_tx, s_ty))];
        const int u = leaf.localIdx(s.r, s.c);
        const int v = leaf.localIdx(t.r, t.c);
        bestGlobal = leaf.D[static_cast<size_t>(u * (b * b) + v)];
    } else {
        // Find the lowest common ancestor level where s and t share a common parent node
        lca_lvl = findLcaLevel(s, t);
        if (lca_lvl > H) {
            // Disconnected in capped hierarchy (e.g., Flat-BRICK baseline where H=1)
            return INF;
        }
    }

    // 3. Upward viewpoint lifting:
    // Lift both endpoints s and t through the hierarchy to extract distance vectors
    // to perimeter ports at all ancestor levels.
    liftPoint(s, X_s);
    liftPoint(t, Y_t);

    // 4. Meet evaluation at common ancestors W >= LCA:
    // For every level from LCA up to root H, s and t map to corresponding child blocks
    // of a common parent node W. Any shortest path traversing exterior seams must pass
    // through interface boundary ports of some ancestor node W.
    // The candidate distance via interface ports p_s and p_t is:
    //   cand = X_s(lvl-1)[p_s] + D_W(v_s, v_t) + Y_t(lvl-1)[p_t]
    for (int lvl = lca_lvl; lvl <= H; ++lvl) {
        const int p_idx_s = X_s.parent_node_idx[static_cast<size_t>(lvl)];
        const int p_idx_t = Y_t.parent_node_idx[static_cast<size_t>(lvl)];
        if (p_idx_s != p_idx_t || p_idx_s < 0) continue;

        const InternalNode& W = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx_s)];
        const int s_offset = X_s.child_offset[static_cast<size_t>(lvl)];
        const int t_offset = Y_t.child_offset[static_cast<size_t>(lvl)];

        const int nV = static_cast<int>(W.V_U.size());
        const auto& X_vec = X_s.X[static_cast<size_t>(lvl - 1)];
        const auto& Y_vec = Y_t.X[static_cast<size_t>(lvl - 1)];

        const auto& active_s = X_s.active[static_cast<size_t>(lvl - 1)];
        const auto& active_t = Y_t.active[static_cast<size_t>(lvl - 1)];
        if (active_s.empty() || active_t.empty()) continue;

        // Bounding check: if the minimum possible port-to-endpoint distance sum
        // already equals or exceeds bestGlobal, no port combination at this level can improve it.
        if (add_dist(X_vec[static_cast<size_t>(active_s[0])],
                     Y_vec[static_cast<size_t>(active_t[0])]) >= bestGlobal) {
            continue;
        }

        for (int p_s : active_s) {
            const dist_t xi = X_vec[static_cast<size_t>(p_s)];
            if (xi >= bestGlobal) break; // Sorted ascending: further ports will also exceed bestGlobal

            const int vi = s_offset + p_s;
            const dist_t* row = &W.D_U[static_cast<size_t>(vi * nV)];

            for (int p_t : active_t) {
                const dist_t yj = Y_vec[static_cast<size_t>(p_t)];
                if (add_dist(xi, yj) >= bestGlobal) break; // EARLY BREAK: active_t is sorted ascending

                const dist_t dij = row[static_cast<size_t>(t_offset + p_t)];
                if (dij < INF) {
                    const dist_t cand = add_dist(xi, dij, yj);
                    if (cand < bestGlobal) {
                        bestGlobal = cand;
                    }
                }
            }
        }
    }

    return bestGlobal;
}

// ============================================================================
// Batched All-Pairs Matrix Query: Vectorized Multi-Level Meet
// ============================================================================

void MinPlusOracle::queryAllPairs(
    const std::vector<Point>& viewpoints,
    const GridMap& grid,
    std::vector<dist_t>& matrix,
    HTripQueryScratch& scratch) const {

    const size_t k = viewpoints.size();
    if (matrix.size() != k * k) {
        matrix.resize(k * k);
    }
    std::fill(matrix.begin(), matrix.end(), INF);

    if (k == 0) return;

    // 1. Hybrid dispatch:
    // If a hybrid threshold is configured and k is sufficiently small, fall back
    // to search baselines. By default, hybrid_k_threshold == 0 (pure H-TRIP is active).
    if (config.hybrid_k_threshold > 0 && static_cast<int>(k) <= config.hybrid_k_threshold) {
        if (grid.hasWeights()) {
            matrix = DijkstraOracle::queryAllPairs(grid, viewpoints);
        } else {
            matrix = BfsOracle::queryAllPairs(grid, viewpoints);
        }
        return;
    }

    // Set self-distances to zero for all traversable viewpoints
    for (size_t i = 0; i < k; ++i) {
        if (grid.isPassable(viewpoints[i])) {
            matrix[i * k + i] = 0;
        }
    }

    const int b = config.tile_size_b;
    const int gw = config.group_w;
    const int gh = config.group_h;
    const int num_leaves = num_leaf_tiles_x * num_leaf_tiles_y;

    // 2. Zero-Heap Allocation Scratch Invariant:
    // Ensure reusable scratch vectors have sufficient capacity for k viewpoints,
    // tree height H, and all leaf/internal tiles without any runtime reallocations.
    scratch.ensureCapacity(k, H, static_cast<size_t>(num_leaves));
    std::fill_n(scratch.closed.data(), k * k, static_cast<uint8_t>(0));
    for (size_t i = 0; i < k; ++i) {
        scratch.closed[i * k + i] = 1;
    }

    for (int l_idx = 0; l_idx < num_leaves; ++l_idx) {
        scratch.leaf_vps[static_cast<size_t>(l_idx)].clear();
    }

    auto& vp_data = scratch.vp_data;

    // 3. Level 0: Leaf tile grouping & port distance extraction
    // Assign each traversable viewpoint to its base leaf tile and extract its distance
    // vector to all perimeter ports from the precomputed intra-leaf matrix D.
    for (size_t i = 0; i < k; ++i) {
        const Point p = viewpoints[i];
        vp_data[i].pt = p;
        if (!grid.isPassable(p)) continue;

        const int tx = p.c / b;
        const int ty = p.r / b;
        if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) continue;

        vp_data[i].leaf_tx = tx;
        vp_data[i].leaf_ty = ty;
        const int l_idx = leafIndex(tx, ty);
        scratch.leaf_vps[static_cast<size_t>(l_idx)].push_back(static_cast<int>(i));

        const LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        const int u = leaf.localIdx(p.r, p.c);
        const int nP = static_cast<int>(leaf.ports.size());

        vp_data[i].X[0].assign(static_cast<size_t>(nP), INF);
        vp_data[i].active[0].clear();

        for (int port = 0; port < nP; ++port) {
            const int pCell = leaf.localIdx(leaf.ports[static_cast<size_t>(port)].r,
                                            leaf.ports[static_cast<size_t>(port)].c);
            const dist_t d = leaf.D[static_cast<size_t>(u * (b * b) + pCell)];
            vp_data[i].X[0][static_cast<size_t>(port)] = d;
            if (d < INF) vp_data[i].active[0].push_back(port);
        }
        std::sort(vp_data[i].active[0].begin(), vp_data[i].active[0].end(),
                  [&](int a, int c) { return vp_data[i].X[0][static_cast<size_t>(a)] < vp_data[i].X[0][static_cast<size_t>(c)]; });
    }

    // 4. Intra-leaf initial pairwise evaluations + Manhattan lower-bound closing:
    // Pairs within the same leaf evaluate exact shortest paths via D(u, v).
    // If the path length matches the Manhattan distance |r1-r2| + |c1-c2|, it is
    // physically impossible for any path leaving the tile and re-entering to be shorter.
    // Such pairs are marked closed, completely bypassing all higher-level ancestor meets.
    for (int l_idx = 0; l_idx < num_leaves; ++l_idx) {
        const auto& vps = scratch.leaf_vps[static_cast<size_t>(l_idx)];
        if (vps.size() < 2) continue;
        const LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        for (size_t a = 0; a < vps.size(); ++a) {
            const int i = vps[a];
            const int u = leaf.localIdx(vp_data[static_cast<size_t>(i)].pt.r, vp_data[static_cast<size_t>(i)].pt.c);
            for (size_t c = a + 1; c < vps.size(); ++c) {
                const int j = vps[c];
                const int v = leaf.localIdx(vp_data[static_cast<size_t>(j)].pt.r, vp_data[static_cast<size_t>(j)].pt.c);
                const dist_t d = leaf.D[static_cast<size_t>(u * (b * b) + v)];
                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                const size_t idx_ji = static_cast<size_t>(j) * k + static_cast<size_t>(i);
                matrix[idx_ij] = d;
                matrix[idx_ji] = d;
                if (d < INF && d == static_cast<dist_t>(vp_data[static_cast<size_t>(i)].pt.manhattan(vp_data[static_cast<size_t>(j)].pt))) {
                    scratch.closed[idx_ij] = 1;
                    scratch.closed[idx_ji] = 1;
                }
            }
        }
    }

    // 5. Bottom-up batched lifting + projected multi-level meet:
    // Ascend the hierarchy level by level from lvl = 1 to H.
    int cur_nx = num_leaf_tiles_x;
    int cur_ny = num_leaf_tiles_y;

    for (int lvl = 1; lvl <= H; ++lvl) {
        const int next_nx = cur_nx / gw;
        const int next_ny = cur_ny / gh;
        const int next_num_nodes = next_nx * next_ny;

        const auto* input_vps = (lvl == 1) ? &scratch.leaf_vps : &scratch.node_vps[(lvl - 2) % 2];
        auto* output_vps = &scratch.node_vps[(lvl - 1) % 2];

        for (int p_idx = 0; p_idx < next_num_nodes; ++p_idx) {
            (*output_vps)[static_cast<size_t>(p_idx)].clear();
        }

        for (int py = 0; py < next_ny; ++py) {
            for (int px = 0; px < next_nx; ++px) {
                const int p_idx = py * next_nx + px;
                const InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];
                const int nExt = static_cast<int>(parent.extPorts.size());

                struct ChildInfo {
                    int quad;
                    int offset;
                    int nPorts;
                    const std::vector<int>* vps;
                };
                std::array<ChildInfo, kMaxActiveChildren> active_children;
                size_t num_active_children = 0;

                for (int cy = 0; cy < gh; ++cy) {
                    for (int cx = 0; cx < gw; ++cx) {
                        const int c_idx = (py * gh + cy) * cur_nx + (px * gw + cx);
                        const auto& c_vps = (*input_vps)[static_cast<size_t>(c_idx)];
                        if (c_vps.empty()) continue;

                        const int quad = cy * gw + cx;
                        const int s_offset = parent.childOffset[static_cast<size_t>(quad)];
                        const int nPorts = static_cast<int>(vp_data[static_cast<size_t>(c_vps[0])].X[static_cast<size_t>(lvl - 1)].size());

                        active_children[num_active_children++] = {quad, s_offset, nPorts, &c_vps};
                    }
                }

                if (num_active_children == 0) continue;

                // Merge viewpoints into parent's viewpoint list for the next hierarchy level
                auto& p_vps = (*output_vps)[static_cast<size_t>(p_idx)];
                for (size_t ac = 0; ac < num_active_children; ++ac) {
                    const auto& ch = active_children[ac];
                    for (int vp_idx : *ch.vps) {
                        p_vps.push_back(vp_idx);
                    }
                }

                // ------------------------------------------------------------
                // 5a. Batched lifting to parent's exterior ports (only if lvl < H)
                // Project child-level distance vectors X(lvl-1) to parent exterior ports
                // using 16-lane AVX2 SIMD min-plus outer products:
                //   X(lvl)[e] = min_{p} { X(lvl-1)[p] + child_to_ext[quad](p, e) }
                // ------------------------------------------------------------
                if (lvl < H) {
                    for (size_t ac = 0; ac < num_active_children; ++ac) {
                        const auto& ch = active_children[ac];
                        const auto& submat = parent.child_to_ext[static_cast<size_t>(ch.quad)];

                        for (int vp_idx : *ch.vps) {
                            auto& x_out = vp_data[static_cast<size_t>(vp_idx)].X[static_cast<size_t>(lvl)];
                            x_out.assign(static_cast<size_t>(nExt), INF);
                            const auto& x_in = vp_data[static_cast<size_t>(vp_idx)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_in = vp_data[static_cast<size_t>(vp_idx)].active[static_cast<size_t>(lvl - 1)];

                            for (int p : act_in) {
                                const dist_t xp = x_in[static_cast<size_t>(p)];
                                const dist_t* row = &submat[static_cast<size_t>(p * nExt)];
                                const __m256i v_xp = _mm256_set1_epi16(static_cast<short>(xp));
                                int e = 0;
                                for (; e + 16 <= nExt; e += 16) {
                                    const __m256i v_row = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row[e]));
                                    const __m256i v_sum = _mm256_adds_epu16(v_xp, v_row);
                                    const __m256i v_cur = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&x_out[static_cast<size_t>(e)]));
                                    const __m256i v_min = _mm256_min_epu16(v_cur, v_sum);
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&x_out[static_cast<size_t>(e)]), v_min);
                                }
                                for (; e < nExt; ++e) {
                                    const dist_t de = row[e];
                                    if (de < INF) {
                                        const dist_t cand = add_dist(xp, de);
                                        if (cand < x_out[static_cast<size_t>(e)]) x_out[static_cast<size_t>(e)] = cand;
                                    }
                                }
                            }

                            auto& act_out = vp_data[static_cast<size_t>(vp_idx)].active[static_cast<size_t>(lvl)];
                            act_out.clear();
                            for (int e = 0; e < nExt; ++e) {
                                if (x_out[static_cast<size_t>(e)] < INF) act_out.push_back(e);
                            }
                            std::sort(act_out.begin(), act_out.end(),
                                      [&](int a, int c) { return x_out[static_cast<size_t>(a)] < x_out[static_cast<size_t>(c)]; });
                        }
                    }
                }

                // ------------------------------------------------------------
                // 5b. Inter-child meets: for each pair of children (c1, c2)
                // Shortest paths between viewpoints in distinct sibling children
                // meet at interface ports within parent node W.
                // ------------------------------------------------------------
                const int nV = static_cast<int>(parent.V_U.size());
                for (size_t ac1 = 0; ac1 < num_active_children; ++ac1) {
                    const auto& ch1 = active_children[ac1];

                    for (size_t ac2 = ac1 + 1; ac2 < num_active_children; ++ac2) {
                        const auto& ch2 = active_children[ac2];

                        // O4 Optimization: skip the entire child pair if every cross-pair is already closed
                        bool block_all_closed = true;
                        for (int i : *ch1.vps) {
                            for (int j : *ch2.vps) {
                                if (!scratch.closed[static_cast<size_t>(i) * k + static_cast<size_t>(j)]) {
                                    block_all_closed = false;
                                    break;
                                }
                            }
                            if (!block_all_closed) break;
                        }
                        if (block_all_closed) continue;

                        if (scratch.proj_buf.size() < static_cast<size_t>(ch2.nPorts)) {
                            scratch.proj_buf.resize(static_cast<size_t>(ch2.nPorts));
                        }

                        // O1 Optimization (Row Compaction):
                        // Singleton child groups read parent D_U rows directly.
                        // When child 1 contains multiple viewpoints, compact the union of
                        // active rows into a contiguous buffer (block_buf) once, allowing all
                        // viewpoints in child 1 to reuse the compacted block with optimal cache locality.
                        const bool compact_rows = ch1.vps->size() >= 2;
                        if (compact_rows) {
                            if (scratch.row_marks.size() < static_cast<size_t>(ch1.nPorts)) {
                                scratch.row_marks.resize(static_cast<size_t>(ch1.nPorts));
                            }
                            std::fill_n(scratch.row_marks.data(), ch1.nPorts, -1);
                            scratch.row_list.clear();
                            for (int i : *ch1.vps) {
                                const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                                for (int p1 : act_i) {
                                    if (scratch.row_marks[static_cast<size_t>(p1)] < 0) {
                                        scratch.row_marks[static_cast<size_t>(p1)] = static_cast<int>(scratch.row_list.size());
                                        scratch.row_list.push_back(p1);
                                    }
                                }
                            }
                            const size_t n_rows = scratch.row_list.size();
                            if (scratch.block_buf.size() < n_rows * static_cast<size_t>(ch2.nPorts)) {
                                scratch.block_buf.resize(n_rows * static_cast<size_t>(ch2.nPorts));
                            }
                            for (size_t ri = 0; ri < n_rows; ++ri) {
                                const int p1 = scratch.row_list[ri];
                                const dist_t* src_row = &parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch2.offset)];
                                std::copy_n(src_row, ch2.nPorts, scratch.block_buf.data() + ri * static_cast<size_t>(ch2.nPorts));
                            }
                        }

                        // Project child 1 viewpoints to child 2 boundary ports
                        for (int i : *ch1.vps) {
                            const auto& X_vec = vp_data[static_cast<size_t>(i)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            if (act_i.empty()) continue;

                            // O2 Manhattan check: skip projection if ALL viewpoints in child 2 are already closed
                            bool all_closed = true;
                            for (int j : *ch2.vps) {
                                if (!scratch.closed[static_cast<size_t>(i) * k + static_cast<size_t>(j)]) {
                                    all_closed = false;
                                    break;
                                }
                            }
                            if (all_closed) continue;

                            std::fill_n(scratch.proj_buf.data(), ch2.nPorts, INF);
                            for (int p1 : act_i) {
                                const dist_t xp = X_vec[static_cast<size_t>(p1)];
                                const dist_t* row_ch2 = compact_rows
                                    ? (scratch.block_buf.data()
                                        + static_cast<size_t>(scratch.row_marks[static_cast<size_t>(p1)]) * static_cast<size_t>(ch2.nPorts))
                                    : (&parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch2.offset)]);

                                const __m256i v_xp = _mm256_set1_epi16(static_cast<short>(xp));
                                int p2 = 0;
                                for (; p2 + 16 <= ch2.nPorts; p2 += 16) {
                                    const __m256i v_row = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_ch2[p2]));
                                    const __m256i v_sum = _mm256_adds_epu16(v_xp, v_row);
                                    const __m256i v_cur = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&scratch.proj_buf[static_cast<size_t>(p2)]));
                                    const __m256i v_min = _mm256_min_epu16(v_cur, v_sum);
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&scratch.proj_buf[static_cast<size_t>(p2)]), v_min);
                                }
                                for (; p2 < ch2.nPorts; ++p2) {
                                    const dist_t dij = row_ch2[p2];
                                    if (dij < INF) {
                                        const dist_t cand = add_dist(xp, dij);
                                        if (cand < scratch.proj_buf[static_cast<size_t>(p2)]) scratch.proj_buf[static_cast<size_t>(p2)] = cand;
                                    }
                                }
                            }

                            // Meet projected ch1 distance vector with actual viewpoint distance vectors in child 2:
                            //   matrix[i, j] = min_{p2} { proj_buf[p2] + Y_j[p2] }
                            for (int j : *ch2.vps) {
                                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                                if (scratch.closed[idx_ij]) continue;

                                dist_t best = matrix[idx_ij];
                                const auto& Y_vec = vp_data[static_cast<size_t>(j)].X[static_cast<size_t>(lvl - 1)];
                                const auto& act_j = vp_data[static_cast<size_t>(j)].active[static_cast<size_t>(lvl - 1)];
                                if (act_j.empty()) continue;

                                if (act_j.size() >= 32) {
                                    // --------------------------------------------------------
                                    // Dense AVX2 Min-Plus Dot Product:
                                    // When the destination viewpoint has >= 32 active ports,
                                    // vectorize the min-plus dot product across all 16-lane blocks:
                                    //   sum = adds_epu16(proj_buf, Y_vec)
                                    //   acc = min_epu16(acc, sum)
                                    // --------------------------------------------------------
                                    __m256i v_acc = _mm256_set1_epi16(static_cast<short>(INF));
                                    const dist_t* pb = scratch.proj_buf.data();
                                    const dist_t* yv = Y_vec.data();
                                    int p2 = 0;
                                    for (; p2 + 16 <= ch2.nPorts; p2 += 16) {
                                        const __m256i v_p = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pb + p2));
                                        const __m256i v_y = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(yv + p2));
                                        const __m256i v_s = _mm256_adds_epu16(v_p, v_y);
                                        v_acc = _mm256_min_epu16(v_acc, v_s);
                                    }
                                    // Horizontal reduction across 16 lanes of __m256i:
                                    // 1. Cast and extract two 128-bit halves and take pairwise min
                                    // 2. Perform 8-byte, 4-byte, and 2-byte shift-and-min reductions
                                    const __m128i v_lo = _mm256_castsi256_si128(v_acc);
                                    const __m128i v_hi = _mm256_extracti128_si256(v_acc, 1);
                                    __m128i v_m = _mm_min_epu16(v_lo, v_hi);
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 8));
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 4));
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 2));
                                    dist_t cand = static_cast<dist_t>(_mm_cvtsi128_si32(v_m) & 0xFFFF);
                                    // Scalar tail for remaining boundary ports (< 16)
                                    for (; p2 < ch2.nPorts; ++p2) {
                                        const dist_t tc = add_dist(pb[static_cast<size_t>(p2)], yv[static_cast<size_t>(p2)]);
                                        if (tc < cand) cand = tc;
                                    }
                                    if (cand < best) best = cand;
                                } else {
                                    // Sparse sorted active list iteration with early break
                                    for (int p2 : act_j) {
                                        const dist_t yj = Y_vec[static_cast<size_t>(p2)];
                                        if (yj >= best) break;
                                        const dist_t cand = add_dist(scratch.proj_buf[static_cast<size_t>(p2)], yj);
                                        if (cand < best) best = cand;
                                    }
                                }

                                if (best < matrix[idx_ij]) {
                                    matrix[idx_ij] = best;
                                    matrix[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = best;
                                }
                                // If the current shortest path matches the Manhattan lower bound,
                                // no subsequent ancestor level can find a shorter path. Mark closed.
                                if (best < INF && best == static_cast<dist_t>(vp_data[static_cast<size_t>(i)].pt.manhattan(vp_data[static_cast<size_t>(j)].pt))) {
                                    scratch.closed[idx_ij] = 1;
                                    scratch.closed[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = 1;
                                }
                            }
                        }
                    }

                    // --------------------------------------------------------
                    // 5c. Intra-child shortcut through parent:
                    // For pairs of viewpoints residing in the SAME child, check whether
                    // a path leaving the child, traversing exterior seams in parent W,
                    // and re-entering the child is strictly shorter than the intra-child path.
                    // --------------------------------------------------------
                    if (ch1.vps->size() >= 2) {
                        // Compact the symmetric child block over the union of active rows
                        if (scratch.row_marks.size() < static_cast<size_t>(ch1.nPorts)) {
                            scratch.row_marks.resize(static_cast<size_t>(ch1.nPorts));
                        }
                        std::fill_n(scratch.row_marks.data(), ch1.nPorts, -1);
                        scratch.row_list.clear();
                        for (int i : *ch1.vps) {
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            for (int p1 : act_i) {
                                if (scratch.row_marks[static_cast<size_t>(p1)] < 0) {
                                    scratch.row_marks[static_cast<size_t>(p1)] = static_cast<int>(scratch.row_list.size());
                                    scratch.row_list.push_back(p1);
                                }
                            }
                        }
                        const size_t n_rows = scratch.row_list.size();
                        if (scratch.block_buf.size() < n_rows * n_rows) {
                            scratch.block_buf.resize(n_rows * n_rows);
                        }
                        for (size_t ri = 0; ri < n_rows; ++ri) {
                            const int p1 = scratch.row_list[ri];
                            const dist_t* src_row = &parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch1.offset)];
                            dist_t* dst_row = scratch.block_buf.data() + ri * n_rows;
                            for (size_t cj = 0; cj < n_rows; ++cj) {
                                dst_row[cj] = src_row[scratch.row_list[cj]];
                            }
                        }

                        for (size_t a = 0; a < ch1.vps->size(); ++a) {
                            const int i = (*ch1.vps)[a];
                            const auto& X_vec = vp_data[static_cast<size_t>(i)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            if (act_i.empty()) continue;

                            for (size_t c = a + 1; c < ch1.vps->size(); ++c) {
                                const int j = (*ch1.vps)[c];
                                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                                if (scratch.closed[idx_ij]) continue;

                                dist_t best = matrix[idx_ij];
                                const auto& Y_vec = vp_data[static_cast<size_t>(j)].X[static_cast<size_t>(lvl - 1)];
                                const auto& act_j = vp_data[static_cast<size_t>(j)].active[static_cast<size_t>(lvl - 1)];
                                if (act_j.empty()) continue;

                                if (add_dist(X_vec[static_cast<size_t>(act_i[0])],
                                             Y_vec[static_cast<size_t>(act_j[0])]) >= best) {
                                    continue;
                                }

                                for (int p_s : act_i) {
                                    const dist_t xi = X_vec[static_cast<size_t>(p_s)];
                                    if (xi >= best) break;
                                    const dist_t* row_s = scratch.block_buf.data()
                                        + static_cast<size_t>(scratch.row_marks[static_cast<size_t>(p_s)]) * n_rows;

                                    for (int p_t : act_j) {
                                        const dist_t yj = Y_vec[static_cast<size_t>(p_t)];
                                        if (add_dist(xi, yj) >= best) break;
                                        const dist_t dij = row_s[scratch.row_marks[static_cast<size_t>(p_t)]];
                                        if (dij < INF) {
                                            const dist_t cand = add_dist(xi, dij, yj);
                                            if (cand < best) best = cand;
                                        }
                                    }
                                }

                                if (best < matrix[idx_ij]) {
                                    matrix[idx_ij] = best;
                                    matrix[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = best;
                                }
                                if (best < INF && best == static_cast<dist_t>(vp_data[static_cast<size_t>(i)].pt.manhattan(vp_data[static_cast<size_t>(j)].pt))) {
                                    scratch.closed[idx_ij] = 1;
                                    scratch.closed[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        cur_nx = next_nx;
        cur_ny = next_ny;
    }
}

std::vector<dist_t> MinPlusOracle::queryAllPairs(
    const std::vector<Point>& viewpoints, const GridMap& grid) const {
    std::vector<dist_t> matrix;
    queryAllPairs(viewpoints, grid, matrix, query_scratch);
    return matrix;
}

// ============================================================================
// Incremental All-Pairs Matrix Query: Cached Multi-Level Reuse
// ============================================================================

uint32_t MinPlusOracle::nodeVersionAtLvl(int lvl, int leaf_tx, int leaf_ty) const noexcept {
    if (lvl < 1 || lvl > H) return 0;
    if (leaf_tx < 0 || leaf_tx >= num_leaf_tiles_x || leaf_ty < 0 || leaf_ty >= num_leaf_tiles_y) return 0;
    int tx = leaf_tx;
    int ty = leaf_ty;
    for (int l = 0; l < lvl; ++l) {
        tx /= config.group_w;
        ty /= config.group_h;
    }
    return levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(tileIndexAtLevel(lvl, tx, ty))].version;
}

void MinPlusOracle::queryAllPairsIncremental(
    const std::vector<Point>& viewpoints, const GridMap& grid,
    std::vector<dist_t>& out_matrix, HTripQueryScratch& scratch,
    HTripIncrementalState& state) const {

    const size_t k = viewpoints.size();
    if (k == 0) {
        out_matrix.clear();
        state.initialized = false;
        return;
    }

    // 1. Hybrid dispatch:
    // If a hybrid threshold is configured and k is sufficiently small, fall back
    // to search baselines. Leave the incremental state uninitialized so stale
    // lifts are never adopted.
    if (config.hybrid_k_threshold > 0 && static_cast<int>(k) <= config.hybrid_k_threshold) {
        if (grid.hasWeights()) {
            out_matrix = DijkstraOracle::queryAllPairs(grid, viewpoints);
        } else {
            out_matrix = BfsOracle::queryAllPairs(grid, viewpoints);
        }
        state.initialized = false;
        return;
    }

    // 2. Consistency guards:
    // Trigger a full recomputation if the state is uninitialized, the hierarchy epoch
    // changed (e.g., dynamic root expansion occurred), or the viewpoint set changed.
    bool full_reset = !state.initialized || state.epoch_seen != epoch || state.vp_pts.size() != k;
    if (!full_reset) {
        for (size_t i = 0; i < k; ++i) {
            if (state.vp_pts[i] != viewpoints[i]) {
                full_reset = true;
                break;
            }
        }
    }

    if (full_reset) {
        // Execute a complete from-scratch all-pairs query
        queryAllPairs(viewpoints, grid, state.matrix, scratch);

        // Take ownership of the freshly computed viewpoint lifts into persistent state cache
        state.vp_cache.assign(scratch.vp_data.begin(), scratch.vp_data.begin() + k);
        state.vp_pts = viewpoints;
        state.epoch_seen = epoch;
        state.initialized = true;
        state.vp_leaf_version.assign(k, 0);
        state.vp_lvl_version.assign(k, std::vector<uint32_t>(static_cast<size_t>(H + 1), 0));
        state.node_version_snap.assign(static_cast<size_t>(H), std::vector<uint32_t>());
        const int b = config.tile_size_b;
        for (size_t i = 0; i < k; ++i) {
            const Point p = viewpoints[i];
            const int tx = p.c / b;
            const int ty = p.r / b;
            if (grid.isPassable(p) && tx >= 0 && tx < num_leaf_tiles_x && ty >= 0 && ty < num_leaf_tiles_y) {
                state.vp_leaf_version[i] = leaves[static_cast<size_t>(leafIndex(tx, ty))].version;
                for (int lvl = 1; lvl <= H; ++lvl) {
                    state.vp_lvl_version[i][static_cast<size_t>(lvl)] = nodeVersionAtLvl(lvl, tx, ty);
                }
            }
        }
        for (int lvl = 1; lvl <= H; ++lvl) {
            const auto& lv = levels[static_cast<size_t>(lvl - 1)];
            state.node_version_snap[static_cast<size_t>(lvl - 1)].resize(lv.size());
            for (size_t n = 0; n < lv.size(); ++n) {
                state.node_version_snap[static_cast<size_t>(lvl - 1)][n] = lv[n].version;
            }
        }
        state.last_update_counter = update_counter;
        out_matrix = state.matrix;
        return;
    }

    // 3. P11 Fast-Path Return:
    // If no hierarchy updates occurred since the previous query, the entire distance
    // matrix is guaranteed identical. Return immediately in O(1) time without visiting
    // any hierarchy nodes or scratch buffers.
    if (state.last_update_counter == update_counter) {
        out_matrix = state.matrix;
        return;
    }

    out_matrix = state.matrix;

    // Self-distances: a viewpoint blocked during reset and reopened must regain diagonal 0
    for (size_t i = 0; i < k; ++i) {
        out_matrix[i * k + i] = grid.isPassable(viewpoints[i]) ? 0 : INF;
    }

    const int b = config.tile_size_b;
    const int gw = config.group_w;
    const int gh = config.group_h;
    const int num_leaves = num_leaf_tiles_x * num_leaf_tiles_y;

    scratch.ensureCapacity(k, H, static_cast<size_t>(num_leaves));
    std::fill_n(scratch.closed.data(), k * k, static_cast<uint8_t>(0));
    for (size_t i = 0; i < k; ++i) {
        scratch.closed[i * k + i] = 1;
    }
    for (int l_idx = 0; l_idx < num_leaves; ++l_idx) {
        scratch.leaf_vps[static_cast<size_t>(l_idx)].clear();
    }

    // Persistent viewpoint lift cache
    if (state.vp_cache.size() < k) {
        state.vp_cache.resize(k);
    }
    auto& vp_data = state.vp_cache;
    std::fill_n(scratch.vp_new.data(), k, static_cast<uint8_t>(0));
    auto& vp_new = scratch.vp_new;

    // 4. Manhattan Lower-Bound Pre-Closing:
    // Pairs already reaching their theoretical Manhattan distance lower bound
    // |r1-r2| + |c1-c2| under monotone edge-weight decrease can never be improved
    // by any upper-level seam meet. Pre-close them to completely prune upper-level meets.
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            const size_t idx_ij = i * k + j;
            const dist_t m = out_matrix[idx_ij];
            if (m < INF && m == static_cast<dist_t>(viewpoints[i].manhattan(viewpoints[j]))) {
                scratch.closed[idx_ij] = 1;
                scratch.closed[j * k + i] = 1;
            }
        }
    }

    // 5. Level 0: Leaf grouping with version-validated X[0] reuse
    // For each viewpoint, check whether its base leaf tile version has changed since
    // the last query. If the leaf version is unchanged, the precomputed X[0] distance
    // vector to perimeter ports remains strictly exact and is reused without computation.
    for (size_t i = 0; i < k; ++i) {
        const Point p = viewpoints[i];
        vp_data[i].pt = p;
        if (!grid.isPassable(p)) continue;
        const int tx = p.c / b;
        const int ty = p.r / b;
        if (tx < 0 || tx >= num_leaf_tiles_x || ty < 0 || ty >= num_leaf_tiles_y) continue;
        const int l_idx = leafIndex(tx, ty);
        scratch.leaf_vps[static_cast<size_t>(l_idx)].push_back(static_cast<int>(i));

        const LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        if (leaf.version == state.vp_leaf_version[i]) continue; // Cached X[0] still exact!

        const int u = leaf.localIdx(p.r, p.c);
        const int nP = static_cast<int>(leaf.ports.size());
        vp_data[i].X[0].assign(static_cast<size_t>(nP), INF);
        vp_data[i].active[0].clear();
        for (int port = 0; port < nP; ++port) {
            const int pCell = leaf.localIdx(leaf.ports[static_cast<size_t>(port)].r,
                                            leaf.ports[static_cast<size_t>(port)].c);
            const dist_t d = leaf.D[static_cast<size_t>(u * (b * b) + pCell)];
            vp_data[i].X[0][static_cast<size_t>(port)] = d;
            if (d < INF) vp_data[i].active[0].push_back(port);
        }
        std::sort(vp_data[i].active[0].begin(), vp_data[i].active[0].end(),
                  [&](int a, int c) { return vp_data[i].X[0][static_cast<size_t>(a)] < vp_data[i].X[0][static_cast<size_t>(c)]; });
        state.vp_leaf_version[i] = leaf.version;
        vp_new[i] = 1; // Mark viewpoint as having updated port distances
    }

    // 6. Intra-leaf initial pairwise evaluations (changed leaves only; monotone min)
    // Only re-evaluate pairs if at least one viewpoint or leaf tile was modified.
    for (int l_idx = 0; l_idx < num_leaves; ++l_idx) {
        const auto& vps = scratch.leaf_vps[static_cast<size_t>(l_idx)];
        if (vps.size() < 2) continue;
        const LeafTile& leaf = leaves[static_cast<size_t>(l_idx)];
        for (size_t a = 0; a < vps.size(); ++a) {
            const int i = vps[a];
            for (size_t c = a + 1; c < vps.size(); ++c) {
                const int j = vps[c];
                if (!vp_new[static_cast<size_t>(i)] && !vp_new[static_cast<size_t>(j)]) continue;
                const int u = leaf.localIdx(vp_data[static_cast<size_t>(i)].pt.r, vp_data[static_cast<size_t>(i)].pt.c);
                const int v = leaf.localIdx(vp_data[static_cast<size_t>(j)].pt.r, vp_data[static_cast<size_t>(j)].pt.c);
                const dist_t d = leaf.D[static_cast<size_t>(u * (b * b) + v)];
                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                const size_t idx_ji = static_cast<size_t>(j) * k + static_cast<size_t>(i);
                if (d < out_matrix[idx_ij]) out_matrix[idx_ij] = d;
                if (d < out_matrix[idx_ji]) out_matrix[idx_ji] = d;
            }
        }
    }

    // 7. Bottom-up batched lifting + projected meet, restricted to changed nodes
    int cur_nx = num_leaf_tiles_x;
    int cur_ny = num_leaf_tiles_y;

    for (int lvl = 1; lvl <= H; ++lvl) {
        const int next_nx = cur_nx / gw;
        const int next_ny = cur_ny / gh;
        const int next_num_nodes = next_nx * next_ny;

        const auto* input_vps = (lvl == 1) ? &scratch.leaf_vps : &scratch.node_vps[(lvl - 2) % 2];
        auto* output_vps = &scratch.node_vps[(lvl - 1) % 2];
        for (int p_idx = 0; p_idx < next_num_nodes; ++p_idx) {
            (*output_vps)[static_cast<size_t>(p_idx)].clear();
        }

        for (int py = 0; py < next_ny; ++py) {
            for (int px = 0; px < next_nx; ++px) {
                const int p_idx = py * next_nx + px;
                const InternalNode& parent = levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];
                const int nExt = static_cast<int>(parent.extPorts.size());

                struct ChildInfo {
                    int quad;
                    int offset;
                    int nPorts;
                    const std::vector<int>* vps;
                };
                std::array<ChildInfo, kMaxActiveChildren> active_children;
                size_t num_active_children = 0;
                bool node_has_new = false;

                for (int cy = 0; cy < gh; ++cy) {
                    for (int cx = 0; cx < gw; ++cx) {
                        const int c_idx = (py * gh + cy) * cur_nx + (px * gw + cx);
                        const auto& c_vps = (*input_vps)[static_cast<size_t>(c_idx)];
                        if (c_vps.empty()) continue;
                        const int quad = cy * gw + cx;
                        const int s_offset = parent.childOffset[static_cast<size_t>(quad)];
                        const int nPorts = static_cast<int>(vp_data[static_cast<size_t>(c_vps[0])].X[static_cast<size_t>(lvl - 1)].size());
                        active_children[num_active_children++] = {quad, s_offset, nPorts, &c_vps};
                        for (int vp_idx : c_vps) {
                            if (vp_new[static_cast<size_t>(vp_idx)]) node_has_new = true;
                        }
                    }
                }
                if (num_active_children == 0) continue;

                auto& p_vps = (*output_vps)[static_cast<size_t>(p_idx)];
                for (size_t ac = 0; ac < num_active_children; ++ac) {
                    for (int vp_idx : *active_children[ac].vps) {
                        p_vps.push_back(vp_idx);
                    }
                }

                const uint32_t node_ver = parent.version;
                const bool node_changed = node_ver != state.node_version_snap[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)];
                const bool node_active = node_changed || node_has_new;
                if (node_active) {
                    state.node_version_snap[static_cast<size_t>(lvl - 1)][static_cast<size_t>(p_idx)] = node_ver;
                }

                // 1. Batched lifting to parent's exterior ports (only if lvl < H)
                if (lvl < H) {
                    for (size_t ac = 0; ac < num_active_children; ++ac) {
                        const auto& ch = active_children[ac];
                        const auto& submat = parent.child_to_ext[static_cast<size_t>(ch.quad)];
                        for (int vp_idx : *ch.vps) {
                            const bool need_lift = vp_new[static_cast<size_t>(vp_idx)]
                                || (node_ver != state.vp_lvl_version[static_cast<size_t>(vp_idx)][static_cast<size_t>(lvl)]);
                            if (!need_lift) {
                                vp_new[static_cast<size_t>(vp_idx)] = 0;
                                continue;
                            }

                            auto& x_out = vp_data[static_cast<size_t>(vp_idx)].X[static_cast<size_t>(lvl)];
                            x_out.assign(static_cast<size_t>(nExt), INF);
                            const auto& x_in = vp_data[static_cast<size_t>(vp_idx)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_in = vp_data[static_cast<size_t>(vp_idx)].active[static_cast<size_t>(lvl - 1)];

                            for (int p : act_in) {
                                dist_t xp = x_in[static_cast<size_t>(p)];
                                const dist_t* row = &submat[static_cast<size_t>(p * nExt)];
                                __m256i v_xp = _mm256_set1_epi16(static_cast<short>(xp));
                                int e = 0;
                                for (; e + 16 <= nExt; e += 16) {
                                    __m256i v_row = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row[e]));
                                    __m256i v_sum = _mm256_adds_epu16(v_xp, v_row);
                                    __m256i v_cur = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&x_out[static_cast<size_t>(e)]));
                                    __m256i v_min = _mm256_min_epu16(v_cur, v_sum);
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&x_out[static_cast<size_t>(e)]), v_min);
                                }
                                for (; e < nExt; ++e) {
                                    dist_t de = row[e];
                                    if (de < INF) {
                                        dist_t cand = add_dist(xp, de);
                                        if (cand < x_out[static_cast<size_t>(e)]) x_out[static_cast<size_t>(e)] = cand;
                                    }
                                }
                            }

                            auto& act_out = vp_data[static_cast<size_t>(vp_idx)].active[static_cast<size_t>(lvl)];
                            act_out.clear();
                            for (int e = 0; e < nExt; ++e) {
                                if (x_out[static_cast<size_t>(e)] < INF) act_out.push_back(e);
                            }
                            std::sort(act_out.begin(), act_out.end(),
                                      [&](int a, int c) { return x_out[static_cast<size_t>(a)] < x_out[static_cast<size_t>(c)]; });

                            state.vp_lvl_version[static_cast<size_t>(vp_idx)][static_cast<size_t>(lvl)] = node_ver;
                            vp_new[static_cast<size_t>(vp_idx)] = 1;
                        }
                    }
                }

                if (!node_active) continue;

                // 2. Inter-child meets: for each pair of children (c1, c2)
                const int nV = static_cast<int>(parent.V_U.size());
                for (size_t ac1 = 0; ac1 < num_active_children; ++ac1) {
                    const auto& ch1 = active_children[ac1];

                    for (size_t ac2 = ac1 + 1; ac2 < num_active_children; ++ac2) {
                        const auto& ch2 = active_children[ac2];

                        bool block_all_closed = true;
                        for (int i : *ch1.vps) {
                            for (int j : *ch2.vps) {
                                if (!scratch.closed[static_cast<size_t>(i) * k + static_cast<size_t>(j)]) {
                                    block_all_closed = false;
                                    break;
                                }
                            }
                            if (!block_all_closed) break;
                        }
                        if (block_all_closed) continue;

                        // Nothing to do unless this node's D_U changed or some
                        // participating viewpoint was re-lifted this call.
                        if (!node_changed) {
                            bool any_new = false;
                            for (int v : *ch1.vps) {
                                if (vp_new[static_cast<size_t>(v)]) {
                                    any_new = true;
                                    break;
                                }
                            }
                            if (!any_new) {
                                for (int v : *ch2.vps) {
                                    if (vp_new[static_cast<size_t>(v)]) {
                                        any_new = true;
                                        break;
                                    }
                                }
                            }
                            if (!any_new) continue;
                        }

                        if (scratch.proj_buf.size() < static_cast<size_t>(ch2.nPorts)) {
                            scratch.proj_buf.resize(static_cast<size_t>(ch2.nPorts));
                        }
                        // O1: singleton child groups read D_U rows directly; wider
                        // groups compact the union of active ch1 rows once and share
                        // them across all viewpoints in child 1.
                        const bool compact_rows = ch1.vps->size() >= 2;
                        if (compact_rows) {
                            if (scratch.row_marks.size() < static_cast<size_t>(ch1.nPorts)) {
                                scratch.row_marks.resize(static_cast<size_t>(ch1.nPorts));
                            }
                            std::fill_n(scratch.row_marks.data(), ch1.nPorts, -1);
                            scratch.row_list.clear();
                            for (int i : *ch1.vps) {
                                const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                                for (int p1 : act_i) {
                                    if (scratch.row_marks[static_cast<size_t>(p1)] < 0) {
                                        scratch.row_marks[static_cast<size_t>(p1)] = static_cast<int>(scratch.row_list.size());
                                        scratch.row_list.push_back(p1);
                                    }
                                }
                            }
                            const size_t n_rows = scratch.row_list.size();
                            if (scratch.block_buf.size() < n_rows * static_cast<size_t>(ch2.nPorts)) {
                                scratch.block_buf.resize(n_rows * static_cast<size_t>(ch2.nPorts));
                            }
                            for (size_t ri = 0; ri < n_rows; ++ri) {
                                const int p1 = scratch.row_list[ri];
                                const dist_t* src_row = &parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch2.offset)];
                                std::copy_n(src_row, ch2.nPorts, scratch.block_buf.data() + ri * static_cast<size_t>(ch2.nPorts));
                            }
                        }

                        for (int i : *ch1.vps) {
                            const auto& X_vec = vp_data[static_cast<size_t>(i)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            if (act_i.empty()) continue;

                            bool all_closed = true;
                            for (int j : *ch2.vps) {
                                if (!scratch.closed[static_cast<size_t>(i) * k + static_cast<size_t>(j)]) {
                                    all_closed = false;
                                    break;
                                }
                            }
                            if (all_closed) continue;

                            std::fill_n(scratch.proj_buf.data(), ch2.nPorts, INF);
                            for (int p1 : act_i) {
                                dist_t xp = X_vec[static_cast<size_t>(p1)];
                                const dist_t* row_ch2 = compact_rows
                                    ? (scratch.block_buf.data()
                                        + static_cast<size_t>(scratch.row_marks[static_cast<size_t>(p1)]) * static_cast<size_t>(ch2.nPorts))
                                    : (&parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch2.offset)]);

                                __m256i v_xp = _mm256_set1_epi16(static_cast<short>(xp));
                                int p2 = 0;
                                for (; p2 + 16 <= ch2.nPorts; p2 += 16) {
                                    __m256i v_row = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&row_ch2[p2]));
                                    __m256i v_sum = _mm256_adds_epu16(v_xp, v_row);
                                    __m256i v_cur = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&scratch.proj_buf[static_cast<size_t>(p2)]));
                                    __m256i v_min = _mm256_min_epu16(v_cur, v_sum);
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&scratch.proj_buf[static_cast<size_t>(p2)]), v_min);
                                }
                                for (; p2 < ch2.nPorts; ++p2) {
                                    dist_t dij = row_ch2[p2];
                                    if (dij < INF) {
                                        dist_t cand = add_dist(xp, dij);
                                        if (cand < scratch.proj_buf[static_cast<size_t>(p2)]) scratch.proj_buf[static_cast<size_t>(p2)] = cand;
                                    }
                                }
                            }

                            for (int j : *ch2.vps) {
                                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                                if (scratch.closed[idx_ij]) continue;
                                if (!node_changed && !vp_new[static_cast<size_t>(i)] && !vp_new[static_cast<size_t>(j)]) continue;

                                dist_t best = out_matrix[idx_ij];
                                const auto& Y_vec = vp_data[static_cast<size_t>(j)].X[static_cast<size_t>(lvl - 1)];
                                const auto& act_j = vp_data[static_cast<size_t>(j)].active[static_cast<size_t>(lvl - 1)];
                                if (act_j.empty()) continue;

                                if (act_j.size() >= 32) {
                                    __m256i v_acc = _mm256_set1_epi16(static_cast<short>(INF));
                                    const dist_t* pb = scratch.proj_buf.data();
                                    const dist_t* yv = Y_vec.data();
                                    int p2 = 0;
                                    for (; p2 + 16 <= ch2.nPorts; p2 += 16) {
                                        __m256i v_p = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pb + p2));
                                        __m256i v_y = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(yv + p2));
                                        __m256i v_s = _mm256_adds_epu16(v_p, v_y);
                                        v_acc = _mm256_min_epu16(v_acc, v_s);
                                    }
                                    __m128i v_lo = _mm256_castsi256_si128(v_acc);
                                    __m128i v_hi = _mm256_extracti128_si256(v_acc, 1);
                                    __m128i v_m = _mm_min_epu16(v_lo, v_hi);
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 8));
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 4));
                                    v_m = _mm_min_epu16(v_m, _mm_srli_si128(v_m, 2));
                                    dist_t cand = static_cast<dist_t>(_mm_cvtsi128_si32(v_m) & 0xFFFF);
                                    // Scalar tail for ports not covered by full 16-lane blocks
                                    for (; p2 < ch2.nPorts; ++p2) {
                                        dist_t tc = add_dist(pb[static_cast<size_t>(p2)], yv[static_cast<size_t>(p2)]);
                                        if (tc < cand) cand = tc;
                                    }
                                    if (cand < best) best = cand;
                                } else {
                                    for (int p2 : act_j) {
                                        dist_t yj = Y_vec[static_cast<size_t>(p2)];
                                        if (yj >= best) break;
                                        dist_t cand = add_dist(scratch.proj_buf[static_cast<size_t>(p2)], yj);
                                        if (cand < best) best = cand;
                                    }
                                }

                                if (best < out_matrix[idx_ij]) {
                                    out_matrix[idx_ij] = best;
                                    out_matrix[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = best;
                                }
                                if (best < INF && best == static_cast<dist_t>(vp_data[static_cast<size_t>(i)].pt.manhattan(vp_data[static_cast<size_t>(j)].pt))) {
                                    scratch.closed[idx_ij] = 1;
                                    scratch.closed[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = 1;
                                }
                            }
                        }
                    }

                    // 3. Intra-child shortcut through parent (for pairs in the SAME child)
                    if (ch1.vps->size() >= 2) {
                        if (scratch.row_marks.size() < static_cast<size_t>(ch1.nPorts)) {
                            scratch.row_marks.resize(static_cast<size_t>(ch1.nPorts));
                        }
                        std::fill_n(scratch.row_marks.data(), ch1.nPorts, -1);
                        scratch.row_list.clear();
                        for (int i : *ch1.vps) {
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            for (int p1 : act_i) {
                                if (scratch.row_marks[static_cast<size_t>(p1)] < 0) {
                                    scratch.row_marks[static_cast<size_t>(p1)] = static_cast<int>(scratch.row_list.size());
                                    scratch.row_list.push_back(p1);
                                }
                            }
                        }
                        const size_t n_rows = scratch.row_list.size();
                        if (scratch.block_buf.size() < n_rows * n_rows) {
                            scratch.block_buf.resize(n_rows * n_rows);
                        }
                        for (size_t ri = 0; ri < n_rows; ++ri) {
                            const int p1 = scratch.row_list[ri];
                            const dist_t* src_row = &parent.D_U[static_cast<size_t>((ch1.offset + p1) * nV + ch1.offset)];
                            dist_t* dst_row = scratch.block_buf.data() + ri * n_rows;
                            for (size_t cj = 0; cj < n_rows; ++cj) {
                                dst_row[cj] = src_row[scratch.row_list[cj]];
                            }
                        }

                        for (size_t a = 0; a < ch1.vps->size(); ++a) {
                            const int i = (*ch1.vps)[a];
                            const auto& X_vec = vp_data[static_cast<size_t>(i)].X[static_cast<size_t>(lvl - 1)];
                            const auto& act_i = vp_data[static_cast<size_t>(i)].active[static_cast<size_t>(lvl - 1)];
                            if (act_i.empty()) continue;

                            for (size_t c = a + 1; c < ch1.vps->size(); ++c) {
                                const int j = (*ch1.vps)[c];
                                const size_t idx_ij = static_cast<size_t>(i) * k + static_cast<size_t>(j);
                                if (scratch.closed[idx_ij]) continue;
                                if (!node_changed && !vp_new[static_cast<size_t>(i)] && !vp_new[static_cast<size_t>(j)]) continue;

                                dist_t best = out_matrix[idx_ij];
                                const auto& Y_vec = vp_data[static_cast<size_t>(j)].X[static_cast<size_t>(lvl - 1)];
                                const auto& act_j = vp_data[static_cast<size_t>(j)].active[static_cast<size_t>(lvl - 1)];
                                if (act_j.empty()) continue;

                                if (add_dist(X_vec[static_cast<size_t>(act_i[0])],
                                             Y_vec[static_cast<size_t>(act_j[0])]) >= best) {
                                    continue;
                                }

                                for (int p_s : act_i) {
                                    dist_t xi = X_vec[static_cast<size_t>(p_s)];
                                    if (xi >= best) break;
                                    const dist_t* row_s = scratch.block_buf.data()
                                        + static_cast<size_t>(scratch.row_marks[static_cast<size_t>(p_s)]) * n_rows;

                                    for (int p_t : act_j) {
                                        dist_t yj = Y_vec[static_cast<size_t>(p_t)];
                                        if (add_dist(xi, yj) >= best) break;
                                        dist_t dij = row_s[scratch.row_marks[static_cast<size_t>(p_t)]];
                                        if (dij < INF) {
                                            dist_t cand = add_dist(xi, dij, yj);
                                            if (cand < best) best = cand;
                                        }
                                    }
                                }

                                if (best < out_matrix[idx_ij]) {
                                    out_matrix[idx_ij] = best;
                                    out_matrix[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = best;
                                }
                                if (best < INF && best == static_cast<dist_t>(vp_data[static_cast<size_t>(i)].pt.manhattan(vp_data[static_cast<size_t>(j)].pt))) {
                                    scratch.closed[idx_ij] = 1;
                                    scratch.closed[static_cast<size_t>(j) * k + static_cast<size_t>(i)] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        cur_nx = next_nx;
        cur_ny = next_ny;
    }

    state.matrix = out_matrix;
    state.last_update_counter = update_counter;
}

// ============================================================================
// Dynamic Hierarchy Expansion: Root Lifting & Memory Footprint Accounting
// ============================================================================

void MinPlusOracle::liftRoot(const GridMap& grid) {
    const int b = config.tile_size_b;
    const int old_S = S;
    const int old_H = H;

    // 1. Scale bounding box side length S and increment hierarchy height H:
    // When exploration expands beyond the current bounding box S x S, H-TRIP doubles
    // or scales the spatial envelope by group_w, introducing a new root level.
    S = old_S * config.group_w;
    H = old_H + 1;

    const int new_leaf_x = S / b;
    const int new_leaf_y = S / b;

    // 2. Expand leaves array:
    // Existing leaf tiles are preserved via std::move, while newly enclosed spatial
    // regions are instantiated and initialized from the grid map.
    std::vector<LeafTile> new_leaves;
    new_leaves.reserve(static_cast<size_t>(new_leaf_x * new_leaf_y));

    for (int ty = 0; ty < new_leaf_y; ++ty) {
        for (int tx = 0; tx < new_leaf_x; ++tx) {
            if (tx < num_leaf_tiles_x && ty < num_leaf_tiles_y) {
                new_leaves.push_back(std::move(leaves[static_cast<size_t>(ty * num_leaf_tiles_x + tx)]));
            } else {
                LeafTile leaf(tx, ty, b);
                leaf.recomputeFromGrid(grid);
                new_leaves.push_back(std::move(leaf));
            }
        }
    }
    leaves = std::move(new_leaves);
    num_leaf_tiles_x = new_leaf_x;
    num_leaf_tiles_y = new_leaf_y;

    // 3. Expand all intermediate levels 1 .. old_H:
    for (int lvl = 1; lvl <= old_H; ++lvl) {
        int side_lvl = b;
        for (int l = 0; l < lvl; ++l) side_lvl *= config.group_w;
        const int old_nx = old_S / side_lvl;
        const int old_ny = old_S / side_lvl;
        const int new_nx = S / side_lvl;
        const int new_ny = S / side_lvl;

        std::vector<InternalNode> new_level;
        new_level.reserve(static_cast<size_t>(new_nx * new_ny));

        for (int ty = 0; ty < new_ny; ++ty) {
            for (int tx = 0; tx < new_nx; ++tx) {
                if (tx < old_nx && ty < old_ny) {
                    new_level.push_back(std::move(levels[static_cast<size_t>(lvl - 1)][static_cast<size_t>(ty * old_nx + tx)]));
                } else {
                    InternalNode node(lvl, tx, ty, b, config.group_w, config.group_h);
                    const int num_children = config.group_w * config.group_h;
                    std::vector<std::span<const dist_t>> child_summaries(static_cast<size_t>(num_children));
                    for (int k = 0; k < num_children; ++k) {
                        const int cy = k / config.group_w;
                        const int cx = k % config.group_w;
                        const int child_tx = tx * config.group_w + cx;
                        const int child_ty = ty * config.group_h + cy;
                        if (lvl == 1) {
                            const int l_idx = child_ty * num_leaf_tiles_x + child_tx;
                            child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(l_idx)].S;
                        } else {
                            const int prev_side = side_lvl / config.group_w;
                            const int prev_nx = S / prev_side;
                            const int p_idx = child_ty * prev_nx + child_tx;
                            child_summaries[static_cast<size_t>(k)] =
                                levels[static_cast<size_t>(lvl - 2)][static_cast<size_t>(p_idx)].S_U;
                        }
                    }
                    node.recomputeFromChildren(child_summaries, grid);
                    new_level.push_back(std::move(node));
                }
            }
        }
        levels[static_cast<size_t>(lvl - 1)] = std::move(new_level);
    }

    // 4. Allocate and compute the new root internal node at level H:
    InternalNode root(H, 0, 0, b, config.group_w, config.group_h);
    const int num_children = config.group_w * config.group_h;
    std::vector<std::span<const dist_t>> root_child_summaries(static_cast<size_t>(num_children));
    const int prev_side = S / config.group_w;
    const int prev_nx = S / prev_side;
    for (int k = 0; k < num_children; ++k) {
        const int cy = k / config.group_w;
        const int cx = k % config.group_w;
        const int p_idx = cy * prev_nx + cx;
        if (old_H == 0) {
            root_child_summaries[static_cast<size_t>(k)] = leaves[static_cast<size_t>(p_idx)].S;
        } else {
            root_child_summaries[static_cast<size_t>(k)] =
                levels[static_cast<size_t>(old_H - 1)][static_cast<size_t>(p_idx)].S_U;
        }
    }
    root.recomputeFromChildren(root_child_summaries, grid);

    std::vector<InternalNode> root_level;
    root_level.push_back(std::move(root));
    levels.push_back(std::move(root_level));

    // 5. Hierarchy Epoch Invalidation:
    // Increment epoch counter to invalidate all viewpoint lift caches across
    // persistent incremental query states.
    epoch++;
}

size_t MinPlusOracle::totalMemoryBytes() const noexcept {
    // Total memory allocated across all leaf intra-tile metric tables, boundary summaries,
    // internal node seam closure matrices, and rectangular child-to-exterior projection tables.
    size_t total = 0;
    for (const auto& l : leaves) {
        total += l.D.size() * sizeof(dist_t);
        total += l.S.size() * sizeof(dist_t);
    }
    for (const auto& lvl : levels) {
        for (const auto& n : lvl) {
            total += n.D_U.size() * sizeof(dist_t);
            total += n.S_U.size() * sizeof(dist_t);
            for (const auto& sub : n.child_to_ext) {
                total += sub.size() * sizeof(dist_t);
            }
        }
    }
    return total;
}

size_t MinPlusOracle::effectiveMemoryBytes(const GridMap& grid) const noexcept {
    // Effective memory footprint: filters out leaf tiles and internal nodes that
    // lie completely outside the physical grid boundaries or contain zero traversable cells.
    size_t total = 0;
    const int b = config.tile_size_b;
    for (const auto& l : leaves) {
        if (l.tx * b >= grid.cols || l.ty * b >= grid.rows) continue;
        bool has_free = false;
        const int max_r = std::min(grid.rows, (l.ty + 1) * b);
        const int max_c = std::min(grid.cols, (l.tx + 1) * b);
        for (int r = l.ty * b; r < max_r && !has_free; ++r) {
            for (int c = l.tx * b; c < max_c; ++c) {
                if (grid.getObserved({r, c}) == CellState::Free) {
                    has_free = true;
                    break;
                }
            }
        }
        if (has_free) {
            total += l.D.size() * sizeof(dist_t);
            total += l.S.size() * sizeof(dist_t);
        }
    }
    for (const auto& lvl : levels) {
        for (const auto& n : lvl) {
            if (n.C0 >= grid.cols || n.R0 >= grid.rows) continue;
            bool has_free = false;
            const int max_r = std::min(grid.rows, n.R0 + n.s);
            const int max_c = std::min(grid.cols, n.C0 + n.s);
            for (int r = n.R0; r < max_r && !has_free; ++r) {
                for (int c = n.C0; c < max_c; ++c) {
                    if (grid.getObserved({r, c}) == CellState::Free) {
                        has_free = true;
                        break;
                    }
                }
            }
            if (has_free) {
                total += n.D_U.size() * sizeof(dist_t);
                total += n.S_U.size() * sizeof(dist_t);
            }
        }
    }
    return total;
}

} // namespace htrip
