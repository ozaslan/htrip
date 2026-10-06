#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include "leaf_tile.hpp"
#include "internal_node.hpp"
#include <vector>

namespace htrip {

/**
 * @class MinPlusOracle
 * @brief Master distance oracle implementing H-TRIP (Hierarchical Tropical Rank-One
 *        Incremental Propagation) for dynamic robotic navigation and multi-goal planning.
 *
 * Mathematical Foundations (ICRA 2027 Paper):
 * 1. Semiring Architecture:
 *    Operates over the tropical min-plus semiring (R_>=0 U {INF}, min, +).
 *    Matrix multiplication represents path concatenation; Kleene star represents metric closure.
 * 2. Spatial Hierarchy:
 *    The 2D grid space of dimension S x S is partitioned into a multi-level tree:
 *    - Level 0 (Leaves): Array of (S/b) x (S/b) base leaf tiles of side b x b cells.
 *      Each leaf maintains dense APSP metric D (dimension b^2 x b^2) and boundary summary S (4b - 4).
 *    - Levels 1 to H (Internal Nodes): Successively aggregate g x g child regions (where g = group_w = group_h).
 *      Each internal node maintains interface closure D_U and exterior boundary summary S_U.
 * 3. Dynamic Revelations (Theorem 1 - Monotone Exactness):
 *    When free space is revealed, the oracle avoids global BFS or Floyd-Warshall recomputation.
 *    Instead, it applies exact algebraic Conway rank-1 min-plus updates:
 *    D' = D (min) [ D(:, K) + D(K, :) ],
 *    propagating only strictly decreased boundary summary metrics up the tree branch in milliseconds.
 * 4. Batched Queries (AVX2 SIMD Meet):
 *    For k viewpoints, computes the full k x k shortest-path distance matrix by lifting viewpoints
 *    to their Lowest Common Ancestor (LCA) levels and evaluating the dense meet via explicit
 *    AVX2 vector intrinsics (_mm256_adds_epu16, _mm256_min_epu16).
 * 5. Incremental Query Cache:
 *    Maintains version-stamped lifted viewpoints, selectively re-lifting only viewpoints whose
 *    enclosing tiles were perturbed, achieving high amortized speedup for online multi-robot exploration.
 */
class MinPlusOracle {
public:
    HierarchyConfig config;   ///< Configuration knobs (tile size b, branching factor g, hybrid threshold)
    int S = 0;                ///< Spatial bounding side length spanned by hierarchy root (S x S)
    int H = 0;                ///< Hierarchy depth (number of internal levels, H >= 0)
    int num_leaf_tiles_x = 0; ///< Grid width in leaf tiles: S / b
    int num_leaf_tiles_y = 0; ///< Grid height in leaf tiles: S / b

    std::vector<LeafTile> leaves; ///< 1D flat array of base leaf tiles (size: num_leaf_tiles_x * num_leaf_tiles_y)

    /**
     * @brief 2D array of internal hierarchy levels: levels[0] is level 1, levels[H - 1] is root level H.
     */
    std::vector<std::vector<InternalNode>> levels;

    MinPlusOracle() = default;

    /**
     * @brief Constructs an oracle with specified hierarchy parameters.
     */
    explicit MinPlusOracle(HierarchyConfig cfg) : config(cfg) {}

    /**
     * @brief Builds the complete hierarchical min-plus distance oracle bottom-up from an occupancy grid.
     *
     * Construction Pipeline:
     * 1. Determines bounding power S and hierarchy depth H.
     * 2. Instantiates and precomputes all base LeafTile objects (local APSP metric D and boundary S).
     * 3. Constructs levels 1 to H bottom-up: each internal node gathers child summaries, inserts touching
     *    seams, and computes transitive metric closure using AVX2 SIMD.
     * 4. Increments epoch to invalidate any previous incremental cache.
     *
     * @param grid Global occupancy grid.
     * @return True if construction succeeded; false if parameters were invalid.
     *
     * @complexity
     * Time: O((S/b)^2 * b^4 + sum_{lvl=1}^H N_lvl * |seams| * |V_U|^2 / 16).
     * Memory: Sized at construction; persistent index storage.
     */
    bool build(const GridMap& grid);

    /**
     * @brief Returns total number of internal aggregation levels: H.
     */
    [[nodiscard]] int depth() const noexcept { return H; }

    /**
     * @brief Maps leaf tile 2D coordinate (tx, ty) to flat 1D index in leaves array.
     */
    [[nodiscard]] int leafIndex(int tx, int ty) const noexcept {
        return ty * num_leaf_tiles_x + tx;
    }

    /**
     * @brief Maps 2D tile coordinate (tx, ty) at level lvl to flat index in levels[lvl - 1].
     */
    [[nodiscard]] int tileIndexAtLevel(int lvl, int tx, int ty) const noexcept;

    /**
     * @brief Handles a dynamic cell revelation (cell transitions from Unknown/Obstacle to Free).
     *
     * Applies Conway rank-1 updates to the enclosing leaf tile and propagates resulting port decreases
     * and touching seam edges up the ancestor branch to the root.
     *
     * @param p Coordinate of newly opened cell.
     * @param grid Global occupancy grid.
     *
     * @complexity
     * Time: O(b^4 + H * |seams| * |V_U|^2 / 16).
     */
    void onCellOpened(Point p, const GridMap& grid);

    /**
     * @brief Handles a dynamic cell blockage (cell transitions from Free to Obstacle).
     *
     * Recomputes the affected leaf tile via local BFS and rebuilds ancestor node closures along the branch.
     * Increments epoch to trigger a full recomputation on the next incremental query.
     *
     * @param p Coordinate of blocked cell.
     * @param grid Global occupancy grid.
     *
     * @complexity
     * Time: O(b^4 + H * |V_U|^3 / 16).
     */
    void onCellBlocked(Point p, const GridMap& grid);

    /**
     * @brief Handles a single cell traversal weight modification.
     *
     * @param p Cell coordinate.
     * @param new_weight New traversal cost.
     * @param grid Global occupancy grid.
     */
    void onCellWeightChanged(Point p, dist_t new_weight, GridMap& grid);

    /**
     * @brief Batched weight modification for multiple cells within a single leaf tile (tx, ty).
     *
     * Recomputes the containing leaf tile ONCE and propagates ancestor closures up to the root ONCE.
     *
     * @param tx Tile column index.
     * @param ty Tile row index.
     * @param cell_weights List of (point, weight) pairs.
     * @param grid Global occupancy grid.
     */
    void onTileWeightsChanged(int tx, int ty,
                             const std::vector<std::pair<Point, dist_t>>& cell_weights,
                             GridMap& grid);

    /**
     * @brief Batched update when multiple leaf tiles have been altered.
     *
     * Recomputes affected leaf metrics and propagates decreased boundary summaries upwards
     * using exact Conway rank-1 star updates.
     *
     * @param leaf_indices List of flat leaf tile indices.
     * @param grid Global occupancy grid.
     */
    void onBatchLeavesChanged(const std::vector<int>& leaf_indices, const GridMap& grid);

    /**
     * @brief Dynamic batched cell revelation update across multiple leaf tiles.
     *
     * Partitions newly opened cells by leaf tile, applies vectorized Conway rank-1 updates to
     * each dirty leaf, and propagates boundary decreases upward in topological bottom-up order.
     *
     * @param newly_opened List of newly revealed free cell coordinates.
     * @param grid Global occupancy grid.
     *
     * @complexity
     * Time: O(|dirty_leaves| * b^4 / 16 + |dirty_ancestors| * |V_U|^2 / 16).
     */
    void onBatchCellsOpened(const std::vector<Point>& newly_opened, const GridMap& grid);

    /**
     * @struct LiftedViewpoint
     * @brief Distance projection vectors from a single viewpoint to boundary ports across all hierarchy levels.
     */
    struct LiftedViewpoint {
        std::vector<std::vector<dist_t>> X;     ///< [lvl][port]: Shortest path distance to port at level lvl
        std::vector<std::vector<int>> active;   ///< [lvl]: Indices of ports with finite distance (X[lvl][p] < INF)
        std::vector<int> parent_node_idx;       ///< [lvl]: Index of ancestor node at level lvl
        std::vector<int> child_quad;            ///< [lvl]: Quadrant index of child in parent node
        std::vector<int> child_offset;          ///< [lvl]: Start offset in parent.V_U
    };

    /**
     * @struct HTripQueryScratch
     * @brief Preallocated, reusable query workspace ensuring ZERO heap allocations during queries.
     */
    struct HTripQueryScratch {
        struct ViewpointLift {
            Point pt;
            int leaf_tx = 0, leaf_ty = 0;
            std::vector<std::vector<dist_t>> X;      ///< [lvl][port]
            std::vector<std::vector<int>> active;    ///< [lvl][idx]
        };

        std::vector<ViewpointLift> vp_data;
        std::vector<std::vector<int>> leaf_vps;
        std::vector<std::vector<int>> node_vps[2];
        std::vector<dist_t> proj_buf;
        std::vector<dist_t> block_buf;
        std::vector<int> row_marks;
        std::vector<int> row_list;
        std::vector<uint8_t> closed;
        std::vector<uint8_t> vp_new;

        /**
         * @brief Ensures internal scratch buffers are sized for k viewpoints, depth H, and num_leaves.
         *
         * Sizing Invariant:
         * Only grows buffers if k or H exceeds current capacity; never shrinks or reallocates during queries.
         */
        void ensureCapacity(size_t k, int H, size_t num_leaves) {
            if (vp_data.size() < k) {
                vp_data.resize(k);
            }
            for (size_t i = 0; i < k; ++i) {
                if (vp_data[i].X.size() < static_cast<size_t>(H + 1)) {
                    vp_data[i].X.resize(static_cast<size_t>(H + 1));
                    vp_data[i].active.resize(static_cast<size_t>(H + 1));
                }
            }
            if (leaf_vps.size() < num_leaves) {
                leaf_vps.resize(num_leaves);
            }
            for (int b = 0; b < 2; ++b) {
                if (node_vps[b].size() < num_leaves) {
                    node_vps[b].resize(num_leaves);
                }
            }
            if (closed.size() < k * k) {
                closed.resize(k * k);
            }
            if (vp_new.size() < k) {
                vp_new.resize(k);
            }
        }
    };

    mutable HTripQueryScratch query_scratch; ///< Reusable workspace for batched matrix queries
    mutable LiftedViewpoint s_lift_scratch;  ///< Reusable workspace for source in single-pair queries
    mutable LiftedViewpoint t_lift_scratch;  ///< Reusable workspace for target in single-pair queries

    uint32_t update_counter = 0; ///< Version stamp incremented on every topological update
    uint64_t epoch = 0;          ///< Structural epoch incremented on build, liftRoot, or blockages

    /**
     * @struct LastBatchStats
     * @brief Telemetry counters recorded during the most recent dynamic update batch.
     */
    struct LastBatchStats {
        int n_touched_leaves = 0;
        int n_dirty_ancestors = 0;
        int n_pivots = 0;
    };
    LastBatchStats last_batch_stats{};

    /**
     * @struct HTripIncrementalState
     * @brief Persistent cache of lifted viewpoints and version stamps for exact incremental matrix refresh.
     */
    struct HTripIncrementalState {
        bool initialized = false;
        uint64_t epoch_seen = 0;
        uint32_t last_update_counter = 0; ///< Fast-path check: returns cached matrix if update_counter unchanged
        std::vector<Point> vp_pts;
        std::vector<uint32_t> vp_leaf_version;              ///< [i]: Leaf version at time of lift
        std::vector<std::vector<uint32_t>> vp_lvl_version;  ///< [i][lvl]: Ancestor node version at time of lift
        std::vector<std::vector<uint32_t>> node_version_snap; ///< [lvl-1][node_idx]: Version snapshot
        std::vector<dist_t> matrix;                          ///< Persistent k x k result matrix
        std::vector<HTripQueryScratch::ViewpointLift> vp_cache; ///< Persistent copy of lifted port distances
    };

    /**
     * @brief Exact incremental all-pairs k x k distance matrix refresh.
     *
     * Algorithmic Strategy:
     * - If update_counter == state.last_update_counter: zero work, returns cached matrix immediately.
     * - Under cell revelations (insertions), inspects version stamps: only viewpoints whose enclosing
     *   leaf or ancestor nodes experienced perturbations are re-lifted.
     * - Full recomputation is triggered automatically if epoch changed or viewpoint count changed.
     *
     * @param viewpoints List of k query coordinates.
     * @param grid Global occupancy grid.
     * @param out_matrix Output flat k x k distance matrix.
     * @param scratch Reusable scratch buffers.
     * @param state Persistent incremental state cache.
     */
    void queryAllPairsIncremental(const std::vector<Point>& viewpoints,
                                  const GridMap& grid,
                                  std::vector<dist_t>& out_matrix,
                                  HTripQueryScratch& scratch,
                                  HTripIncrementalState& state) const;

    /**
     * @brief Computes shortest-path distance between source s and target t.
     *
     * Algorithm:
     * 1. Same-leaf fast path: If s and t lie in the same leaf tile, looks up D[s, t] directly in O(1).
     * 2. Hierarchical LCA meet: Lifts s and t to their Lowest Common Ancestor level lca_lvl.
     *    Evaluates min-plus meet: d(s, t) = min_{p, q in V_LCA} (d(s, p) + D_LCA(p, q) + d(q, t)).
     *
     * @param s Source coordinate.
     * @param t Target coordinate.
     * @param grid Global occupancy grid.
     * @return Shortest path distance, or INF if disconnected.
     *
     * @complexity
     * Time: O(1) same leaf; O(p_leaf^2 + H * p_child * p_ext + |V_LCA|^2) across tiles.
     * Memory: Uses internal mutable scratch (zero heap allocation).
     */
    [[nodiscard]] dist_t queryDistance(Point s, Point t, const GridMap& grid) const;

    /**
     * @brief Zero-allocation single-pair query using caller-supplied lift scratch buffers.
     */
    [[nodiscard]] dist_t queryDistance(Point s, Point t, const GridMap& grid,
                                       LiftedViewpoint& X_s, LiftedViewpoint& Y_t) const;

    /**
     * @brief Computes full k x k all-pairs distance matrix using preallocated scratch (zero heap allocations).
     *
     * Vectorized AVX2 Min-Plus Meet:
     * Lifts all k viewpoints through the hierarchy into scratch.vp_data.
     * For each pair of viewpoints (i, j), resolves at their LCA node using AVX2 SIMD min-plus dot products:
     * _mm256_adds_epu16 and _mm256_min_epu16 with zero-branching minimum accumulation.
     *
     * @param viewpoints List of k query points.
     * @param grid Global occupancy grid.
     * @param out_matrix Output k x k distance matrix.
     * @param scratch Preallocated query workspace.
     *
     * @complexity
     * Time: O(k * H * p_child * p_ext + k^2 * |V_LCA| / 16) with AVX2 SIMD.
     * Memory: Zero runtime heap allocations.
     */
    void queryAllPairs(
        const std::vector<Point>& viewpoints,
        const GridMap& grid,
        std::vector<dist_t>& out_matrix,
        HTripQueryScratch& scratch) const;

    /**
     * @brief High-throughput all-pairs distance matrix computation (convenience wrapper).
     *
     * Allocates the result vector and uses internal mutable scratch.
     */
    [[nodiscard]] std::vector<dist_t> queryAllPairs(
        const std::vector<Point>& viewpoints, const GridMap& grid) const;

    /**
     * @brief Lifts a point p from its leaf tile to all ancestor levels in the hierarchy.
     *
     * Level 0: copies distances from p to all leaf perimeter ports from D.
     * Level lvl: projects distances from level lvl - 1 onto exterior boundary ports using child_to_ext.
     * Precomputes list of active (finite) port indices to accelerate sparse meets.
     *
     * @param p Point coordinate.
     * @param lifted Output LiftedViewpoint structure.
     */
    void liftPoint(Point p, LiftedViewpoint& lifted) const;

    /**
     * @brief Determines the lowest common ancestor (LCA) level of two points in the hierarchy.
     *
     * @param p1 First point.
     * @param p2 Second point.
     * @return LCA level in [1, H + 1]. Level 1 means sibling leaf tiles sharing a level-1 parent.
     */
    [[nodiscard]] int findLcaLevel(Point p1, Point p2) const noexcept;

    /**
     * @brief Dynamically expands root bounding box: S -> S * g and depth H -> H + 1.
     *
     * Used when exploration reveals territory beyond the initial bounding box.
     * Wraps existing tree as quadrant 0 of a newly instantiated root level.
     *
     * @param grid Global occupancy grid.
     */
    void liftRoot(const GridMap& grid = GridMap());

    /**
     * @brief Total memory consumed by persistent index tables in bytes.
     */
    [[nodiscard]] size_t totalMemoryBytes() const noexcept;

    /**
     * @brief Effective memory consumed by explored non-empty tiles.
     */
    [[nodiscard]] size_t effectiveMemoryBytes(const GridMap& grid) const noexcept;

private:
    struct TileDecreases {
        int tx = 0;
        int ty = 0;
        std::vector<PortDecrease> decreases;
    };

    void propagateUpward(int leaf_tx, int leaf_ty,
                         const std::vector<PortDecrease>& leafDecreases,
                         const std::vector<SeamEdge>& seamEdges,
                         const GridMap& grid);

    void propagateBatchUpward(std::vector<std::pair<int, int>> dirty_tiles,
                              const GridMap& grid);

    [[nodiscard]] uint32_t nodeVersionAtLvl(int lvl, int leaf_tx, int leaf_ty) const noexcept;
};

/// Type alias for the official method name: H-TRIP
using HTripOracle = MinPlusOracle;

} // namespace htrip

