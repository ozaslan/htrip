#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>
#include <span>

namespace htrip {

/**
 * @class InternalNode
 * @brief Multi-level hierarchical aggregation node at level lvl >= 1 in the H-TRIP tree.
 *
 * Mathematical and Structural Role:
 * An internal node aggregates g x g child regions (where g = group_w = group_h) into a
 * composite super-tile covering an s x s grid region:
 * s = b * g^lvl (e.g., for b=16, g=2: lvl 1 has s=32, lvl 2 has s=64, lvl 3 has s=128).
 *
 * Key Component Relations:
 * 1. Interface Vertices (V_U): The union of all perimeter ports of its g x g children.
 *    For child side child_s, each child has 4*child_s - 4 ports. Total |V_U| = g^2 * (4*child_s - 4).
 * 2. Exterior Boundary Ports (extPorts): The subset of V_U that lies on the outer perimeter
 *    of the composite super-tile (exactly 4s - 4 ports).
 * 3. Super-Tile Metric Closure (D_U): Dense APSP matrix of size |V_U| x |V_U|, capturing all-pairs
 *    shortest paths between child boundary ports routed through child interiors and inter-child seams.
 * 4. Exterior Summary (S_U): Projection of D_U onto exterior ports extPorts (size: |extPorts| x |extPorts|).
 * 5. Child-to-Exterior Projection Tables (child_to_ext): Pre-packed, cache-aligned rectangular
 *    submatrices mapping each child's boundary ports directly to exterior boundary ports,
 *    enabling O(p_child * p_ext) hierarchical query lifting during shortest-path queries.
 *
 * Tropical Algebra:
 * Initial composition and dynamic updates operate via the tropical min-plus semiring:
 * - Direct edges are loaded from child summaries S_child and 1-hop inter-child seams.
 * - Metric closure is computed via Conway star relaxation with AVX2 SIMD intrinsics.
 */
class InternalNode {
public:
    int level = 1;       ///< Hierarchy level (1, 2, ..., H)
    int tx = 0;          ///< Node column index at this hierarchy level
    int ty = 0;          ///< Node row index at this hierarchy level
    int R0 = 0;          ///< Global grid row coordinate of top-left cell: ty * s
    int C0 = 0;          ///< Global grid column coordinate of top-left cell: tx * s
    int s = 32;          ///< Composite super-tile side length in cells: b * g^lvl
    int group_w = 2;     ///< Horizontal branching factor (default: 2)
    int group_h = 2;     ///< Vertical branching factor (default: 2)
    uint32_t version = 0; ///< Version counter incremented whenever D_U changes

    std::vector<Point> V_U;          ///< Concatenated child boundary port coordinates
    std::vector<int> childOffset;    ///< Prefix offsets: child k's ports span [childOffset[k], childOffset[k+1])
    std::vector<Point> extPorts;     ///< Outer perimeter boundary ports of this super-tile (size: 4s - 4)
    std::vector<int> extPortIdx;     ///< Indices in V_U corresponding to each exterior port in extPorts

    std::vector<dist_t> D_U;         ///< APSP distance metric over all child interface ports: |V_U| x |V_U|
    std::vector<dist_t> S_U;         ///< Exterior boundary summary matrix: |extPorts| x |extPorts|
    std::vector<std::vector<dist_t>> child_to_ext; ///< Direct child-to-exterior submatrices: [quad][p * nExt + e]
    bool child_to_ext_dirty = false; ///< Dirty flag indicating D_U changed and child_to_ext requires refresh

    std::vector<int> active_i_scratch;        ///< Reusable scratch buffer for active row indices in rank-1 updates
    std::vector<int> active_seam_scratch;     ///< Reusable scratch buffer for dirty seam port indices
    std::vector<uint8_t> in_active_seam_scratch; ///< Fast lookup bitmask for active seam port membership
    std::vector<int> reach_k_scratch;         ///< Reusable scratch buffer for reachable vertices during Conway relaxation

    /**
     * @struct SeamCandidate
     * @brief Statically precomputed touching boundary port pair between adjacent sibling children.
     */
    struct SeamCandidate {
        int u; ///< Index in V_U belonging to first child
        int v; ///< Index in V_U belonging to second child
    };
    std::vector<SeamCandidate> seam_candidates; ///< List of geometric touching port pairs (Manhattan distance = 1)

    InternalNode() = default;

    /**
     * @brief Constructs an internal aggregation node at level lvl.
     *
     * Precomputes child partitions, gathers child ports V_U, identifies exterior boundary ports,
     * builds the static touching seam candidate index, and initializes D_U and S_U to tropical identity.
     *
     * @param lvl Hierarchy level (1 <= lvl <= H).
     * @param tx_ Grid column index at level lvl.
     * @param ty_ Grid row index at level lvl.
     * @param b_ Base leaf tile dimension.
     * @param gw Horizontal branching factor (must equal gh).
     * @param gh Vertical branching factor (must equal gw).
     */
    InternalNode(int lvl, int tx_, int ty_, int b_, int gw = 2, int gh = 2);

    /**
     * @brief Refreshes cached child-to-exterior rectangular distance submatrices from D_U.
     *
     * For each child quadrant q, copies the row slice corresponding to child q's ports
     * projected onto extPortIdx into contiguous row-major memory:
     * child_to_ext[q][p * nExt + e] = D_U[(offset_q + p) * nV + extPortIdx[e]].
     *
     * @complexity
     * Time: O(|V_U| * |extPorts|)
     * Memory: O(|V_U| * |extPorts|)
     */
    void refreshChildToExt();

    /**
     * @brief Returns total number of child interface vertices: |V_U|.
     */
    [[nodiscard]] int numInterfaceVertices() const noexcept {
        return static_cast<int>(V_U.size());
    }

    /**
     * @brief Returns number of exterior boundary ports: |extPorts| = 4s - 4.
     */
    [[nodiscard]] int numExtPorts() const noexcept {
        return static_cast<int>(extPorts.size());
    }

    /**
     * @brief Maps a global grid coordinate p to its corresponding port index in V_U.
     *
     * Determines which child quadrant contains p, computes local port index along that child's
     * perimeter, and adds childOffset[quad]. Returns -1 if p does not lie on a child perimeter.
     *
     * @param p Point coordinate.
     * @return Port index in [0, |V_U|), or -1 if not found.
     *
     * @complexity
     * Time: O(1) arithmetic coordinate calculation.
     * Memory: O(1).
     */
    [[nodiscard]] int findVertexIndex(Point p) const noexcept;

    /**
     * @brief Constructs or recomputes D_U bottom-up from child summaries and inter-child seam edges.
     *
     * Algorithmic Pipeline:
     * 1. Loads intra-child metric edges from child_summaries into diagonal blocks of D_U.
     * 2. Inspects static seam_candidates against the grid to insert passable inter-child seam edges.
     * 3. Executes Conway star Kleene closure over active seam endpoints using AVX2 SIMD min-plus relaxations.
     * 4. Projects D_U onto exterior boundary ports to produce S_U and refreshes child_to_ext tables.
     *
     * @param child_summaries Array of boundary summaries S from child nodes (zero-copy spans).
     * @param grid Global occupancy grid.
     *
     * @complexity
     * Time: O(|seams| * |V_U|^2 / 16) with AVX2 SIMD.
     * Memory: O(|V_U|) scratch space.
     */
    void recomputeFromChildren(const std::vector<std::span<const dist_t>>& child_summaries,
                               const GridMap& grid);

    /**
     * @brief Performs an algebraic Conway rank-1 min-plus update on D_U for edge (u, v) with weight w.
     *
     * Formula: D'_U = D_U (min) [ D_U(:, u) + w + D_U(v, :) ] (min) [ D_U(:, v) + w + D_U(u, :) ].
     *
     * @param u Index in V_U.
     * @param v Index in V_U.
     * @param w Edge weight.
     *
     * @complexity
     * Time: O(|V_active|^2) <= O(|V_U|^2).
     */
    void updateRank1(int u, int v, dist_t w);

    /**
     * @brief Inserts an inter-child seam edge between points p1 and p2 into D_U via rank-1 update.
     *
     * @param p1 First seam endpoint.
     * @param p2 Second seam endpoint.
     * @param w Seam traversal cost (default: UNIT_STEP).
     * @return True iff both endpoints exist in V_U and the edge was inserted.
     */
    bool insertSeamEdge(Point p1, Point p2, dist_t w = UNIT_STEP);

    /**
     * @brief Incremental upward propagation of child distance decreases and dynamic seam edges.
     *
     * Algorithmic Pipeline:
     * 1. Relaxes direct edges from child_decreases and newly opened seam_edges into D_U.
     * 2. Gathers unique altered vertices into pivot list K.
     * 3. Executes batched Conway star relaxation over K using AVX2 SIMD min-plus intrinsics.
     * 4. If not root, refreshes exterior summary S_U and returns decreased exterior port pairs
     *    for propagation to parent levels.
     *
     * @param child_decreases List of decreased port distances emitted by child subregions.
     * @param seam_edges List of newly opened inter-child seam edges enclosed at this level.
     * @param is_root True if this node is the hierarchy root (skips S_U propagation).
     * @return List of decreased exterior port pairs (PortDecrease) propagating upwards.
     *
     * @complexity
     * Time: O(|K| * |V_U|^2 / 16) with AVX2 SIMD.
     * Memory: O(|V_U|) scratch buffers.
     */
    std::vector<PortDecrease> applyConwayPerturbations(
        const std::vector<ChildDecrease>& child_decreases,
        const std::vector<SeamEdge>& seam_edges,
        bool is_root = false);

    /**
     * @brief Refreshes exterior summary S_U from D_U and returns strictly decreased port pairs.
     *
     * @return List of PortDecrease records for exterior boundary ports.
     *
     * @complexity
     * Time: O(|extPorts|^2).
     */
    std::vector<PortDecrease> refreshBoundarySummary();
};

} // namespace htrip

