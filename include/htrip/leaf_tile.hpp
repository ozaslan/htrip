#pragma once

#include "types.hpp"
#include "grid_map.hpp"
#include <vector>

namespace htrip {

/**
 * @brief Generates the canonical perimeter boundary port coordinates for a square tile.
 *
 * Perimeter Layout Convention:
 * Ports are extracted along the 4 boundaries of the tile [R0, R0+s) x [C0, C0+s)
 * in the following strictly deterministic order:
 * 1. Top row: (R0, C0 + c) for c in [0, s) (s ports),
 * 2. Bottom row: (R0 + s - 1, C0 + c) for c in [0, s) (s ports),
 * 3. Left column (excluding corners): (R0 + r, C0) for r in [1, s - 1) (s - 2 ports),
 * 4. Right column (excluding corners): (R0 + r, C0 + s - 1) for r in [1, s - 1) (s - 2 ports).
 *
 * Total port count: |ports| = 2s + 2(s - 2) = 4s - 4.
 *
 * INVARIANT:
 * This canonical ordering is identical across all leaf tiles and internal super-tiles,
 * enabling direct index-based seam alignment between adjacent tiles.
 *
 * @param R0 Top-left row coordinate in global grid.
 * @param C0 Top-left column coordinate in global grid.
 * @param s Tile side length in cells.
 * @return Vector of boundary Point coordinates of size 4s - 4.
 *
 * @complexity
 * Time: O(s)
 * Memory: O(s)
 */
std::vector<Point> getTilePorts(int R0, int C0, int s);

/**
 * @class LeafTile
 * @brief Base-level spatial tile of dimension b x b cells in the H-TRIP hierarchy.
 *
 * Mathematical and Algorithmic Role:
 * A leaf tile represents a local subgrid T = [R0, R0+b) x [C0, C0+b).
 * It encapsulates:
 * 1. The local cell vertices V_T (where |V_T| = b^2).
 * 2. The boundary ports P_T on the tile perimeter (where |P_T| = 4b - 4).
 * 3. The dense intra-tile all-pairs shortest path (APSP) metric D of dimension (b^2) x (b^2).
 * 4. The boundary summary matrix S of dimension |P_T| x |P_T|, formed by projecting D onto P_T:
 *    S_{ij} = D(port_i, port_j).
 *
 * Dynamic Perturbations (Conway Star):
 * When cells inside T become traversable, D is updated purely algebraically using
 * Conway rank-1 min-plus updates (D = min(D, D[:, k] + D[k, :])) with AVX2 SIMD acceleration.
 * This completely avoids running BFS or Dijkstra during online dynamic exploration.
 *
 * Storage Complexity:
 * - D: b^4 * sizeof(dist_t) bytes (e.g., 256^2 * 2 = 128 KB for b = 16).
 * - S: (4b - 4)^2 * sizeof(dist_t) bytes (e.g., 60^2 * 2 = 7.2 KB for b = 16).
 */
class LeafTile {
public:
    int tx = 0;               ///< Tile grid column index (0 <= tx < num_leaf_tiles_x)
    int ty = 0;               ///< Tile grid row index (0 <= ty < num_leaf_tiles_y)
    int R0 = 0;               ///< Global row offset of top-left cell: ty * b
    int C0 = 0;               ///< Global column offset of top-left cell: tx * b
    int b = 16;               ///< Tile side length in cells (typically 16)
    uint32_t version = 0;     ///< Monotonically increasing version counter, incremented when D changes
    int last_pivot_count = 0; ///< Number of active pivot vertices |K| relaxed in the last update

    std::vector<Point> ports;  ///< Canonical perimeter boundary ports (size: 4b - 4)
    std::vector<dist_t> D;     ///< Intra-tile APSP matrix: (b^2) x (b^2) stored row-major
    std::vector<dist_t> S;     ///< Boundary port summary matrix: |ports| x |ports| stored row-major

    LeafTile() = default;

    /**
     * @brief Constructs a leaf tile at grid position (tx, ty) with side length b.
     *
     * Initializes D and S to tropical identity: diagonal elements are 0, off-diagonals are INF.
     *
     * @param tx_ Tile column index.
     * @param ty_ Tile row index.
     * @param b_ Tile side length.
     */
    LeafTile(int tx_, int ty_, int b_);

    /**
     * @brief Maps global grid coordinates (r, c) to local tile cell index in [0, b^2).
     *
     * Local formula: (r - R0) * b + (c - C0).
     */
    [[nodiscard]] int localIdx(int r, int c) const noexcept {
        return (r - R0) * b + (c - C0);
    }

    /**
     * @brief Maps a global Point coordinate to local tile cell index in [0, b^2).
     */
    [[nodiscard]] int localIdx(Point p) const noexcept {
        return localIdx(p.r, p.c);
    }

    /**
     * @brief Maps a local cell index in [0, b^2) to its global Point coordinate.
     */
    [[nodiscard]] Point localPoint(int idx) const noexcept {
        return {R0 + idx / b, C0 + idx % b};
    }

    /**
     * @brief Returns the number of perimeter boundary ports: 4b - 4.
     */
    [[nodiscard]] int numPorts() const noexcept {
        return static_cast<int>(ports.size());
    }

    /**
     * @brief Returns total number of cells in the tile: b^2.
     */
    [[nodiscard]] int numCells() const noexcept {
        return b * b;
    }

    /**
     * @brief Recomputes the complete intra-tile APSP matrix D and boundary summary S from the grid.
     *
     * Algorithm:
     * - If grid is uniform unit-weight: executes b^2 breadth-first searches (one from each passable cell).
     * - If grid has non-unit weights: executes b^2 Dijkstra searches using a min-heap.
     * - Projects resulting cell distances D onto boundary ports to produce summary S.
     *
     * @param grid Global occupancy grid map.
     *
     * @complexity
     * Time: O(b^4) for unit-weight BFS; O(b^4 * log b) for weighted Dijkstra.
     * Memory: O(b^2) auxiliary search queue/heap memory.
     */
    void recomputeFromGrid(const GridMap& grid);

    /**
     * @brief Performs an algebraic Conway rank-1 min-plus update for a single edge (u, v) of weight w.
     *
     * Mathematical Operation:
     * Updates D via the symmetric rank-1 tropical outer product:
     * D' = D (min) [ D(:, u) + w + D(v, :) ] (min) [ D(:, v) + w + D(u, :) ].
     *
     * @param u Local source cell index in [0, b^2).
     * @param v Local target cell index in [0, b^2).
     * @param w Edge weight.
     *
     * @complexity
     * Time: O(|V_active|^2) <= O(b^4), where V_active are cells connected to u or v.
     * Memory: O(b^2) preallocated scratch buffer.
     */
    void updateRank1(int u, int v, dist_t w);

    /**
     * @brief Updates boundary summary S from D, returning any strictly decreased port-to-port distances.
     *
     * Checks each entry S_{ij} against D(port_i, port_j). If strictly smaller, updates S_{ij}
     * and records a PortDecrease event.
     *
     * @return List of PortDecrease records for upward propagation in the hierarchy.
     *
     * @complexity
     * Time: O(p^2) = O(b^2), where p = 4b - 4.
     * Memory: O(p^2) worst-case return vector.
     */
    std::vector<PortDecrease> refreshBoundarySummary();

    /// Reusable scratch vector for storing active cells during rank-1 updates (avoids heap reallocations)
    std::vector<int> active_cells_scratch;

    /**
     * @brief Handles a dynamic cell revelation (cell transitions to Free).
     *
     * Connects cell p to its passable 4-neighbors via rank-1 tropical updates,
     * and returns the resulting boundary port decreases.
     *
     * @param p Coordinate of newly opened cell.
     * @param grid Global occupancy grid.
     * @return List of boundary port decreases.
     *
     * @complexity
     * Time: O(deg * b^4) <= O(b^4), where deg <= 4.
     */
    std::vector<PortDecrease> onCellOpened(Point p, const GridMap& grid);

    /**
     * @brief Batched dynamic cell revelation using vectorized Conway star relaxation.
     *
     * Algorithmic Pipeline:
     * 1. Inserts direct edges between all newly opened cells and their passable 4-neighbors.
     * 2. Gathers the set of distinct altered vertices into pivot list K.
     * 3. Executes an in-place batched Conway star relaxation: for each k in K, row_i is relaxed
     *    using AVX2 SIMD min-plus vector intrinsics (_mm256_adds_epu16, _mm256_min_epu16).
     * 4. Updates boundary summary S and returns decreased port pairs.
     *
     * Zero-Search Invariant:
     * Neither BFS nor Dijkstra is executed during this operation.
     *
     * @param cells List of newly opened cell coordinates belonging to this tile.
     * @param grid Global occupancy grid.
     * @return List of boundary port decreases.
     *
     * @complexity
     * Time: O(|K| * b^4 / 16) with AVX2 SIMD, where |K| <= 4 * |cells|.
     * Memory: O(b^2) auxiliary bitmask and pivot list.
     */
    std::vector<PortDecrease> onBatchCellsOpened(const std::vector<Point>& cells, const GridMap& grid);

    /**
     * @brief Handles a dynamic cell blockage (cell transitions from Free to Obstacle).
     *
     * Because edge deletions are non-monotone in the tropical min-plus semiring,
     * the local metric D is recomputed from the grid via local BFS.
     *
     * @param p Coordinate of blocked cell (unused; recomputeFromGrid inspects grid).
     * @param grid Global occupancy grid.
     * @return List of boundary port decreases (typically empty or altered).
     *
     * @complexity
     * Time: O(b^4).
     */
    std::vector<PortDecrease> onCellBlocked(Point p, const GridMap& grid);
};

} // namespace htrip

