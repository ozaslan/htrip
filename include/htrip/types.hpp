#pragma once

#include <cstdint>
#include <vector>
#include <compare>
#include <cmath>

namespace htrip {

/**
 * @file types.hpp
 * @brief Fundamental mathematical types, tropical semiring arithmetic, spatial coordinates,
 *        and hierarchy configuration for the H-TRIP distance oracle.
 *
 * Mathematical Foundations:
 * H-TRIP operates over the tropical min-plus semiring (R_>=0 U {INF}, min, +), where:
 * - The additive operation is the minimum: a (+) b = min(a, b), with identity INF.
 * - The multiplicative operation is standard addition: a (*) b = a + b, with identity 0.
 *
 * In this semiring, matrix multiplication corresponds to finding shortest paths via intermediate
 * boundary ports: (A (*) B)_{ij} = min_k (A_{ik} + B_{kj}).
 * The reflexive transitive closure (Kleene star) corresponds to all-pairs shortest paths:
 * A^* = I (+) A (+) A^2 (+) ... = min_{m >= 0} A^m.
 */

/**
 * @brief Discrete path length and distance metric representation.
 *
 * Representation Rationale:
 * A 16-bit unsigned integer (uint16_t) is chosen to maximize SIMD packing density.
 * A 256-bit AVX2 register can pack exactly 16 distance words simultaneously, enabling
 * 16-wide vectorized min-plus row relaxations via _mm256_adds_epu16 and _mm256_min_epu16.
 */
using dist_t = uint16_t;

/**
 * @brief Tropical infinity (unreachable vertex or impassable edge).
 *
 * Numerical Rationale:
 * Set to 60000 to remain strictly below UINT16_MAX (65535). This provides a safety margin
 * of 5535 units, preventing integer overflow during intermediate additions before saturation
 * checks clamp the result back to INF. A maximum distance of 60,000 steps comfortably exceeds
 * the diameter of realistic robotic exploration grid workspaces (e.g., 2048 x 2048 grids).
 */
static constexpr dist_t INF = 60000;

/**
 * @brief Default edge traversal cost between 4-connected grid neighbors.
 */
static constexpr dist_t UNIT_STEP = 1;

/**
 * @brief Saturated addition of two distances in the tropical semiring.
 *
 * Implements the semiring product a (*) b = a + b with clamping at INF to prevent
 * arithmetic wrap-around.
 *
 * @param a First distance operand.
 * @param b Second distance operand.
 * @return Saturated sum clamped at INF.
 *
 * @complexity
 * Time: O(1)
 * Memory: O(1)
 */
inline constexpr dist_t add_dist(dist_t a, dist_t b) noexcept {
    if (a >= INF || b >= INF) return INF;
    uint32_t s = static_cast<uint32_t>(a) + static_cast<uint32_t>(b);
    return static_cast<dist_t>(s >= INF ? INF : s);
}

/**
 * @brief Saturated addition of three distances in the tropical semiring.
 *
 * Computes a + b + c clamped at INF, commonly used during three-segment path composition
 * (source-to-port, port-to-port, port-to-target).
 *
 * @param a First distance operand.
 * @param b Second distance operand.
 * @param c Third distance operand.
 * @return Saturated sum clamped at INF.
 *
 * @complexity
 * Time: O(1)
 * Memory: O(1)
 */
inline constexpr dist_t add_dist(dist_t a, dist_t b, dist_t c) noexcept {
    if (a >= INF || b >= INF || c >= INF) return INF;
    uint32_t s = static_cast<uint32_t>(a) + static_cast<uint32_t>(b) + static_cast<uint32_t>(c);
    return static_cast<dist_t>(s >= INF ? INF : s);
}

/**
 * @brief Ternary occupancy state of a discrete grid cell.
 *
 * Represents the robot's belief state during exploration:
 * - Unknown: Unexplored fog-of-war (treated as impassable by the navigation oracle).
 * - Free: Confirmed traversable free space.
 * - Obstacle: Confirmed solid obstacle / wall.
 */
enum class CellState : uint8_t {
    Unknown = 0,
    Free = 1,
    Obstacle = 2
};

/**
 * @brief Discrete 2D integer coordinate in grid space (row, column).
 *
 * Uses row-major ordering where r corresponds to the y-axis (vertical) and
 * c corresponds to the x-axis (horizontal).
 */
struct Point {
    int r = 0;  ///< Row index (0 <= r < rows)
    int c = 0;  ///< Column index (0 <= c < cols)

    constexpr Point() noexcept = default;
    constexpr Point(int r_, int c_) noexcept : r(r_), c(c_) {}

    constexpr bool operator==(const Point& o) const noexcept {
        return r == o.r && c == o.c;
    }

    constexpr bool operator!=(const Point& o) const noexcept {
        return !(*this == o);
    }

    /**
     * @brief Lexicographical ordering (row first, column second).
     */
    constexpr auto operator<=>(const Point& o) const noexcept {
        if (r != o.r) return r <=> o.r;
        return c <=> o.c;
    }

    /**
     * @brief Computes Manhattan distance (L1 norm) between two grid points.
     *
     * @param o Target point.
     * @return |r - o.r| + |c - o.c|
     */
    [[nodiscard]] constexpr int manhattan(const Point& o) const noexcept {
        return std::abs(r - o.r) + std::abs(c - o.c);
    }
};

/**
 * @brief Record of a decreased shortest-path distance between boundary ports of a tile.
 *
 * Emitted when dynamic cell revelations decrease the interior distance between
 * boundary port p_a and boundary port p_b.
 */
struct PortDecrease {
    int p_a = 0;      ///< Index of first boundary port in tile port list
    int p_b = 0;      ///< Index of second boundary port in tile port list
    dist_t w = INF;   ///< New decreased shortest-path distance between p_a and p_b
};

/**
 * @brief Record of a child port distance decrease routed to an internal tree node.
 *
 * Encapsulates the child quadrant/subregion index along with the specific port pair
 * and the newly discovered shorter path weight.
 */
struct ChildDecrease {
    int quad = 0;     ///< Linear child index within parent node (0 <= quad < group_w * group_h)
    int p_a = 0;      ///< Local port index within the child's boundary port set
    int p_b = 0;      ///< Local port index within the child's boundary port set
    dist_t w = INF;   ///< New decreased distance
};

/**
 * @brief Boundary seam edge crossing the boundary between adjacent child subregions.
 *
 * Represents direct 4-connected grid adjacency between a boundary port in one child
 * and an adjacent boundary port in a neighboring child tile.
 */
struct SeamEdge {
    Point p1;             ///< Grid coordinate of port in first child
    Point p2;             ///< Grid coordinate of touching port in second child
    int lca_lvl = 1;      ///< Lowest hierarchy level containing both children
    dist_t w = UNIT_STEP; ///< Edge traversal weight (default: 1)
};

/**
 * @brief Upper bound on the number of child subtrees per internal node.
 *
 * Used both by config validation and by fixed-size child tables in the query scratch.
 * Supports branching factors up to g = 16 (16 x 16 = 256 children).
 */
inline constexpr int kMaxActiveChildren = 256;

/**
 * @brief Geometric hierarchy configuration parameters.
 *
 * Governs the multi-level spatial decomposition of the grid map into recursive
 * super-tiles.
 */
struct HierarchyConfig {
    int tile_size_b = 16;       ///< Side length of base leaf tiles (b in {8, 16, 32, 64})
    int group_w = 2;            ///< Horizontal branching factor (default: 2 for quadtree)
    int group_h = 2;            ///< Vertical branching factor (default: 2 for quadtree)
    int hybrid_k_threshold = 0; ///< Threshold below which queries fall back to k x BFS (default: 0 = pure H-TRIP)

    /**
     * @brief Checks if the configuration represents a standard 2x2 quadtree hierarchy.
     */
    [[nodiscard]] constexpr bool isQuadtree() const noexcept {
        return group_w == 2 && group_h == 2;
    }

    /**
     * @brief Validates geometric parameters and branching invariants.
     *
     * INVARIANT:
     * - tile_size_b > 0
     * - group_w >= 2 and group_h >= 2
     * - group_w == group_h (square aggregation branching is mandatory)
     * - group_w * group_h <= kMaxActiveChildren
     */
    [[nodiscard]] constexpr bool isValid() const noexcept {
        return tile_size_b > 0 && group_w >= 2 && group_h >= 2 && group_w == group_h
            && group_w * group_h <= kMaxActiveChildren;
    }
};

} // namespace htrip

