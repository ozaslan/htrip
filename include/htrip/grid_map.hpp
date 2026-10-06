#pragma once

#include "types.hpp"
#include <vector>
#include <string>
#include <random>

namespace htrip {

/**
 * @class GridMap
 * @brief Discrete 2D occupancy grid representing both ground-truth reality and the robot's
 *        partially observed belief state.
 *
 * Mathematical and Structural Role:
 * Represents the 4-connected spatial lattice G = (V, E) embedded in R^2.
 * Vertices correspond to integer grid coordinates (r, c).
 * Directed edges exist between 4-connected spatial neighbors with weights determined by terrain cost.
 *
 * The map maintains dual state layers:
 * 1. Ground Truth (ground_truth): The underlying physical environment (Free vs. Obstacle).
 * 2. Observed Belief (observed): The robot's online perceived map at time t, containing:
 *    - Unknown (fog-of-war, unexplored),
 *    - Free (confirmed traversable),
 *    - Obstacle (confirmed solid boundary).
 *
 * Memory Layout:
 * All grids are stored as flat 1D contiguous vectors of length rows * cols in row-major order:
 * index(r, c) = r * cols + c.
 */
class GridMap {
public:
    int rows = 0;  ///< Vertical dimension (height) of the grid
    int cols = 0;  ///< Horizontal dimension (width) of the grid

    /**
     * @brief Complete ground truth environment map (Free or Obstacle).
     *
     * Used by the simulation harness to reveal cells as the robot's sensors advance.
     */
    std::vector<CellState> ground_truth;

    /**
     * @brief The robot's current online belief state (Unknown, Free, or Obstacle).
     *
     * The H-TRIP distance oracle builds its graph topology strictly from this observed layer.
     */
    std::vector<CellState> observed;

    /**
     * @brief Costmap traversal weights for weighted shortest-path planning.
     *
     * Default value is 1 for uniform terrain. Non-unit values model rough terrain, slopes,
     * or risk zones. Impassable cells have effective weight INF.
     */
    std::vector<dist_t> weights;

    GridMap() = default;

    /**
     * @brief Constructs a grid map with specified dimensions and uniform initial state.
     *
     * @param r Number of rows.
     * @param c Number of columns.
     * @param initial_state Initial state for all cells in the observed grid (default: Unknown).
     */
    GridMap(int r, int c, CellState initial_state = CellState::Unknown);

    /**
     * @brief Checks whether discrete coordinates fall strictly within the grid boundaries.
     *
     * @param r Row index.
     * @param c Column index.
     * @return True iff 0 <= r < rows and 0 <= c < cols.
     *
     * @complexity
     * Time: O(1)
     * Memory: O(1)
     */
    [[nodiscard]] bool inBounds(int r, int c) const noexcept {
        return r >= 0 && r < rows && c >= 0 && c < cols;
    }

    /**
     * @brief Checks whether a Point coordinate falls strictly within the grid boundaries.
     *
     * @param p Point coordinate.
     * @return True iff p is within bounds.
     */
    [[nodiscard]] bool inBounds(Point p) const noexcept {
        return inBounds(p.r, p.c);
    }

    /**
     * @brief Maps a 2D coordinate (r, c) to a flat 1D vector index.
     *
     * Row-major index formula: idx = r * cols + c.
     *
     * @param r Row index.
     * @param c Column index.
     * @return 1D array index.
     */
    [[nodiscard]] int index(int r, int c) const noexcept {
        return r * cols + c;
    }

    /**
     * @brief Maps a Point coordinate to a flat 1D vector index.
     *
     * @param p Point coordinate.
     * @return 1D array index.
     */
    [[nodiscard]] int index(Point p) const noexcept {
        return index(p.r, p.c);
    }

    /**
     * @brief Reconstructs 2D grid coordinates from a flat 1D array index.
     *
     * @param idx Flat array index.
     * @return Point with r = idx / cols, c = idx % cols.
     */
    [[nodiscard]] Point point(int idx) const noexcept {
        return {idx / cols, idx % cols};
    }

    /**
     * @brief Evaluates whether a cell is traversable by the robot.
     *
     * INVARIANT:
     * In H-TRIP, only confirmed Free cells in the observed grid are passable.
     * Unknown cells (unexplored space) and Obstacle cells are strictly impassable.
     *
     * @param r Row index.
     * @param c Column index.
     * @return True iff cell is in bounds and observed as CellState::Free.
     *
     * @complexity
     * Time: O(1)
     * Memory: O(1)
     */
    [[nodiscard]] bool isPassable(int r, int c) const noexcept {
        if (!inBounds(r, c)) return false;
        return observed[index(r, c)] == CellState::Free;
    }

    /**
     * @brief Evaluates whether a Point coordinate is traversable.
     *
     * @param p Point coordinate.
     * @return True iff p is passable.
     */
    [[nodiscard]] bool isPassable(Point p) const noexcept {
        return isPassable(p.r, p.c);
    }

    /**
     * @brief Computes discrete obstacle clearance in the observed grid.
     *
     * Clearance values:
     * - 0: Impassable (Obstacle, Unknown, or out-of-bounds).
     * - 1: Immediately adjacent to an obstacle (inscribed collision zone).
     * - 2: Buffer safety margin (within 2-hop Chebyshev distance of an obstacle).
     * - 3: Wide open corridor (at least 2 cells away from nearest obstacle).
     *
     * @param r Row index.
     * @param c Column index.
     * @return Discrete clearance tier in {0, 1, 2, 3}.
     *
     * @complexity
     * Time: O(1) bounded neighborhood inspection (up to 24 cells).
     * Memory: O(1).
     */
    [[nodiscard]] uint8_t getClearance(int r, int c) const noexcept {
        if (!isPassable(r, c)) return 0;

        // Check 1-hop 8-neighbors (clearance <= 1)
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                if (dr == 0 && dc == 0) continue;
                int nr = r + dr;
                int nc = c + dc;
                if (!inBounds(nr, nc) || observed[static_cast<size_t>(index(nr, nc))] == CellState::Obstacle) {
                    return 1;
                }
            }
        }

        // Check 2-hop Chebyshev neighborhood (clearance <= 2)
        for (int dr = -2; dr <= 2; ++dr) {
            for (int dc = -2; dc <= 2; ++dc) {
                if (std::abs(dr) < 2 && std::abs(dc) < 2) continue;
                int nr = r + dr;
                int nc = c + dc;
                if (!inBounds(nr, nc) || observed[static_cast<size_t>(index(nr, nc))] == CellState::Obstacle) {
                    return 2;
                }
            }
        }

        return 3;
    }

    /**
     * @brief Computes discrete obstacle clearance for a Point coordinate.
     *
     * @param p Point coordinate.
     * @return Discrete clearance tier in {0, 1, 2, 3}.
     */
    [[nodiscard]] uint8_t getClearance(Point p) const noexcept {
        return getClearance(p.r, p.c);
    }

    /**
     * @brief Accesses the observed cell state at (r, c).
     *
     * @param r Row index.
     * @param c Column index.
     * @return Observed CellState, or CellState::Obstacle if out of bounds.
     */
    [[nodiscard]] CellState getObserved(int r, int c) const noexcept {
        if (!inBounds(r, c)) return CellState::Obstacle;
        return observed[index(r, c)];
    }

    /**
     * @brief Accesses the observed cell state at Point p.
     *
     * @param p Point coordinate.
     * @return Observed CellState.
     */
    [[nodiscard]] CellState getObserved(Point p) const noexcept {
        return getObserved(p.r, p.c);
    }

    /**
     * @brief Sets the observed cell state at (r, c).
     *
     * @param r Row index.
     * @param c Column index.
     * @param s New observed state.
     */
    void setObserved(int r, int c, CellState s) noexcept {
        if (inBounds(r, c)) observed[index(r, c)] = s;
    }

    /**
     * @brief Sets the observed cell state at Point p.
     *
     * @param p Point coordinate.
     * @param s New observed state.
     */
    void setObserved(Point p, CellState s) noexcept {
        setObserved(p.r, p.c, s);
    }

    /**
     * @brief Checks if non-unit terrain weights have been configured.
     *
     * @return True if costmap weights array is populated.
     */
    [[nodiscard]] bool hasWeights() const noexcept {
        return !weights.empty();
    }

    /**
     * @brief Retrieves the terrain traversal weight for a cell.
     *
     * @param r Row index.
     * @param c Column index.
     * @return Cell traversal cost, or INF if out of bounds or an obstacle.
     */
    [[nodiscard]] dist_t getCellWeight(int r, int c) const noexcept {
        if (!inBounds(r, c)) return INF;
        if (observed[static_cast<size_t>(index(r, c))] == CellState::Obstacle) return INF;
        if (weights.empty()) return 1;
        return weights[static_cast<size_t>(index(r, c))];
    }

    /**
     * @brief Retrieves the terrain traversal weight for a Point coordinate.
     *
     * @param p Point coordinate.
     * @return Cell traversal cost.
     */
    [[nodiscard]] dist_t getCellWeight(Point p) const noexcept {
        return getCellWeight(p.r, p.c);
    }

    /**
     * @brief Sets the terrain traversal cost for a cell.
     *
     * @param r Row index.
     * @param c Column index.
     * @param w Traversal cost.
     */
    void setCellWeight(int r, int c, dist_t w) noexcept {
        if (!inBounds(r, c)) return;
        if (weights.size() < static_cast<size_t>(rows * cols)) {
            weights.assign(static_cast<size_t>(rows * cols), 1);
        }
        weights[static_cast<size_t>(index(r, c))] = w;
    }

    /**
     * @brief Sets the terrain traversal cost for a Point coordinate.
     *
     * @param p Point coordinate.
     * @param w Traversal cost.
     */
    void setCellWeight(Point p, dist_t w) noexcept {
        setCellWeight(p.r, p.c, w);
    }

    /**
     * @brief Computes the undirected edge weight bridging two adjacent 4-connected cells.
     *
     * Mathematical Formula:
     * The edge traversal cost between adjacent cells u and v is computed as the arithmetic
     * average of their individual cell weights: w(u, v) = max(1, (w_u + w_v) / 2).
     * If either cell is impassable (weight >= INF), returns INF.
     *
     * @param u First endpoint.
     * @param v Second endpoint.
     * @return Symmetric traversal cost between u and v.
     */
    [[nodiscard]] dist_t edgeWeight(Point u, Point v) const noexcept {
        dist_t wu = getCellWeight(u);
        dist_t wv = getCellWeight(v);
        if (wu >= INF || wv >= INF) return INF;
        uint32_t sum = static_cast<uint32_t>(wu) + static_cast<uint32_t>(wv);
        dist_t avg = static_cast<dist_t>(sum / 2);
        return avg == 0 ? 1 : avg;
    }

    /**
     * @brief Loads a standard MovingAI .map benchmark file into ground truth.
     *
     * Parses MovingAI headers ("type octile", "height H", "width W", "map").
     * Recognizes '.', 'G', 'S' as Free cells, and all other characters ('@', 'T', 'O')
     * as Obstacle cells. Sets observed grid to all-Unknown.
     *
     * @param filepath Path to the .map file.
     * @return True iff file was successfully opened and parsed.
     *
     * @complexity
     * Time: O(rows * cols)
     * Memory: O(rows * cols)
     */
    bool loadMovingAI(const std::string& filepath);

    /**
     * @brief Procedurally generates a grid world composed of regular rooms separated by walls.
     *
     * Carves room interiors of dimension room_size x room_size, leaving 1-cell thick separating
     * walls, and inserts random doorways connecting adjacent rooms.
     * Sets observed grid equal to ground truth.
     *
     * @param r Grid rows.
     * @param c Grid columns.
     * @param room_size Room dimension in cells (minimum 4).
     * @param seed Random generator seed.
     *
     * @complexity
     * Time: O(r * c)
     * Memory: O(r * c)
     */
    void generateProceduralRooms(int r, int c, int room_size = 16, uint64_t seed = 42);

    /**
     * @brief Generates a realistic partially explored exploration snapshot at time t.
     *
     * Executes a bounded breadth-first search from start_point until reaching target_coverage
     * of all free cells in ground truth. Explored free space and adjacent obstacle boundaries
     * are marked observed, while the remainder of the environment remains Unknown (fog-of-war).
     *
     * @param start_point Exploration origin coordinate.
     * @param target_coverage Fraction of free cells to uncover (e.g., 0.45 for 45% explored).
     * @param island_prob Probability of skipping corridors to create unobserved interior islands.
     * @param seed Random generator seed.
     * @return True if successfully initialized.
     *
     * @complexity
     * Time: O(rows * cols)
     * Memory: O(rows * cols)
     */
    bool initExplorationSnapshot(Point start_point, double target_coverage = 0.45,
                                double island_prob = 0.1, uint64_t seed = 12345);

    /**
     * @brief Dynamically opens a doorway cell (revelation event).
     *
     * Transitions observed state from Obstacle/Unknown to Free.
     * Used by benchmarks to simulate dynamic shortcut discoveries.
     *
     * @param p Coordinate of cell to open.
     * @return True iff cell was in bounds.
     */
    bool openDoorway(Point p);

    /**
     * @brief Dynamically blocks a corridor cell (blockage event).
     *
     * Transitions observed state from Free to Obstacle.
     * Used by benchmarks to simulate dynamic obstacle blockages.
     *
     * @param p Coordinate of cell to block.
     * @return True iff cell was in bounds.
     */
    bool blockCorridor(Point p);

    /**
     * @brief Reveals the true underlying state of a cell from ground truth.
     *
     * Transitions observed state from Unknown to ground_truth[p].
     *
     * @param p Coordinate of cell to reveal.
     * @return True iff cell was in bounds.
     */
    bool revealCell(Point p);
};

} // namespace htrip

