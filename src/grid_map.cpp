#include "htrip/grid_map.hpp"
#include <fstream>
#include <sstream>
#include <queue>
#include <algorithm>
#include <iostream>

namespace htrip {

GridMap::GridMap(int r, int c, CellState initial_state)
    : rows(r), cols(c),
      ground_truth(static_cast<size_t>(r * c), CellState::Free),
      observed(static_cast<size_t>(r * c), initial_state),
      weights(static_cast<size_t>(r * c), 1) {}

bool GridMap::loadMovingAI(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) return false;

    std::string line;
    int h = 0, w = 0;

    // Parse MovingAI metadata header
    // Expected format:
    //   type octile
    //   height <H>
    //   width <W>
    //   map
    while (std::getline(file, line)) {
        // Strip trailing line delimiters (\r, \n, whitespace)
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (line.rfind("height", 0) == 0) {
            std::stringstream ss(line);
            std::string key;
            ss >> key >> h;
        } else if (line.rfind("width", 0) == 0) {
            std::stringstream ss(line);
            std::string key;
            ss >> key >> w;
        } else if (line == "map") {
            // Header finished; grid ASCII rows begin immediately below
            break;
        }
    }

    if (h <= 0 || w <= 0) return false;

    rows = h;
    cols = w;
    const size_t nCells = static_cast<size_t>(rows * cols);
    ground_truth.assign(nCells, CellState::Obstacle);
    observed.assign(nCells, CellState::Unknown);
    weights.assign(nCells, 1);

    // Read grid rows character-by-character
    // Passable terrain characters:
    //   '.' : normal terrain
    //   'G' : grass / open space
    //   'S' : swamp / traversable with cost
    // Impassable terrain:
    //   '@', 'O', 'T', 'W' : trees, out-of-bounds, walls, deep water
    for (int r = 0; r < rows; ++r) {
        if (!std::getline(file, line)) break;
        for (int c = 0; c < cols && c < static_cast<int>(line.size()); ++c) {
            char ch = line[static_cast<size_t>(c)];
            if (ch == '.' || ch == 'G' || ch == 'S') {
                ground_truth[index(r, c)] = CellState::Free;
            } else {
                ground_truth[index(r, c)] = CellState::Obstacle;
            }
        }
    }
    return true;
}

void GridMap::generateProceduralRooms(int r, int c, int room_size, uint64_t seed) {
    if (room_size < 4) room_size = 4;
    rows = r;
    cols = c;
    const size_t nCells = static_cast<size_t>(rows * cols);
    ground_truth.assign(nCells, CellState::Obstacle);
    observed.assign(nCells, CellState::Unknown);
    weights.assign(nCells, 1);

    std::mt19937_64 rng(seed);

    const int num_rooms_y = rows / room_size;
    const int num_rooms_x = cols / room_size;

    // 1. Carve rectangular room interiors.
    // Each room spans coordinates [ry * room_size + 1, (ry + 1) * room_size - 1),
    // leaving a 1-cell border between room boundaries.
    for (int ry = 0; ry < num_rooms_y; ++ry) {
        for (int rx = 0; rx < num_rooms_x; ++rx) {
            int r_start = ry * room_size + 1;
            int r_end = std::min((ry + 1) * room_size - 1, rows - 1);
            int c_start = rx * room_size + 1;
            int c_end = std::min((rx + 1) * room_size - 1, cols - 1);

            for (int y = r_start; y < r_end; ++y) {
                for (int x = c_start; x < c_end; ++x) {
                    ground_truth[index(y, x)] = CellState::Free;
                }
            }
        }
    }

    // 2. Carve doorways connecting adjacent rooms.
    // Because rooms are separated by 2-cell thick walls (the boundary cells of both rooms),
    // a doorway must pierce both wall layers (col1 and col2, or row1 and row2) to form
    // a connected passage.
    for (int ry = 0; ry < num_rooms_y; ++ry) {
        for (int rx = 0; rx < num_rooms_x; ++rx) {
            // Horizontal doorway to east adjacent room
            if (rx + 1 < num_rooms_x) {
                int col1 = (rx + 1) * room_size - 1;
                int col2 = (rx + 1) * room_size;
                int door_row = ry * room_size + 1 + static_cast<int>(rng() % static_cast<uint64_t>(room_size - 2));
                if (door_row < rows) {
                    if (col1 < cols) ground_truth[index(door_row, col1)] = CellState::Free;
                    if (col2 < cols) ground_truth[index(door_row, col2)] = CellState::Free;
                }
            }
            // Vertical doorway to south adjacent room
            if (ry + 1 < num_rooms_y) {
                int row1 = (ry + 1) * room_size - 1;
                int row2 = (ry + 1) * room_size;
                int door_col = rx * room_size + 1 + static_cast<int>(rng() % static_cast<uint64_t>(room_size - 2));
                if (door_col < cols) {
                    if (row1 < rows) ground_truth[index(row1, door_col)] = CellState::Free;
                    if (row2 < rows) ground_truth[index(row2, door_col)] = CellState::Free;
                }
            }
        }
    }

    // For synthetic room benchmarks, observed is initially initialized identical to ground truth
    observed = ground_truth;
}

bool GridMap::initExplorationSnapshot(Point start_point, double target_coverage,
                                     double island_prob, uint64_t seed) {
    // Ensure starting point is valid and traversable
    if (!inBounds(start_point) || ground_truth[index(start_point)] != CellState::Free) {
        bool found = false;
        for (int r = 0; r < rows && !found; ++r) {
            for (int c = 0; c < cols && !found; ++c) {
                if (ground_truth[index(r, c)] == CellState::Free) {
                    start_point = {r, c};
                    found = true;
                }
            }
        }
        if (!found) return false;
    }

    // Initialize the observed belief layer to complete fog-of-war (Unknown)
    observed.assign(static_cast<size_t>(rows * cols), CellState::Unknown);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Floodfill BFS wavefront simulation:
    // Simulates a robot expanding its explored territory contiguously from start_point.
    std::queue<Point> q;
    std::vector<bool> visited(static_cast<size_t>(rows * cols), false);

    q.push(start_point);
    visited[static_cast<size_t>(index(start_point))] = true;
    observed[static_cast<size_t>(index(start_point))] = CellState::Free;

    int total_free = 0;
    for (CellState s : ground_truth) {
        if (s == CellState::Free) total_free++;
    }
    const int max_to_explore = std::max(10, static_cast<int>(total_free * target_coverage));
    int explored_free = 1;

    // 4-connected cardinal neighbor offsets: North, South, West, East
    static const int dr[4] = {-1, 1, 0, 0};
    static const int dc[4] = {0, 0, -1, 1};

    while (!q.empty() && explored_free < max_to_explore) {
        Point curr = q.front();
        q.pop();

        for (int i = 0; i < 4; ++i) {
            int nr = curr.r + dr[i];
            int nc = curr.c + dc[i];
            if (!inBounds(nr, nc)) continue;

            int nidx = index(nr, nc);
            if (visited[static_cast<size_t>(nidx)]) continue;
            visited[static_cast<size_t>(nidx)] = true;

            CellState gt = ground_truth[static_cast<size_t>(nidx)];
            if (gt == CellState::Obstacle) {
                // Obstacles forming the boundary of visited free cells become revealed
                observed[static_cast<size_t>(nidx)] = CellState::Obstacle;
            } else if (gt == CellState::Free) {
                // Stochastic corridor suppression:
                // Allows certain corridors to remain unexplored, producing realistic unreached "islands"
                if (island_prob > 0.0 && dist(rng) < island_prob && explored_free > 50) {
                    continue; // Skip propagation; leaves downstream space Unknown
                }
                observed[static_cast<size_t>(nidx)] = CellState::Free;
                explored_free++;
                q.push({nr, nc});
                if (explored_free >= max_to_explore) break;
            }
        }
    }

    return true;
}

bool GridMap::openDoorway(Point p) {
    if (!inBounds(p)) return false;
    observed[static_cast<size_t>(index(p))] = CellState::Free;
    return true;
}

bool GridMap::blockCorridor(Point p) {
    if (!inBounds(p)) return false;
    observed[static_cast<size_t>(index(p))] = CellState::Obstacle;
    return true;
}

bool GridMap::revealCell(Point p) {
    if (!inBounds(p)) return false;
    observed[static_cast<size_t>(index(p))] = ground_truth[static_cast<size_t>(index(p))];
    return true;
}

} // namespace htrip

