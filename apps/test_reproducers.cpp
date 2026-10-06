#include "htrip/minplus_oracle.hpp"
#include "htrip/bfs_oracle.hpp"
#include "htrip/dijkstra_oracle.hpp"
#include "htrip/fuel_oracle.hpp"
#include "htrip/grid_map.hpp"
#include <iostream>
#include <vector>
#include <cassert>

using namespace htrip;

// Color helpers
static const char* GREEN = "\033[32m";
static const char* RED   = "\033[31m";
static const char* RESET = "\033[0m";

// -----------------------------------------------------------------------------
// Test 1.A: Weighted cell opening must use grid.edgeWeight, NOT UNIT_STEP
// -----------------------------------------------------------------------------
static bool testBug1A_WeightedCellOpen() {
    std::cout << "[Test 1.A: Weighted Cell Opening (onCellOpened)]\n";
    // 8x8 corridor with uniform weight 5, b=2, g=2
    GridMap grid(8, 8);
    for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < 8; ++c) {
            grid.setObserved(r, c, CellState::Free);
            grid.setCellWeight(r, c, 5);
        }
    }

    HierarchyConfig hc;
    hc.tile_size_b = 2;
    hc.group_w = 2;
    hc.group_h = 2;

    MinPlusOracle oracle(hc);
    if (!oracle.build(grid)) return false;

    // Block cell (0, 3) then reopen it
    Point p(0, 3);
    grid.setObserved(p, CellState::Obstacle);
    oracle.onCellBlocked(p, grid);

    // Re-open cell (0, 3)
    grid.setObserved(p, CellState::Free);
    oracle.onCellOpened(p, grid);

    // Query distance between (0, 0) and (0, 5)
    // Distance should be 5 steps across weight-5 cells = 5 * 5 = 25.
    Point s(0, 0), t(0, 5);
    dist_t d_htrip = oracle.queryDistance(s, t, grid);
    dist_t d_dijkstra = DijkstraOracle::queryDistance(grid, s, t);

    bool ok = (d_htrip == d_dijkstra && d_htrip == 25);
    std::cout << "  (0,0) -> (0,5): H-TRIP=" << d_htrip << ", Dijkstra=" << d_dijkstra
              << " -> " << (ok ? GREEN : RED) << (ok ? "PASSED" : "FAILED (Bug 1.A Confirmed: unit step used)") << RESET << "\n";
    return ok;
}

// -----------------------------------------------------------------------------
// Test 1.B: Early break in onCellBlocked / onCellWeightChanged leaves stale seam edges
// -----------------------------------------------------------------------------
static bool testBug1B_ExteriorSeamEarlyBreak() {
    std::cout << "\n[Test 1.B: Exterior Seam Early Break in Ancestor Propagation]\n";
    // 32x32 grid, b=8, g=2. Bridge cell at (15, 15) connects quadrant 1 to quadrant 2.
    // Inside Node(0,0), (15, 15) is disconnected from other exterior ports,
    // so parent.S_U == old_S_U causes premature break, leaving stale seam edge in grandparent.
    GridMap grid(32, 32);
    for (int r = 0; r < 32; ++r) {
        for (int c = 0; c < 32; ++c) {
            grid.setObserved(r, c, CellState::Free);
        }
    }
    grid.setObserved(14, 15, CellState::Obstacle);
    grid.setObserved(15, 14, CellState::Obstacle);
    grid.setObserved(16, 16, CellState::Obstacle);

    for (int c = 17; c < 32; ++c) grid.setObserved(15, c, CellState::Obstacle);
    for (int r = 17; r < 32; ++r) grid.setObserved(r, 15, CellState::Obstacle);

    HierarchyConfig hc;
    hc.tile_size_b = 8;
    hc.group_w = 2;
    hc.group_h = 2;

    MinPlusOracle oracle(hc);
    if (!oracle.build(grid)) return false;

    // Initial path between (16, 15) and (15, 16) goes through (15, 15): distance 2
    Point p1(16, 15), p2(15, 16);
    dist_t d_init = oracle.queryDistance(p1, p2, grid);

    // Block the bridge at (15, 15)
    Point bridge(15, 15);
    grid.setObserved(bridge, CellState::Obstacle);
    oracle.onCellBlocked(bridge, grid);

    dist_t d_htrip = oracle.queryDistance(p1, p2, grid);
    dist_t d_bfs = BfsOracle::queryDistance(grid, p1, p2);

    bool ok = (d_htrip == d_bfs && d_htrip > d_init);
    std::cout << "  After blocking bridge (15,15): (16,15)->(15,16) H-TRIP=" << d_htrip
              << ", BFS=" << d_bfs
              << " -> " << (ok ? GREEN : RED) << (ok ? "PASSED" : "FAILED (Bug 1.B Confirmed: stale seam edge in ancestor)") << RESET << "\n";
    return ok;
}

// -----------------------------------------------------------------------------
// Test 1.C: Hybrid fallback must be weight-aware
// -----------------------------------------------------------------------------
static bool testBug1C_HybridFallbackWeighted() {
    std::cout << "\n[Test 1.C: Hybrid Fallback on Weighted Grid]\n";
    GridMap grid(16, 16);
    for (int r = 0; r < 16; ++r) {
        for (int c = 0; c < 16; ++c) {
            grid.setObserved(r, c, CellState::Free);
            grid.setCellWeight(r, c, 5);
        }
    }

    HierarchyConfig hc;
    hc.tile_size_b = 4;
    hc.group_w = 2;
    hc.group_h = 2;
    hc.hybrid_k_threshold = 4; // Trigger fallback for k <= 4

    MinPlusOracle oracle(hc);
    if (!oracle.build(grid)) return false;

    std::vector<Point> vps = {{0, 0}, {0, 3}, {0, 6}}; // k = 3 <= 4
    std::vector<dist_t> htrip_mat = oracle.queryAllPairs(vps, grid);
    std::vector<dist_t> dijkstra_mat = DijkstraOracle::queryAllPairs(grid, vps);

    bool ok = (htrip_mat == dijkstra_mat);
    std::cout << "  hybrid_k_threshold=4, k=3: H-TRIP(0,2)=" << htrip_mat[2]
              << ", Dijkstra(0,2)=" << dijkstra_mat[2]
              << " -> " << (ok ? GREEN : RED) << (ok ? "PASSED" : "FAILED (Bug 1.C Confirmed: BFS hop count returned on weighted grid)") << RESET << "\n";
    return ok;
}

// -----------------------------------------------------------------------------
// Test 1.D: FuelOracle preserves surviving pairs when k changes
// -----------------------------------------------------------------------------
static bool testBug1D_FuelPreserveSurvivingPairs() {
    std::cout << "\n[Test 1.D: FuelOracle Cache Preservation on k/Viewpoint List Change]\n";
    GridMap grid(32, 32);
    for (int r = 0; r < 32; ++r) {
        for (int c = 0; c < 32; ++c) {
            grid.setObserved(r, c, CellState::Free);
        }
    }

    FuelCacheOracle fuel(2); // sensor radius = 2
    std::vector<Point> vps1 = {{5, 5}, {5, 15}, {5, 25}};
    fuel.build(grid, vps1);

    dist_t initial_dist = fuel.queryDistance(0, 1); // dist between (5,5) and (5,15)

    // Now k changes: a 4th point is added at (30,30), but (5,5) and (5,15) are outside the sensor box at (30,30)
    std::vector<Point> vps2 = {{5, 5}, {5, 15}, {5, 25}, {30, 30}};
    // Sensor update at (30,30) with box_radius=2
    fuel.onSensorUpdate(Point(30, 30), 2, grid, vps2);

    dist_t preserved_dist = fuel.queryDistance(0, 1);
    bool ok = (preserved_dist == initial_dist && preserved_dist < INF);

    std::cout << "  k changed from 3 to 4: cached d(P0, P1) was " << initial_dist
              << ", now " << preserved_dist
              << " -> " << (ok ? GREEN : RED) << (ok ? "PASSED" : "FAILED (Bug 1.D Confirmed: cache wiped on k change)") << RESET << "\n";
    return ok;
}

int main() {
    std::cout << "=================================================================\n";
    std::cout << "  ICRA 2027: REPRODUCER TESTS FOR REPORTED CODE REVIEW BUGS\n";
    std::cout << "=================================================================\n\n";

    bool b1a = testBug1A_WeightedCellOpen();
    bool b1b = testBug1B_ExteriorSeamEarlyBreak();
    bool b1c = testBug1C_HybridFallbackWeighted();
    bool b1d = testBug1D_FuelPreserveSurvivingPairs();

    int failed = (!b1a) + (!b1b) + (!b1c) + (!b1d);
    std::cout << "\n-----------------------------------------------------------------\n";
    std::cout << "  Reproducer Summary: " << failed << " bugs confirmed.\n";
    std::cout << "=================================================================\n";
    return failed == 0 ? 0 : 1;
}
