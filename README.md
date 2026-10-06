# HTRIP

Hierarchical Tropical Rank-one Incremental Propagation

![HTRIP clay letters standing on a folded map, with routes drawn across the letters](images/htrip-banner.png)

HTRIP (Hierarchical Tropical Rank-one Incremental Propagation) is a C++20 library for exact many-to-many shortest paths on a growing occupancy grid. It partitions the known map into tiles, stores port-to-port distances in the min-plus semiring, and inserts each newly revealed free cell as a rank-one update that propagates only decreased summaries to ancestor tiles. A shared hierarchical lift then fills the dense distance matrix among the active robot poses and frontier viewpoints.

This repository is that library. It contains the sources, the tests, the command-line benchmarks, the six Moving AI maps used in the evaluation, and five exploration recordings. The index kernels require AVX2. The library compiles on its own, without the H-BRICK sources.

## Using the library

There is no installed package. Add this repository as a subdirectory and link `htrip`.

```cmake
add_subdirectory(htrip)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE htrip)
```

The public headers live under `include/htrip/` and the namespace is `htrip`. A Moving AI file stores the ground-truth floorplan and leaves every cell unobserved. The fully revealed experiments copy that floorplan into the observed grid before building the index. `htrip::campaign::loadFullyRevealed` does both steps.

Load Boston at the published leaf side 16 and group size 2, ask for the distance between two cells, and compare it with one BFS:

```cpp
#include "htrip/bfs_oracle.hpp"
#include "htrip/campaign_support.hpp"
#include "htrip/minplus_oracle.hpp"

htrip::GridMap grid;
htrip::campaign::loadFullyRevealed(
    "datasets/movingai/street/maps/Boston_0_256.map", grid);

htrip::HierarchyConfig config;
config.tile_size_b = 16;
config.group_w = 2;
config.group_h = 2;

htrip::MinPlusOracle index(config);
index.build(grid);

const std::vector<htrip::Point> viewpoints{
    htrip::Point{20, 20},
    htrip::Point{80, 40},
};
const std::vector<htrip::dist_t> matrix = index.queryAllPairs(viewpoints, grid);
const htrip::dist_t indexed = matrix[0 * viewpoints.size() + 1];
const htrip::dist_t bfs = htrip::BfsOracle::queryDistance(
    grid, viewpoints[0], viewpoints[1]);
```

`matrix` is row-major, so entry `i * k + j` is the distance from viewpoint `i` to viewpoint `j`. `indexed` and `bfs` are equal when both cells lie on free space. The unreachable value is `htrip::INF`, which is 60000.

The same call fills the whole matrix. For k viewpoints the length of `matrix` is `k * k`:

```cpp
const std::vector<htrip::dist_t> all = index.queryAllPairs(viewpoints, grid);
const std::vector<htrip::dist_t> reference =
    htrip::BfsOracle::queryAllPairs(grid, viewpoints);
```

A hot query can pass a caller-owned buffer so the call does not allocate. `HTripQueryScratch` is sized for the viewpoint count and the hierarchy depth:

```cpp
htrip::MinPlusOracle::HTripQueryScratch scratch;
std::vector<htrip::dist_t> out;
index.queryAllPairs(viewpoints, grid, out, scratch);
```

When exploration opens a cell, write it into the observed grid and then tell the index. The update touches the leaf that contains the cell and propagates decreased port distances to its ancestors:

```cpp
htrip::GridMap room(32, 32, htrip::CellState::Obstacle);
for (int r = 1; r < 31; ++r) {
    for (int c = 1; c < 16; ++c) {
        room.setObserved(r, c, htrip::CellState::Free);
    }
}

htrip::HierarchyConfig room_config;
room_config.tile_size_b = 8;
room_config.group_w = 2;
room_config.group_h = 2;
htrip::MinPlusOracle room_index(room_config);
room_index.build(room);

const htrip::Point opened{8, 20};
room.setObserved(opened, htrip::CellState::Free);
room_index.onCellOpened(opened, room);
```

`onBatchCellsOpened` applies the same update to several cells at once. That is the call timed in the opening experiment.

## Maps

The evaluation uses six Moving AI occupancy grids. The files are already unpacked at `datasets/movingai/<set>/maps/<file>.map`.

| Paper id | Set | Map file | Leaf b | Group g |
|----------|-----|----------|--------|---------|
| M1 | street | Boston_0_256.map | 16 | 2 |
| M2 | street | Boston_0_1024.map | 32 | 2 |
| M3 | maze | maze512-1-0.map | 16 | 2 |
| M4 | maze | maze512-16-0.map | 16 | 2 |
| M5 | street | Paris_0_512.map | 8 | 2 |
| M6 | street | Berlin_0_1024.map | 8 | 2 |

Passable-cell counts and the street-map credit are in [`datasets/movingai/README.md`](datasets/movingai/README.md). The two maze files and the Moving AI citation for that set are in [`datasets/movingai/maze/README.md`](datasets/movingai/maze/README.md).

## Requirements

- A C++20 compiler with AVX2 (GCC 11 or newer, or Clang 14 or newer)
- CMake 3.20 or newer

The Release preset compiles with `-march=native`. The min-plus kernels refuse to compile when `__AVX2__` is absent. Google Benchmark is downloaded the first time CMake configures `bench_micro`.

## Build

From the repository root:

```bash
cmake --preset dev
cmake --build --preset dev
```

Binaries land in `build/dev/`. The measurement programs are:

- `bench_query_vs_k` — many-to-many query time against k times BFS on a fully revealed map
- `bench_update_latency` — index update after opening free cells
- `run_multi_robot_sim` — sixteen-robot exploration with paired query timings and no video unless `--record` is passed
- `campaign_e0` and `campaign_exactness` — build and exactness checks
- `bench_frontier_oracle`, `bench_weighted`, `bench_systematic_bfs`, `bench_micro` — additional checks

## Tests

The default build compiles every test. CTest runs them from the build directory:

```bash
ctest --preset dev --output-on-failure
```

A single case:

```bash
ctest --test-dir build/dev --output-on-failure -R test_exactness
```

The registered cases are `test_exactness`, `test_fuel_baseline`, `test_property_based`, `test_weighted`, `test_reproducers`, and `test_robot_simulator`. They construct their own grids. They do not need the Moving AI files.

## Manuscript measurements

`scripts/run_manuscript_benchmarks.sh` reruns the fully revealed tables. It uses the leaf side and group size in the map table above. The k = 25 pass times H-TRIP and k times BFS on every map, and also pairwise A* on M1, M2, and M3. The speedup pass times H-TRIP and k times BFS at k = 2, 4, 16, 32, 64, and 256. The update pass times the opening cells in the evaluation: same-leaf openings of 1 and 16 cells on M1 and M4, and a 64-cell multi-leaf opening on M1 and M3. Each timed value is the median of ten repeats after three warm-ups. Viewpoint samples are the homogeneous sector sample. Every matrix is checked against BFS.

The files are:

- `campaigns/manuscript/query_k25.csv` and `full_build.csv` for the k = 25 table, including index size and full-build time
- `campaigns/manuscript/query_vs_k.csv` for the speedup table
- `campaigns/manuscript/update_latency.csv` for the opening table

`scripts/run_closed_loop.sh` reruns the exploration table. Each row has sixteen robots, fifty frontier viewpoints, 5 cm cells, speed 2 m/s, at most 500 simulation steps, and stops when free-cell coverage reaches 95 percent. The lidar and leaf side are the six rows of that table. The script writes `bfs_compare.csv` and `mission_summary.json` under `campaigns/closed_loop/` and does not record video.

How a CSV column becomes the number printed in the evaluation is written in [`scripts/README.md`](scripts/README.md).

A one-map smoke run, Boston at k = 25 with one repeat:

```bash
scripts/run_manuscript_benchmarks.sh --quick
```

## Exploration videos

Five sixteen-robot recordings are in [`videos/`](videos/README.md). Each clip shows the fleet opening free space and replanning on the occupancy revealed so far. The player is the smaller upload. The link under it is the full-size file in this repository.

**video-01.** Paris_0_512, leaf side 64, 256 frontier viewpoints, lidar 5 m. Coarse tiles on the radial boulevards, with a large distance matrix at every planning cycle.

https://github.com/user-attachments/assets/7d9afdfe-2027-43fe-8a18-52bc48654f19

Full-size file: [video-01.mp4](videos/video-01.mp4)

**video-02.** The same Paris map and 256 frontiers, with leaf side 32. The boulevards cross more tiles, and an opening updates a shorter ancestor chain.

https://github.com/user-attachments/assets/34d3664b-374b-4582-be09-8a81a305544f

Full-size file: [video-02.mp4](videos/video-02.mp4)

**video-03.** Paris_0_512, leaf side 64, 50 frontiers, lidar 5 m. The hierarchy stays coarse and the query matrix is the smaller frontier set from the exploration table.

https://github.com/user-attachments/assets/163187bd-fe27-4eed-8909-619dd706fbdb

Full-size file: [video-03.mp4](videos/video-03.mp4)

**video-04.** NewYork_0_512, the orthogonal Manhattan grid, leaf side 8, 50 frontiers, lidar 5 m. Long avenues sit in a fine tiling. This street map is a demonstration; it is not one of the six evaluation grids under `datasets/`.

https://github.com/user-attachments/assets/26f35ff0-8f19-4272-865a-baa02246fd18

Full-size file: [video-04.mp4](videos/video-04.mp4)

**video-05.** maze512-16-0, corridors 16 cells wide, leaf side 16, lidar 4 m at 90 degrees, 50 frontiers. This is the narrow-corridor exploration row for paper map M4. The grid is `datasets/movingai/maze/maps/maze512-16-0.map`.

https://github.com/user-attachments/assets/c9761dc5-ac40-4d06-8c8c-3786b615c79b

Full-size file: [video-05.mp4](videos/video-05.mp4)
