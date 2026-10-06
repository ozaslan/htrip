# Manuscript benchmarks

`run_manuscript_benchmarks.sh` and `run_closed_loop.sh` rerun the measurements reported in the H-TRIP evaluation. Both scripts load the six unpacked maps under `datasets/movingai/`. Build first:

```bash
cmake --preset dev
cmake --build --preset dev
```

The query and update binary checks every H-TRIP matrix against BFS and stops if any entry differs. Distances are 16-bit values with unreachable stored as 60000. Min-plus kernels run with AVX2. The machine that produced the published medians was an Intel Core i7-8650U, so another CPU will not reprint those milliseconds. The map, the leaf size, the group size, the viewpoint count, the warm-up count, and the repeat count are the same.

## Fully revealed queries

Every ground-truth passable cell is marked observed free before the index is built. Viewpoints are a homogeneous sector sample. Each timed call is the median of ten repeats after three warm-ups. `OMP_NUM_THREADS` is 1.

| Paper id | File | b | g |
|----------|------|---|---|
| M1 | street/maps/Boston_0_256.map | 16 | 2 |
| M2 | street/maps/Boston_0_1024.map | 32 | 2 |
| M3 | maze/maps/maze512-1-0.map | 16 | 2 |
| M4 | maze/maps/maze512-16-0.map | 16 | 2 |
| M5 | street/maps/Paris_0_512.map | 8 | 2 |
| M6 | street/maps/Berlin_0_1024.map | 8 | 2 |

`scripts/run_manuscript_benchmarks.sh` writes four CSV files under `campaigns/manuscript/`.

`query_k25.csv` is the k = 25 table. One row is one method. `median_ms` for `HTRIP` and `KBFS` are the two query columns. On M1, M2, and M3 the file also contains `ASTAR`, which the text compares with `KBFS` at this same k. `full_build.csv` holds `build_ms` and `index_bytes` for that index. The published index size is `index_bytes` divided by 1024 squared, in mebibytes. The published break-even count is `build_ms` divided by the difference between the `KBFS` median and the `HTRIP` median at k = 25.

`query_vs_k.csv` is the speedup table. k is 2, 4, 16, 32, 64, and 256. Methods are `HTRIP` and `KBFS`. The published ratio is the `KBFS` median divided by the `HTRIP` median on the same map and the same k.

`update_latency.csv` is the opening table. The script runs only the cells that appear there:

| Map | Pattern | Cells opened, m |
|-----|---------|-----------------|
| M1 Boston_0_256 | same-leaf | 1 and 16 |
| M1 Boston_0_256 | multi-leaf | 64 |
| M3 maze512-1-0 | multi-leaf | 64 |
| M4 maze512-16-0 | same-leaf | 1 and 16 |

Each cell is repeated at ten withheld locations. `htrip_upd_median_ms` is the update. `kbfs_dt_median_ms` is a fresh k times BFS fill of the distance matrix after the opening. The published update number is the median of `htrip_upd_median_ms` across those ten locations, and the published refill number is the median of `kbfs_dt_median_ms`. The published update-plus-query number is that update median plus the k = 25 `HTRIP` median from `query_k25.csv` on the same map.

A smoke run times only Boston at k = 25, with one warm-up and one repeat, and skips the sweep and the update cells:

```bash
scripts/run_manuscript_benchmarks.sh --quick
```

## Exploration

`scripts/run_closed_loop.sh` runs the six exploration rows. Each row uses sixteen robots, fifty frontier viewpoints, 5 cm cells, speed 2 m/s, at most 500 simulation steps, and stops earlier when free-cell coverage reaches 95 percent. H-TRIP and k times BFS are timed on the occupancy revealed so far. The script does not record video or coverage stills.

| Paper id | Lidar | b | Log directory |
|----------|-------|---|---------------|
| M1 | 6 m, 180 degrees, 45 rays | 16 | campaigns/closed_loop/boston256_b16_r6 |
| M4 | 4 m, 90 degrees, 30 rays | 16 | campaigns/closed_loop/maze16_b16_r4 |
| M4 | 8 m, 180 degrees, 60 rays | 16 | campaigns/closed_loop/maze16_b16_r8 |
| M5 | 5 m, 90 degrees, 30 rays | 8 | campaigns/closed_loop/paris512_b8_r5 |
| M5 | 5 m, 90 degrees, 30 rays | 32 | campaigns/closed_loop/paris512_b32_r5 |
| M6 | 5 m, 90 degrees, 30 rays | 8 | campaigns/closed_loop/berlin1024_b8_r5 |

Group size is 2. Each directory contains `bfs_compare.csv` and `mission_summary.json`. The published query columns are the paired H-TRIP and k times BFS times in that comparison file, and the coverage column is the final free-cell coverage in the summary.
