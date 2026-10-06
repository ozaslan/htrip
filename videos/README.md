# Exploration videos

These five clips are sixteen-robot runs on Moving AI grids. Each robot carries a lidar, the fleet keeps a set of frontier viewpoints, and the planner reads a dense distance matrix on the occupancy revealed so far. Each recording is 0:25.

## t17 — Paris, leaf 64, 256 frontiers

Paris_0_512, radial boulevards. Sixteen robots, lidar 5 m at 90 degrees with 30 rays, leaf side 64, group size 2, and 256 frontier viewpoints. The coarse leaves cover a long stretch of boulevard, and each planning cycle fills a large matrix.

https://github.com/ozaslan/htrip/blob/main/videos/t17-paris512-k256-b64.mp4

## t16 — Paris, leaf 32, 256 frontiers

The same Paris map, lidar, and 256 frontiers, with leaf side 32. Boulevards are split across more tiles, so an opening updates a shorter ancestor chain while the frontier set stays dense.

https://github.com/ozaslan/htrip/blob/main/videos/t16-paris512-k256-b32.mp4

## t13 — Paris, leaf 64, 50 frontiers

Paris_0_512 again, leaf side 64, with the frontier count used in the exploration table: 50 viewpoints and a 5 m lidar. The hierarchy is coarse and the query matrix is the smaller one.

https://github.com/ozaslan/htrip/blob/main/videos/t13-paris512-k50-b64.mp4

## t08 — New York, leaf 8, 50 frontiers

NewYork_0_512, an orthogonal Manhattan grid from the same Moving AI street collection. Sixteen robots, lidar 5 m at 90 degrees with 30 rays, leaf side 8, and 50 frontiers. Avenues run for many cells inside a fine tiling. This grid is a demonstration recording; it is not one of the six evaluation maps shipped under `datasets/`.

https://github.com/ozaslan/htrip/blob/main/videos/t08-newyork512-k50-b8.mp4

## t03 — maze512-16-0, 4 m lidar

maze512-16-0, corridors 16 cells wide. Sixteen robots, lidar 4 m at 90 degrees with 30 rays, leaf side 16, and 50 frontiers. This is the narrow-corridor exploration row for paper map M4. The map file is `datasets/movingai/maze/maps/maze512-16-0.map`.

https://github.com/ozaslan/htrip/blob/main/videos/t03-maze512-16-k50-r4.mp4
