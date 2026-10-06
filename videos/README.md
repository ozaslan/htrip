# Exploration videos

These five clips are sixteen-robot runs on Moving AI grids. Each robot carries a lidar, the fleet keeps a set of frontier viewpoints, and the planner reads a dense distance matrix on the occupancy revealed so far. Each recording is 0:25. The player uses the smaller upload. The link under it is the full-size file stored in this folder.

## video-01 — Paris, leaf 64, 256 frontiers

Paris_0_512, radial boulevards. Sixteen robots, lidar 5 m at 90 degrees with 30 rays, leaf side 64, group size 2, and 256 frontier viewpoints. The coarse leaves cover a long stretch of boulevard, and each planning cycle fills a large matrix.

https://github.com/user-attachments/assets/7d9afdfe-2027-43fe-8a18-52bc48654f19

Full-size file: [video-01.mp4](video-01.mp4)

## video-02 — Paris, leaf 32, 256 frontiers

The same Paris map, lidar, and 256 frontiers, with leaf side 32. Boulevards are split across more tiles, so an opening updates a shorter ancestor chain while the frontier set stays dense.

https://github.com/user-attachments/assets/34d3664b-374b-4582-be09-8a81a305544f

Full-size file: [video-02.mp4](video-02.mp4)

## video-03 — Paris, leaf 64, 50 frontiers

Paris_0_512 again, leaf side 64, with the frontier count used in the exploration table: 50 viewpoints and a 5 m lidar. The hierarchy is coarse and the query matrix is the smaller one.

https://github.com/user-attachments/assets/163187bd-fe27-4eed-8909-619dd706fbdb

Full-size file: [video-03.mp4](video-03.mp4)

## video-04 — New York, leaf 8, 50 frontiers

NewYork_0_512, an orthogonal Manhattan grid from the same Moving AI street collection. Sixteen robots, lidar 5 m at 90 degrees with 30 rays, leaf side 8, and 50 frontiers. Avenues run for many cells inside a fine tiling. This grid is a demonstration recording; it is not one of the six evaluation maps shipped under `datasets/`.

https://github.com/user-attachments/assets/26f35ff0-8f19-4272-865a-baa02246fd18

Full-size file: [video-04.mp4](video-04.mp4)

## video-05 — maze512-16-0, 4 m lidar

maze512-16-0, corridors 16 cells wide. Sixteen robots, lidar 4 m at 90 degrees with 30 rays, leaf side 16, and 50 frontiers. This is the narrow-corridor exploration row for paper map M4. The map file is `datasets/movingai/maze/maps/maze512-16-0.map`.

https://github.com/user-attachments/assets/c9761dc5-ac40-4d06-8c8c-3786b615c79b

Full-size file: [video-05.mp4](video-05.mp4)
