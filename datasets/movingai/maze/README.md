# Maze maps

The two files in this folder are Moving AI maze grids used as M3 and M4 in the H-TRIP evaluation. They are already unpacked. H-TRIP reads them as occupancy grids: `.`, `G`, and `S` are passable, and `@`, `O`, `T`, and `W` are blocked. A passable cell is a vertex, and each 4-neighbor step between passable cells is an edge of length 1.

| Paper id | File | Grid | Passable cells | Leaf b, group g |
|----------|------|------|----------------|-----------------|
| M3 | maze512-1-0.map | 512 by 512, corridors one cell wide | 131,071 | 16, 2 |
| M4 | maze512-16-0.map | 512 by 512, corridors 16 cells wide | 246,016 | 16, 2 |

M3 is a tree of one-cell corridors, so a shortest path between distant cells is a long unique route. M4 uses the same generator with corridors 16 cells wide. Both files belong to the maze benchmark set published with the Moving AI grid collection. The street maps used beside them live in [`../street/maps/`](../street/maps/). The six-map list is in [`../README.md`](../README.md).

Sixteen robots exploring M4 with a 4 m lidar, leaf side 16, and 50 frontiers are recorded as t03:

<video src="../../../videos/t03-maze512-16-k50-r4.mp4" controls width="720"></video>

The corridor is 16 cells wide. The fleet opens a short stretch of maze at each scan, and the planner keeps 50 frontier viewpoints on the occupancy revealed so far. The other four recordings are in [`../../../videos/README.md`](../../../videos/README.md).

## Citation

If you use these mazes, cite the benchmark paper that published the collection:

Sturtevant, N. R. Benchmarks for Grid-Based Pathfinding. IEEE Transactions on Computational Intelligence and AI in Games, 4(2), 144–148, 2012. https://doi.org/10.1109/TCIAIG.2012.2197681

```bibtex
@article{sturtevant2012benchmarks,
  author  = {Sturtevant, Nathan R.},
  title   = {Benchmarks for Grid-Based Pathfinding},
  journal = {IEEE Transactions on Computational Intelligence and AI in Games},
  volume  = {4},
  number  = {2},
  pages   = {144--148},
  year    = {2012},
  doi     = {10.1109/TCIAIG.2012.2197681}
}
```

The files are redistributed from that collection: https://movingai.com/benchmarks/grids.html

The grids remain subject to the terms of the Moving AI Lab collection.
