# Moving AI maps

These six grids are the maps used in the H-TRIP evaluation. They are already unpacked:

```
datasets/movingai/<set>/maps/<file>.map
```

| Paper id | Set | File | Grid | Leaf b, group g |
|----------|-----|------|------|-----------------|
| M1 | street | Boston_0_256.map | Boston streets, 256 by 256, 47,768 passable cells | 16, 2 |
| M2 | street | Boston_0_1024.map | Boston streets, 1024 by 1024, 796,896 passable cells | 32, 2 |
| M3 | maze | maze512-1-0.map | 512 by 512 maze, corridors one cell wide, 131,071 passable cells | 16, 2 |
| M4 | maze | maze512-16-0.map | 512 by 512 maze, corridors 16 cells wide, 246,016 passable cells | 16, 2 |
| M5 | street | Paris_0_512.map | Paris streets, 512 by 512, 196,567 passable cells | 8, 2 |
| M6 | street | Berlin_0_1024.map | Berlin streets, 1024 by 1024, 794,748 passable cells | 8, 2 |

Each file is an occupancy grid. H-TRIP treats passable cells as free space and builds a shortest-path index on the 4-neighbor grid. Distances are stored as 16-bit tropical values.

The two maze files, and the Moving AI citation for that maze set, are in [`maze/README.md`](maze/README.md). The manuscript benchmark loads these six files. `scripts/run_manuscript_benchmarks.sh` times fully revealed many-to-many queries and index updates. `scripts/run_closed_loop.sh` times the same query while sixteen robots explore. How each CSV column becomes a published number is in [`../../scripts/README.md`](../../scripts/README.md).

## Citation

If you use these grids, cite the benchmark paper that published the collection:

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

The street maps in the collection were contributed by Konstantin Yakovlev and Anton Andreychuk. The grids remain subject to the terms of the Moving AI Lab collection.
