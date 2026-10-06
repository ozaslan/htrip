#!/usr/bin/env bash
# Exploration rows from the H-TRIP evaluation.
#
# Sixteen robots, fifty frontier viewpoints, 5 cm cells, 2 m/s, at most
# 500 simulation steps, stop at 95 percent free-cell coverage. Paired
# H-TRIP and k x BFS. No video and no coverage stills.
#
# Logs go to campaigns/closed_loop/<name>/.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN="${HTRIP_BIN:-$ROOT/build/dev/run_multi_robot_sim}"
OUT="${HTRIP_OUT:-$ROOT/campaigns/closed_loop}"
DATA="$ROOT/datasets/movingai"

if [[ ! -x "$BIN" ]]; then
    echo "Build the project first: cmake --preset dev && cmake --build --preset dev" >&2
    exit 1
fi

mkdir -p "$OUT"
export OMP_NUM_THREADS=1

run_one() {
    local name="$1" map="$2" b="$3" range="$4" span="$5" rays="$6"
    local dir="$OUT/$name"
    mkdir -p "$dir"
    echo "=== $name b=$b lidar=${range}m/${span}deg/${rays}rays ==="
    "$BIN" \
        --map "$map" \
        --robots 16 \
        --frontiers 50 \
        --range "$range" \
        --span "$span" \
        --rays "$rays" \
        --v-max 2.0 \
        --b "$b" --g 2 \
        --cell_size_cm 5 \
        --dt 0.05 \
        --bfs 1 \
        --astar 0 \
        --async 1 \
        --steps 500 \
        --stop-coverage 0.95 \
        --stills 0 \
        --log-dir "$dir"
}

run_one boston256_b16_r6 "$DATA/street/maps/Boston_0_256.map" 16 6.0 180 45
run_one maze16_b16_r4 "$DATA/maze/maps/maze512-16-0.map" 16 4.0 90 30
run_one maze16_b16_r8 "$DATA/maze/maps/maze512-16-0.map" 16 8.0 180 60
run_one paris512_b8_r5 "$DATA/street/maps/Paris_0_512.map" 8 5.0 90 30
run_one paris512_b32_r5 "$DATA/street/maps/Paris_0_512.map" 32 5.0 90 30
run_one berlin1024_b8_r5 "$DATA/street/maps/Berlin_0_1024.map" 8 5.0 90 30

echo "Closed-loop logs: $OUT"
