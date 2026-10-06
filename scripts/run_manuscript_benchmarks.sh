#!/usr/bin/env bash
# Reproduce the fully revealed tables in the H-TRIP evaluation.
#
# k = 25, sector viewpoints, three warm-ups, ten repeats, on all six maps.
# M1, M2, and M3 also time pairwise A* at that k.
# Speedup sweep: k in {2, 4, 16, 32, 64, 256}, H-TRIP and k x BFS.
# Update cells are the (map, pattern, m) rows of the opening table,
# ten withheld locations each.
#
# Pass --quick to time only Boston_0_256 at k = 25 with one repeat.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN="${HTRIP_BIN:-$ROOT/build/dev}"
OUT="${HTRIP_OUT:-$ROOT/campaigns/manuscript}"
DATA="$ROOT/datasets/movingai"

if [[ ! -x "$BIN/bench_query_vs_k" || ! -x "$BIN/bench_update_latency" ]]; then
    echo "Build the project first: cmake --preset dev && cmake --build --preset dev" >&2
    exit 1
fi

QUICK=0
if [[ "${1:-}" == "--quick" ]]; then
    QUICK=1
fi

WARM=3
REP=10
if [[ "$QUICK" == 1 ]]; then
    WARM=1
    REP=1
fi

mkdir -p "$OUT"
export OMP_NUM_THREADS=1

# id|map|b|g|astar_at_k25 (1 or 0)
maps=(
    "Boston_0_256|$DATA/street/maps/Boston_0_256.map|16|2|1"
    "Boston_0_1024|$DATA/street/maps/Boston_0_1024.map|32|2|1"
    "maze512-1-0|$DATA/maze/maps/maze512-1-0.map|16|2|1"
    "maze512-16-0|$DATA/maze/maps/maze512-16-0.map|16|2|0"
    "Paris_0_512|$DATA/street/maps/Paris_0_512.map|8|2|0"
    "Berlin_0_1024|$DATA/street/maps/Berlin_0_1024.map|8|2|0"
)

run_k25() {
    local id="$1" map="$2" b="$3" g="$4" with_astar="$5"
    local methods="HTRIP,KBFS"
    if [[ "$with_astar" == "1" ]]; then
        methods="HTRIP,KBFS,ASTAR"
    fi
    echo "=== k=25 $id b=$b g=$g methods=$methods ==="
    "$BIN/bench_query_vs_k" \
        --map "$map" --map-id "$id" --b "$b" --g "$g" \
        --paper-e5 --methods "$methods" \
        --warmup "$WARM" --repeats "$REP" \
        --csv "$OUT/query_k25.csv" \
        --repeats-csv "$OUT/query_k25_repeats.csv" \
        --build-csv "$OUT/full_build.csv" \
        --rss-csv "$OUT/rss_k25.csv"
}

run_sweep() {
    local id="$1" map="$2" b="$3" g="$4"
    echo "=== sweep $id b=$b g=$g k=2,4,16,32,64,256 ==="
    "$BIN/bench_query_vs_k" \
        --map "$map" --map-id "$id" --b "$b" --g "$g" \
        --geometry sector --k 2,4,16,32,64,256 --methods HTRIP,KBFS \
        --warmup "$WARM" --repeats "$REP" --no-build --no-rss \
        --csv "$OUT/query_vs_k.csv" \
        --repeats-csv "$OUT/query_vs_k_repeats.csv"
}

run_update() {
    local id="$1" map="$2" b="$3" pattern="$4" mlist="$5"
    echo "=== update $id $pattern m=$mlist ==="
    "$BIN/bench_update_latency" \
        --map "$map" --map-id "$id" --b "$b" --g 2 \
        --pattern "$pattern" --m "$mlist" \
        --locations 10 --warmup "$WARM" --repeats "$REP" \
        --csv "$OUT/update_latency.csv"
}

for spec in "${maps[@]}"; do
    IFS='|' read -r id map b g with_astar <<< "$spec"
    if [[ "$QUICK" == 1 && "$id" != "Boston_0_256" ]]; then
        continue
    fi
    run_k25 "$id" "$map" "$b" "$g" "$with_astar"
    if [[ "$QUICK" == 1 ]]; then
        continue
    fi
    run_sweep "$id" "$map" "$b" "$g"
done

if [[ "$QUICK" == 0 ]]; then
    run_update "Boston_0_256" "$DATA/street/maps/Boston_0_256.map" 16 same-leaf 1,16
    run_update "Boston_0_256" "$DATA/street/maps/Boston_0_256.map" 16 multi-leaf 64
    run_update "maze512-1-0" "$DATA/maze/maps/maze512-1-0.map" 16 multi-leaf 64
    run_update "maze512-16-0" "$DATA/maze/maps/maze512-16-0.map" 16 same-leaf 1,16
fi

echo "k = 25:     $OUT/query_k25.csv"
echo "Full build: $OUT/full_build.csv"
if [[ "$QUICK" == 0 ]]; then
    echo "k sweep:    $OUT/query_vs_k.csv"
    echo "Updates:    $OUT/update_latency.csv"
fi
