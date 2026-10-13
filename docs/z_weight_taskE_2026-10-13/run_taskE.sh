#!/usr/bin/env bash
# Task E (CLOUD_PROMPT_Z_weight_taskE_2026-10-13.md): six Model VII runs, one process per core, OMP_NUM_THREADS=1.
# Usage: run_taskE.sh <ogs binary> <deck dir> <output root>
set -u
OGS=$1; DECK=$2; OUT=$3
export OMP_NUM_THREADS=1
unset OGS_ASM_THREADS   # unset = 1 assembly thread (BaseLib/OgsAsmThreads.cpp)
cd "$DECK"
for name in ms33_modelVII_freeswelling ms33_modelVII_freeswelling_level ms33_modelVII_freeswelling_levelz \
            ms33_modelVII_freeswelling_ladder550 ms33_modelVII_freeswelling_ladder550_level ms33_modelVII_freeswelling_ladder550_levelz; do
  mkdir -p "$OUT/$name"
  (
    t0=$(date +%s)
    "$OGS" "$name.prj" -o "$OUT/$name" > "$OUT/$name/$name.log" 2>&1
    rc=$?
    t1=$(date +%s)
    echo "rc=$rc wall_s=$((t1-t0))" > "$OUT/$name/wall.txt"
  ) &
done
wait
