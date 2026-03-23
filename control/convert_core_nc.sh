#!/bin/bash
#
# Batch convert MATCH CCM history files to netCDF using ccm2nc.
#
# Usage: convert_core_nc.sh [run_dir]
#   run_dir defaults to ~/Data/MATCH/2026/sarb4_hourly_202602_core

set -euo pipefail

CCM2NC="/home/fillmore/EarthSystem/MATCH/utils/ccm2nc"
export LD_LIBRARY_PATH="/home/fillmore/miniconda3/envs/sarb/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

RUN_DIR="${1:-$HOME/Data/MATCH/2026/sarb4_hourly_202602_core}"

if [[ ! -d "$RUN_DIR" ]]; then
    echo "Error: directory $RUN_DIR not found" >&2
    exit 1
fi

cd "$RUN_DIR"

count=0
skip=0
fail=0

for hfile in h[0-9][0-9][0-9][0-9]; do
    ncfile="${hfile}.nc"
    if [[ -f "$ncfile" ]]; then
        skip=$((skip + 1))
        continue
    fi
    if "$CCM2NC" -C "$hfile" "$ncfile" 2>&1; then
        count=$((count + 1))
    else
        echo "FAILED: $hfile" >&2
        fail=$((fail + 1))
    fi
done

echo "Done: converted=$count skipped=$skip failed=$fail"
