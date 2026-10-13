#!/usr/bin/env bash
# Driver of the DSM element Jacobian FD verification
# (docs/dsm_jacobian_fd_verification_2026-10.md).
# Usage: Tests/Python/dsm_jacobian_fd_run.sh <build dir> <output dir>
# Runs the Google test RichardsMechanics.DSMElementJacobianFD with the CSV and
# matrix dumps enabled and writes the markdown tables to <output dir>/tables.md.
set -euo pipefail
build=${1:?build dir}
out=${2:?output dir}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out/mat"
OMP_NUM_THREADS=1 \
OGS_DSM_JACFD_CSV="$out/results.csv" \
OGS_DSM_JACFD_MATRICES="$out/mat" \
    "$build/bin/testrunner" --gtest_filter='RichardsMechanics.DSMElementJacobianFD' -l error \
    | tee "$out/testrunner.log"
python3 "$here/dsm_jacobian_fd_analysis.py" "$out" > "$out/tables.md"
echo "tables: $out/tables.md"
