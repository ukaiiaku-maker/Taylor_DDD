#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/opendis311/bin/python}"
BUILD_DIR="${BUILD_DIR:-core/exadis/build-audit-llvm}"
CONFIG="${ARRHENIUS_CONFIG:-configs/exadis_arrhenius_physics_forward.json}"
OUTROOT="${OUTROOT:-results/exadis_native_arrhenius_physics_sweep}"
MAX_STRAIN="${MAX_STRAIN:-2.2e-7}"
TEMPERATURES="${TEMPERATURES:-300 500 700 900 1100}"
STRAIN_RATES="${STRAIN_RATES:-1e1 1e2 1e3 1e4}"

mkdir -p "$OUTROOT"
validation_args=()
for temperature in $TEMPERATURES; do
  for rate in $STRAIN_RATES; do
    label="T${temperature}_rate${rate}"
    outdir="$OUTROOT/$label"
    if [[ "${REUSE_COMPLETED:-0}" == "1" && -f "$outdir/audit_enabled/final_summary.json" ]]; then
      echo "reusing completed native physics case: $label"
    else
      PYTHON_BIN="$PYTHON_BIN" \
      BUILD_DIR="$BUILD_DIR" \
      ARRHENIUS_CONFIG="$CONFIG" \
      ARRHENIUS_TEMPERATURE_K="$temperature" \
      STRAIN_RATE="$rate" \
      MAX_STRAIN="$MAX_STRAIN" \
      MINIMUM_STEPS="${MINIMUM_STEPS:-5}" \
      AUDIT_STRIDE="${AUDIT_STRIDE:-5}" \
      AUDIT_ENABLED_ONLY=1 \
      PRINT_FREQ=1 \
      OUTDIR="$outdir" \
        bash scripts/run_exadis_native_arrhenius_A3.sh
    fi
    validation_args+=(--case "${temperature},${rate},${outdir}")
  done
done

"$PYTHON_BIN" exadis_calibration/native_physics_trend_validation.py \
  "${validation_args[@]}" \
  --output "$OUTROOT/physics_trend_validation.json"
