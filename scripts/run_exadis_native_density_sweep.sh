#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/opendis311/bin/python}"
BUILD_DIR="${BUILD_DIR:-core/exadis/build-audit-llvm}"
CONFIG="${ARRHENIUS_CONFIG:-configs/exadis_arrhenius_physics_forward.json}"
OUTROOT="${OUTROOT:-results/exadis_native_arrhenius_density_sweep}"
MAX_STRAIN="${MAX_STRAIN:-2.0e-6}"
TEMPERATURE="${TEMPERATURE:-300}"
STRAIN_RATE="${STRAIN_RATE:-1e3}"
DENSITY_FACTORS="${DENSITY_FACTORS:-0.5 1.0 2.0}"

mkdir -p "$OUTROOT"
validation_args=()
for factor in $DENSITY_FACTORS; do
  label="rho${factor}"
  outdir="$OUTROOT/$label"
  if [[ "${REUSE_COMPLETED:-0}" == "1" && -f "$outdir/audit_enabled/final_summary.json" ]]; then
    echo "reusing completed native density case: $label"
  else
    PYTHON_BIN="$PYTHON_BIN" \
    BUILD_DIR="$BUILD_DIR" \
    ARRHENIUS_CONFIG="$CONFIG" \
    ARRHENIUS_TEMPERATURE_K="$TEMPERATURE" \
    STRAIN_RATE="$STRAIN_RATE" \
    DENSITY_FACTOR="$factor" \
    MAX_STRAIN="$MAX_STRAIN" \
    MINIMUM_STEPS="${MINIMUM_STEPS:-5}" \
    AUDIT_STRIDE="${AUDIT_STRIDE:-5}" \
    AUDIT_ENABLED_ONLY=1 \
    PRINT_FREQ=1 \
    OUTDIR="$outdir" \
      bash scripts/run_exadis_native_arrhenius_A3.sh
  fi
  validation_args+=(--case "${factor},${outdir}")
done

"$PYTHON_BIN" exadis_calibration/native_density_trend_validation.py \
  "${validation_args[@]}" \
  --output "$OUTROOT/physics_density_validation.json"
