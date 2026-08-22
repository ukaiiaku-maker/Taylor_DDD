#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/opendis311/bin/python}"
BUILD_DIR="${BUILD_DIR:-core/exadis/build-audit-llvm}"
CONFIG="${ARRHENIUS_CONFIG:-configs/exadis_arrhenius_physics_forward.json}"
OUTROOT="${OUTROOT:-results/exadis_native_arrhenius_density_sweep}"
MAX_STRAIN="${MAX_STRAIN:-1.0e-5}"
TEMPERATURE="${TEMPERATURE:-300}"
STRAIN_RATE="${STRAIN_RATE:-1e3}"
DENSITY_FACTORS="${DENSITY_FACTORS:-0.25 0.5 1.0 2.0 4.0 8.0}"
INTERACTION_PARAMETER_SET="${INTERACTION_PARAMETER_SET:-taylor_scaling_test_barrier}"

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
    INTERACTION_PARAMETER_SET="$INTERACTION_PARAMETER_SET" \
    AUDIT_STRIDE="${AUDIT_STRIDE:-5}" \
    AUDIT_ENABLED_ONLY=1 \
    PRINT_FREQ=1 \
    OUTDIR="$outdir" \
      bash scripts/run_exadis_native_arrhenius_A3.sh
  fi
  validation_args+=(--case "${factor},${outdir}")
done

status=0
"$PYTHON_BIN" exadis_calibration/native_density_trend_validation.py \
  "${validation_args[@]}" \
  --output "$OUTROOT/physics_density_validation.json" || status=$?

diagnostic_cases=()
audit_files=()
for factor in $DENSITY_FACTORS; do
  diagnostic_cases+=("$OUTROOT/rho${factor}")
  audit_files+=("$OUTROOT/rho${factor}/audit_enabled/event_audit.jsonl")
done
"$PYTHON_BIN" analysis/native_taylor_density_diagnostics.py \
  "${diagnostic_cases[@]}" \
  --output "$OUTROOT/taylor_density_diagnostics.json" || status=1

"$PYTHON_BIN" analysis/check_taylor_force_work_double_counting.py \
  "${audit_files[@]}" \
  --output "$OUTROOT/double_counting_validation.json" || status=1

exit "$status"
