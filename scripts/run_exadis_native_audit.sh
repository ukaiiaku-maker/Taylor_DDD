#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"
EXADIS_ROOT="${EXADIS_ROOT:-core/exadis}"
BUILD_DIR="${BUILD_DIR:-$EXADIS_ROOT/build-audit}"
EXADIS_DATA="${EXADIS_DATA:-$EXADIS_ROOT/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data}"
OUTDIR="${OUTDIR:-results/exadis_native_audit}"
MAX_STRAIN="${MAX_STRAIN:-1.0e-5}"
STRAIN_RATE="${STRAIN_RATE:-1.0e3}"
AUDIT_STRIDE="${AUDIT_STRIDE:-1}"
PRINT_FREQ="${PRINT_FREQ:-1}"
TOLERANCE="${TOLERANCE:-1.0e-12}"
CROSS_SLIP="${CROSS_SLIP:-0}"
REQUIRE_CANDIDATE_LABELS="${REQUIRE_CANDIDATE_LABELS:-}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_NUM_THREADS

if [[ "$BUILD_DIR" = /* ]]; then
  PYEXADIS_DIR="$BUILD_DIR/python"
else
  PYEXADIS_DIR="$ROOT_DIR/$BUILD_DIR/python"
fi
if [[ "$EXADIS_ROOT" = /* ]]; then
  EXADIS_PYTHON_DIR="$EXADIS_ROOT/python"
else
  EXADIS_PYTHON_DIR="$ROOT_DIR/$EXADIS_ROOT/python"
fi

args=(
  --exadis-root "$EXADIS_ROOT"
  --exadis-data "$EXADIS_DATA"
  --outdir "$OUTDIR"
  --max-strain "$MAX_STRAIN"
  --strain-rate "$STRAIN_RATE"
  --audit-stride "$AUDIT_STRIDE"
  --print-freq "$PRINT_FREQ"
  --tolerance "$TOLERANCE"
)
if [[ "$CROSS_SLIP" == "1" ]]; then
  args+=(--cross-slip)
fi
if [[ -n "$REQUIRE_CANDIDATE_LABELS" ]]; then
  IFS=',' read -r -a required_mechanisms <<< "$REQUIRE_CANDIDATE_LABELS"
  for mechanism in "${required_mechanisms[@]}"; do
    args+=(--require-candidate-labels "$mechanism")
  done
fi

PYTHONPATH="$PYEXADIS_DIR:$EXADIS_PYTHON_DIR:${PYTHONPATH:-}" \
  "$PYTHON_BIN" exadis_audit/native_event_audit.py "${args[@]}"

echo "native invariance report: $OUTDIR/invariance_summary.json"
