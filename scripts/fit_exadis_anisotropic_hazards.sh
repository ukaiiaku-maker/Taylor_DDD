#!/usr/bin/env bash
set -euo pipefail

# Fit equivalent anisotropic Arrhenius hazards against an ExaDiS audit JSONL.
# This script assumes the event audit was produced by
# scripts/run_exadis_binding_event_audit.sh or by a future native audit path.

ROOT=${ROOT:-results/exadis_binding_event_audit/audit_enabled}
AUDIT_JSONL=${AUDIT_JSONL:-$ROOT/event_audit.jsonl}
STRESS_STRAIN=${STRESS_STRAIN:-$ROOT/stress_strain_dens.dat}
OUTDIR=${OUTDIR:-results/exadis_anisotropic_hazard_fit}
PYTHON_BIN=${PYTHON_BIN:-python3}
TEMPERATURE_K=${TEMPERATURE_K:-900}
BURGERS_M=${BURGERS_M:-2.55e-10}
STRAIN_RATE_S=${STRAIN_RATE_S:-1.0e3}

if [[ ! -f "$AUDIT_JSONL" ]]; then
  echo "error: audit JSONL not found: $AUDIT_JSONL" >&2
  echo "run scripts/run_exadis_binding_event_audit.sh first, or set AUDIT_JSONL=/path/to/event_audit.jsonl" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

echo "=== fitting anisotropic equivalent hazards ==="
echo "audit:       $AUDIT_JSONL"
echo "stress data: $STRESS_STRAIN"
echo "outdir:      $OUTDIR"
echo "python:      $PYTHON_BIN"
echo "T:           $TEMPERATURE_K K"
echo "strain rate: $STRAIN_RATE_S s^-1"

args=(
  exadis_calibration/anisotropic_hazard_fit.py
  "$AUDIT_JSONL"
  --outdir "$OUTDIR"
  --temperature-K "$TEMPERATURE_K"
  --burgers-m "$BURGERS_M"
  --strain-rate-s "$STRAIN_RATE_S"
)

if [[ -f "$STRESS_STRAIN" ]]; then
  args+=(--stress-strain-dens "$STRESS_STRAIN")
fi

"$PYTHON_BIN" "${args[@]}" | tee "$OUTDIR/fit_stdout.log"

echo "=== outputs ==="
ls -lh "$OUTDIR"
