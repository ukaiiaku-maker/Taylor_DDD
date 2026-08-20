#!/usr/bin/env bash
set -euo pipefail

# Binding-level stock ExaDiS audit runner.
#
# This runs the stock FCC Cu strain-hardening module stack twice through the
# Python stepping driver:
#   1. audit disabled
#   2. audit enabled
# and then compares final strain, stress, density, node count, and segment count.
#
# It does not connect Arrhenius hazards.  It is an instrumentation-only test.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"
EXADIS_ROOT="${EXADIS_ROOT:-core/exadis}"
EXADIS_DATA="${EXADIS_DATA:-$EXADIS_ROOT/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data}"
OUTDIR="${OUTDIR:-results/exadis_binding_event_audit}"
MAX_STRAIN="${MAX_STRAIN:-1.0e-5}"
AUDIT_STRIDE="${AUDIT_STRIDE:-1}"
PRINT_FREQ="${PRINT_FREQ:-1}"
WRITE_FREQ="${WRITE_FREQ:-0}"
TOLERANCE="${TOLERANCE:-1.0e-10}"
CROSS_SLIP_FLAG="${CROSS_SLIP_FLAG:-0}"

args=(
  --exadis-root "$EXADIS_ROOT"
  --exadis-data "$EXADIS_DATA"
  --outdir "$OUTDIR"
  --max-strain "$MAX_STRAIN"
  --audit-stride "$AUDIT_STRIDE"
  --print-freq "$PRINT_FREQ"
  --write-freq "$WRITE_FREQ"
  --tolerance "$TOLERANCE"
  --mode both
)

if [[ "$CROSS_SLIP_FLAG" == "1" ]]; then
  args+=(--cross-slip)
fi

"$PYTHON_BIN" exadis_audit/binding_event_audit.py "${args[@]}"

cat <<EOF

Wrote:
  $OUTDIR/audit_disabled/final_summary.json
  $OUTDIR/audit_enabled/final_summary.json
  $OUTDIR/audit_enabled/event_audit.jsonl
  $OUTDIR/audit_invariance_summary.json
EOF
