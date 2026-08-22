#!/usr/bin/env bash
set -euo pipefail

# Audit launcher for the stock ExaDiS FCC strain-hardening example.
#
# This script intentionally uses the adapter in two stages:
#   1. dry-run: always available; writes mechanism config and analytical peak table
#   2. stock ExaDiS run: optional; requires pyexadis and the stock example data file
#
# The native Arrhenius hooks are not linked into ExaDiS here.  This launcher is
# for auditing stress/force-work scales and producing a reproducible local run
# directory before native mobility/topology/cross-slip/collision integration.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"
OUTDIR="${OUTDIR:-results/exadis_arrhenius_stock_audit}"
EXADIS_ROOT="${EXADIS_ROOT:-core/exadis}"
EXADIS_EXAMPLE_DIR="$EXADIS_ROOT/examples/22_fcc_Cu_15um_1e3"
EXADIS_DATA="${EXADIS_DATA:-$EXADIS_EXAMPLE_DIR/180chains_16.10e.data}"
MAX_STRAIN="${MAX_STRAIN:-1e-4}"
STRAIN_RATE="${STRAIN_RATE:-1e3}"

mkdir -p "$OUTDIR"

echo "Using Python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import sys
print(sys.executable)
try:
    import numpy as np
    print('numpy', np.__version__)
except Exception as exc:
    raise SystemExit(f'numpy import failed: {exc}')
PY

echo "=== Stage 1: dry-run analytical/audit config ==="
"$PYTHON_BIN" scripts/exadis_arrhenius_strain_hardening_adapter.py \
  --dry-run \
  --outdir "$OUTDIR" \
  --strain-rate "$STRAIN_RATE" \
  --max-strain "$MAX_STRAIN"

if [[ ! -f "$EXADIS_DATA" ]]; then
  cat <<EOF

Stock ExaDiS data file was not found:
  $EXADIS_DATA

Dry-run artifacts were written under:
  $OUTDIR

To run the stock ExaDiS audit, clone or symlink ExaDiS so that the example data is available, for example:
  mkdir -p core
  git clone https://github.com/LLNL/exadis core/exadis

Then rebuild ExaDiS with Python bindings and rerun this script with PYTHONPATH pointing to the ExaDiS python directory.
EOF
  exit 0
fi

echo "=== Stage 2: stock ExaDiS strain-hardening audit run ==="
(
  cd "$EXADIS_EXAMPLE_DIR"
  PYTHONPATH="$ROOT_DIR:$ROOT_DIR/$EXADIS_ROOT/python:${PYTHONPATH:-}" \
  "$PYTHON_BIN" "$ROOT_DIR/scripts/exadis_arrhenius_strain_hardening_adapter.py" \
    --run-stock-exadis \
    --outdir "$ROOT_DIR/$OUTDIR" \
    --exadis-data "$ROOT_DIR/$EXADIS_DATA" \
    --strain-rate "$STRAIN_RATE" \
    --max-strain "$MAX_STRAIN"
)

echo "wrote audit outputs to $OUTDIR"
