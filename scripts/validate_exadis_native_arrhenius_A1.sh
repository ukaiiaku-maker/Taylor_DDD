#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/opendis311/bin/python}"
STOCK_DIR="${STOCK_DIR:-results/exadis_native_audit_post_kernel}"
ARRHENIUS_DIR="${ARRHENIUS_DIR:-results/exadis_native_arrhenius_A1_long}"
OUTPUT="${OUTPUT:-results/exadis_native_arrhenius_A1_long/stage_validation.json}"

"$PYTHON_BIN" exadis_calibration/native_stage_validation.py \
  --stock "$STOCK_DIR" \
  --arrhenius "$ARRHENIUS_DIR" \
  --output "$OUTPUT" \
  --adaptive-event-integration
