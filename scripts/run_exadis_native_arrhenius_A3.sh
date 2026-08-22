#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export ARRHENIUS_MOBILITY="${ARRHENIUS_MOBILITY:-peierls}"
export ARRHENIUS_CONFIG="${ARRHENIUS_CONFIG:-configs/exadis_arrhenius_physics_forward.json}"
export ARRHENIUS_TEMPERATURE_K="${ARRHENIUS_TEMPERATURE_K:-900}"
export ARRHENIUS_TOPOLOGY="on"
export ARRHENIUS_CROSS_SLIP="on"
export ARRHENIUS_COLLISION="off"
export CROSS_SLIP="1"
export OUTDIR="${OUTDIR:-results/exadis_native_arrhenius_A3}"

bash scripts/run_exadis_native_audit.sh
