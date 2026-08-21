#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Two steps are required: the first starts at zero applied stress, while the
# second provides nonzero FCC_0 force/velocity rows for schema validation.
OUTDIR="${OUTDIR:-results/exadis_native_candidate_smoke}" \
MAX_STRAIN="${MAX_STRAIN:-2.2e-7}" \
CROSS_SLIP=1 \
REQUIRE_CANDIDATE_LABELS="cross_slip,collision,topology_split" \
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
bash "$ROOT_DIR/scripts/run_exadis_native_audit.sh"
