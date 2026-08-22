#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/opendis311/bin/python}"
FIT_PYTHON_BIN="${FIT_PYTHON_BIN:-/Users/sdillon/Taylor_DDD/.venv-opendis/bin/python}"
BUILD_DIR="${BUILD_DIR:-core/exadis/build-audit-llvm}"
OUTDIR="${OUTDIR:-results/exadis_native_calibration_campaign}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_NUM_THREADS

rates=(1e1 1e2 1e3 1e4)
strains=(2.2e-9 2.2e-8 2.2e-7 2.2e-6)
stock_data="core/exadis/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data"

for index in "${!rates[@]}"; do
  rate="${rates[$index]}"
  strain="${strains[$index]}"
  run_dir="$OUTDIR/rate_${rate}"
  if [[ -f "$run_dir/audit_enabled/event_audit.jsonl" ]]; then
    continue
  fi
  if [[ "$rate" == "1e3" && -f results/exadis_native_candidate_smoke/audit_enabled/event_audit.jsonl ]]; then
    continue
  fi
  PYTHON_BIN="$PYTHON_BIN" BUILD_DIR="$BUILD_DIR" OUTDIR="$run_dir" \
    MAX_STRAIN="$strain" STRAIN_RATE="$rate" CROSS_SLIP=1 \
    AUDIT_ENABLED_ONLY=1 REQUIRE_CANDIDATE_LABELS="" \
    bash scripts/run_exadis_native_audit.sh
done

# Produce a second native initial configuration by evolving the stock network
# one deterministic step, then use that complete ParaDiS configuration as an
# out-of-sample network for every rate.  The fitter holds state_2 out wholesale.
seed_dir="$OUTDIR/state_2_seed"
state_2_data="$seed_dir/audit_enabled/config.1.data"
if [[ ! -f "$state_2_data" ]]; then
  PYTHON_BIN="$PYTHON_BIN" BUILD_DIR="$BUILD_DIR" OUTDIR="$seed_dir" \
    EXADIS_DATA="$stock_data" MAX_STRAIN=1e-7 STRAIN_RATE=1e3 CROSS_SLIP=1 \
    OUTPUT_FREQUENCY=1 AUDIT_ENABLED_ONLY=1 REQUIRE_CANDIDATE_LABELS="" \
    bash scripts/run_exadis_native_audit.sh
fi
if [[ ! -f "$state_2_data" ]]; then
  echo "second native initial state was not written: $state_2_data" >&2
  exit 2
fi

for index in "${!rates[@]}"; do
  rate="${rates[$index]}"
  strain="${strains[$index]}"
  run_dir="$OUTDIR/state_2/rate_${rate}"
  if [[ -f "$run_dir/audit_enabled/event_audit.jsonl" ]]; then
    continue
  fi
  PYTHON_BIN="$PYTHON_BIN" BUILD_DIR="$BUILD_DIR" OUTDIR="$run_dir" \
    EXADIS_DATA="$state_2_data" MAX_STRAIN="$strain" STRAIN_RATE="$rate" \
    CROSS_SLIP=1 AUDIT_ENABLED_ONLY=1 REQUIRE_CANDIDATE_LABELS="" \
    bash scripts/run_exadis_native_audit.sh
done

manifest="$OUTDIR/campaign_manifest.json"
mkdir -p "$OUTDIR"
"$FIT_PYTHON_BIN" - "$ROOT_DIR" "$OUTDIR" "$manifest" <<'PY'
import json
import sys
from pathlib import Path

root, outdir, output = map(Path, sys.argv[1:])
runs = []
for state in ("state_1", "state_2"):
    for rate in ("1e1", "1e2", "1e3", "1e4"):
        if state == "state_1" and rate == "1e3":
            audit = root / "results/exadis_native_candidate_smoke/audit_enabled/event_audit.jsonl"
        elif state == "state_1":
            audit = outdir / f"rate_{rate}/audit_enabled/event_audit.jsonl"
        else:
            audit = outdir / state / f"rate_{rate}/audit_enabled/event_audit.jsonl"
        runs.append({
            "run_id": f"{state}:rate_{rate}",
            "initial_state_id": state,
            "strain_rate_s": float(rate),
            "audit_jsonl": str(audit),
        })
payload = {
    "temperatures_K": [300, 500, 700, 900, 1100],
    "adaptive_event_integration": True,
    "runs": runs,
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

PYTHONPATH="$ROOT_DIR/exadis_calibration:${PYTHONPATH:-}" \
  "$FIT_PYTHON_BIN" exadis_calibration/native_arrhenius_campaign_fit.py \
  "$manifest" --outdir "$OUTDIR/fit"

PYTHONPATH="$ROOT_DIR/exadis_calibration:${PYTHONPATH:-}" \
  "$FIT_PYTHON_BIN" exadis_calibration/native_mobility_config_gate.py \
  "$OUTDIR/fit/mobility_fit_observed_vs_predicted.csv" \
  configs/exadis_arrhenius_mobility_gate_passed.json \
  --output "$OUTDIR/fit/production_config_gate_summary.json" \
  --adaptive-event-integration

set +e
PYTHONPATH="$ROOT_DIR/exadis_calibration:${PYTHONPATH:-}" \
  "$FIT_PYTHON_BIN" exadis_calibration/native_discrete_hazard_fit.py \
  "$manifest" --outdir results/exadis_native_discrete_hazard_fit
discrete_status=$?
set -e
if [[ $discrete_status -ne 0 && $discrete_status -ne 2 ]]; then
  exit "$discrete_status"
fi
echo "discrete mechanism fit gate exit status: $discrete_status (2 means safely rejected)"
