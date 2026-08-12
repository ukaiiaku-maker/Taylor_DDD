#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v13.py}"
T="${T:-1100}"
TARGET_STRAIN="${TARGET_STRAIN:-0.006}"
MAX_FREE_DX="${MAX_FREE_DX:-0.5}"
FAC="${FAC:-1.0}"

# Densities chosen to bracket the predicted observable peak:
# rate 0.5: peak likely near ~1e16
# rate 1.0: peak likely near ~1.5e16-3e16
RHO_LIST="${RHO_LIST:-1e15 3e15 6e15 1e16 1.5e16 2e16 3e16 6e16}"

for RATE in 0.5 1.0; do
  if [[ "$RATE" == "0.5" ]]; then
    DT="2e-7"
    RTAG="0p5"
  elif [[ "$RATE" == "1.0" ]]; then
    DT="1e-7"
    RTAG="1p0"
  else
    echo "Unknown RATE=$RATE" >&2
    exit 2
  fi

  ROOT="results/v13_actualSwept_FAC1_rate${RTAG}_T1100_peakSweep"
  echo
  echo "============================================================"
  echo "Running RATE=$RATE s^-1, DT=$DT"
  echo "ROOT=$ROOT"
  echo "RHO_LIST=$RHO_LIST"
  echo "============================================================"

  rm -rf "$ROOT"

  DRIVER="$DRIVER" \
  ROOT="$ROOT" \
  RHO_LIST="$RHO_LIST" \
  BS_LIST="on" \
  STRAIN_RATE="$RATE" \
  TARGET_STRAIN="$TARGET_STRAIN" \
  DT="$DT" \
  MAX_FREE_DX="$MAX_FREE_DX" \
  CROSS_FORCE_SCALE_FACTOR="$FAC" \
  bash run_v11_density_forceWork.sh

  python3 analyze_v6_results.py \
    --root "$ROOT" \
    --show-table \
    --make-plots
done

python3 - <<'PY'
from pathlib import Path
import pandas as pd

roots = [
    Path("results/v13_actualSwept_FAC1_rate0p5_T1100_peakSweep"),
    Path("results/v13_actualSwept_FAC1_rate1p0_T1100_peakSweep"),
]

frames = []
for root in roots:
    f = root / "analysis" / "summary_v6.csv"
    if f.exists():
        df = pd.read_csv(f)
        df["root"] = str(root)
        frames.append(df)

if frames:
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["T_K", "rho_m2"])
    out.to_csv("results/v13_actualSwept_FAC1_rate_peakSweep_combined_summary.csv", index=False)
    print("\nCombined summary:")
    cols = [
        "root", "T_K", "rho_m2",
        "tau_tail_median_MPa",
        "epsp_over_epstotal_final",
        "n_crossed_total_final",
        "n_live_pins_tail_median",
        "tau_local_median_tail_MPa",
    ]
    print(out[[c for c in cols if c in out.columns]].to_string(index=False))
    print("\nWrote results/v13_actualSwept_FAC1_rate_peakSweep_combined_summary.csv")
else:
    print("No summary files found.")
PY
