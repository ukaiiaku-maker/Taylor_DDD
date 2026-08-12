#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

# Make sure the runner actually forwards EXTRA_ARGS if you use it elsewhere.
# This v13 diagnostic does not require EXTRA_ARGS; actual swept strain is default.
DRIVER=clean_arrhenius_taylor_ddd_v13.py
RHO_LIST="3e16"
BS_LIST="on"
MAX_FREE_DX=0.5

for FAC in 1.0 0.5 0.25; do
  TAG=${FAC/./p}
  ROOT=results/v13_actualSwept_crossForceFactor_${TAG}_T1100_rho3e16

  echo "=== Running CROSS_FORCE_SCALE_FACTOR=${FAC}, ROOT=${ROOT} ==="
  DRIVER="$DRIVER" \
  ROOT="$ROOT" \
  RHO_LIST="$RHO_LIST" \
  BS_LIST="$BS_LIST" \
  MAX_FREE_DX="$MAX_FREE_DX" \
  CROSS_FORCE_SCALE_FACTOR="$FAC" \
  bash run_v11_density_forceWork.sh

  python3 analyze_v6_results.py \
    --root "$ROOT" \
    --show-table
done

echo
echo "Check v13 actual swept-strain columns:"
for FAC in 1.0 0.5 0.25; do
  TAG=${FAC/./p}
  ROOT=results/v13_actualSwept_crossForceFactor_${TAG}_T1100_rho3e16
  HIST=$(find "$ROOT" -name "single_glider_history.csv" | head -1)
  echo "--- $ROOT"
  python3 - <<PY
import pandas as pd
hist = "$HIST"
df = pd.read_csv(hist)
last = df.iloc[-1]
cols = [
    "eps_total", "eps_plastic", "eps_plastic_actual",
    "eps_plastic_book_free_glide", "eps_plastic_book_total",
    "n_crossed_total", "n_depin",
]
print(df[cols].tail(1).to_string(index=False))
print("sum d_eps_p_actual =", df["d_eps_p_actual"].sum())
print("sum d_eps_p_book_free_glide =", df["d_eps_p_book_free_glide"].sum())
print("sum n_depin =", df["n_depin"].sum())
PY
done
