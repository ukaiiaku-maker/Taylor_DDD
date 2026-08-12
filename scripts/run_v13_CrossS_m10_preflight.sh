#!/usr/bin/env bash
set -euo pipefail

# Quick preflight for CrossS=-10 at the predicted peak window.
# Run this before the full overnight sweep if time permits.

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v13.py}"
ROOT="results/v13_CrossS_m10_preflight_rate0p3"
rm -rf "$ROOT"

RATE=0.30
DT=$(python3 - <<PY
rate=float("$RATE")
print("{:.12g}".format(1e-7/rate))
PY
)

for RHO in 3e15 5e15 1e16; do
  OUTDIR="$ROOT/bs_on/T1100_rho${RHO}"
  mkdir -p "$OUTDIR"
  python3 "$DRIVER" \
    --outdir "$OUTDIR" \
    --temperature-K 1100 \
    --strain-rate "$RATE" \
    --target-strain 0.006 \
    --dt "$DT" \
    --forest-rho-m2 "$RHO" \
    --backstress-mobility on \
    --backstress-com-projection external_drive \
    --capture-mode swept_crossing \
    --snap-swept-capture-to-obstacle \
    --crossing-drive-mode force_work \
    --cross-force-scale-mode line_tension \
    --cross-force-scale-factor 1.0 \
    --max-free-dx-reduced 0.5 \
    --plastic-strain-source actual \
    --expfit-peierls-entropy-kB 0.0 \
    --expfit-cross-entropy-kB -10.0
done

python3 analyze_v6_results.py --root "$ROOT" --show-table --make-plots
