#!/usr/bin/env bash
set -euo pipefail
cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER=${DRIVER:-clean_arrhenius_taylor_ddd_v16.py}
ROOT=results/v16_avalanche_sanity_CrossS_m9p25_floor0p50_rate0p45_rho5e15_FAC0p5
OUTDIR="$ROOT/bs_on/T1100_rho5e15"
rm -rf "$ROOT"
mkdir -p "$OUTDIR"

python3 "$DRIVER" \
  --outdir "$OUTDIR" \
  --temperature-K 1100 \
  --strain-rate 0.45 \
  --target-strain 0.0002 \
  --dt 1e-7 \
  --forest-rho-m2 5e15 \
  --backstress-mobility on \
  --backstress-com-projection external_drive \
  --capture-mode swept_crossing \
  --snap-swept-capture-to-obstacle \
  --crossing-drive-mode force_work \
  --cross-force-scale-mode line_tension \
  --cross-force-scale-factor 0.5 \
  --max-free-dx-reduced 0.5 \
  --plastic-strain-source actual \
  --expfit-peierls-entropy-kB 0.0 \
  --expfit-cross-entropy-kB -9.25 \
  --expfit-cross-floor-frac 0.50 \
  --expfit-peierls-floor-frac 0.0 \
  --avalanche-quiet-steps 2 \
  --avalanche-window-steps 100 \
  --avalanche-min-events 2

cat "$OUTDIR/run_summary.txt"
echo
cat "$OUTDIR/avalanche_summary.txt"
echo
head -20 "$OUTDIR/avalanche_events.csv"
