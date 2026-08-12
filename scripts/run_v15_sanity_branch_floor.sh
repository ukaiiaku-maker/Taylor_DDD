#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

ROOT=results/v15_sanity_CrossS_m9p25_floor0p50_rate0p45_rho5e15
rm -rf "$ROOT"

python3 clean_arrhenius_taylor_ddd_v15.py \
  --outdir "$ROOT/bs_on/T1100_rho5e15" \
  --temperature-K 1100 \
  --strain-rate 0.45 \
  --target-strain 0.006 \
  --dt 1e-7 \
  --forest-rho-m2 5e15 \
  --backstress-mobility on \
  --backstress-com-projection external_drive \
  --capture-mode swept_crossing \
  --snap-swept-capture-to-obstacle \
  --crossing-drive-mode force_work \
  --cross-force-scale-mode line_tension \
  --cross-force-scale-factor 1.0 \
  --max-free-dx-reduced 0.5 \
  --plastic-strain-source actual \
  --expfit-cross-entropy-kB -9.25 \
  --expfit-peierls-entropy-kB 0.0 \
  --expfit-cross-floor-frac 0.50 \
  --expfit-peierls-floor-frac 0.0 \
  --preflight-rdt-warn 0.2 \
  --preflight-rdt-stop 1.0

cat "$ROOT/bs_on/T1100_rho5e15/preflight_diagnostics.txt"
echo
cat "$ROOT/bs_on/T1100_rho5e15/run_summary.txt"
