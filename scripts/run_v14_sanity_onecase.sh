#!/usr/bin/env bash
set -euo pipefail
cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true
DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v14.py}"
ROOT="${ROOT:-results/v14_sanity_CrossS_m9p5_rate0p45_rho5e15}"
rm -rf "$ROOT"
OUTDIR="$ROOT/bs_on/T1100_rho5e15"
mkdir -p "$OUTDIR"
python3 "$DRIVER" \
  --outdir "$OUTDIR" \
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
  --expfit-peierls-entropy-kB 0.0 \
  --expfit-cross-entropy-kB -9.5
python3 analyze_v6_results.py --root "$ROOT" --show-table || true
cat "$OUTDIR/run_summary.txt" || true
