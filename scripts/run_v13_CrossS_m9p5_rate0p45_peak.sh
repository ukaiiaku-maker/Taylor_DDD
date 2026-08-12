#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER=clean_arrhenius_taylor_ddd_v13.py
ROOT=results/v13_CrossS_m9p5_FAC1_T1100_rate0p45_peak
rm -rf "$ROOT"

T=1100
RATE=0.45
DT=1e-7
TARGET_STRAIN=0.006
MAX_FREE_DX=0.5
FAC=1.0
CROSS_S=-9.5
PEIERLS_S=0.0

# Dense around expected peak ~7e15–1e16, sparse outside.
RHO_LIST="1e14 3e14 1e15 3e15 5e15 7e15 1e16 1.5e16 2e16 3e16 1e17 3e17"

echo "============================================================"
echo "v13 peak search"
echo "CrossS=$CROSS_S"
echo "PeierlsS=$PEIERLS_S"
echo "RATE=$RATE s^-1"
echo "DT=$DT"
echo "ROOT=$ROOT"
echo "RHO_LIST=$RHO_LIST"
echo "============================================================"

for RHO in $RHO_LIST; do
  OUTDIR="$ROOT/bs_on/T${T}_rho${RHO}"
  mkdir -p "$OUTDIR"

  echo "[$(date '+%H:%M:%S')] rho=$RHO -> $OUTDIR"

  python3 "$DRIVER" \
    --outdir "$OUTDIR" \
    --temperature-K "$T" \
    --strain-rate "$RATE" \
    --target-strain "$TARGET_STRAIN" \
    --dt "$DT" \
    --forest-rho-m2 "$RHO" \
    --backstress-mobility on \
    --backstress-com-projection external_drive \
    --capture-mode swept_crossing \
    --snap-swept-capture-to-obstacle \
    --crossing-drive-mode force_work \
    --cross-force-scale-mode line_tension \
    --cross-force-scale-factor "$FAC" \
    --max-free-dx-reduced "$MAX_FREE_DX" \
    --plastic-strain-source actual \
    --expfit-peierls-entropy-kB "$PEIERLS_S" \
    --expfit-cross-entropy-kB "$CROSS_S"

  echo "   ok"
done

python3 analyze_v6_results.py \
  --root "$ROOT" \
  --show-table \
  --make-plots

echo
echo "DONE:"
echo "  $ROOT/analysis/summary_v6.csv"
