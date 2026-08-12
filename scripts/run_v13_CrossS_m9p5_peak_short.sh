#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v13.py}"
OUTROOT="${OUTROOT:-results/v13_CrossS_m9p5_FAC1_T1100_peakShort}"

T=1100
TARGET_STRAIN=0.006
MAX_FREE_DX=0.5
FAC=1.0
CROSS_S=-9.5
PEIERLS_S=0.0

# Broad enough to see the peak/turnover if it appears, but avoid 1e18 for this first corrected run.
RHO_LIST="1e14 3e14 1e15 3e15 5e15 7e15 1e16 3e16 1e17 3e17"

for RATE in 0.30 0.60; do
  if [[ "$RATE" == "0.30" ]]; then
    RTAG=0p30
    DT=1.666666667e-7
  elif [[ "$RATE" == "0.60" ]]; then
    RTAG=0p60
    DT=8.333333333e-8
  else
    echo "Unknown RATE=$RATE" >&2
    exit 2
  fi

  ROOT="$OUTROOT/rate${RTAG}"
  rm -rf "$ROOT"

  echo
  echo "============================================================"
  echo "RATE=$RATE s^-1"
  echo "DT=$DT"
  echo "CrossS=$CROSS_S"
  echo "ROOT=$ROOT"
  echo "RHO_LIST=$RHO_LIST"
  echo "============================================================"

  for RHO in $RHO_LIST; do
    OUTDIR="$ROOT/bs_on/T${T}_rho${RHO}"
    mkdir -p "$OUTDIR"

    echo "[$(date '+%H:%M:%S')] rate=$RATE rho=$RHO -> $OUTDIR"

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
    --make-plots || true
done

echo
echo "DONE. Results under:"
echo "  $OUTROOT"
