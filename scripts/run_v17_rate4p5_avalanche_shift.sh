#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER=clean_arrhenius_taylor_ddd_v17.py
ROOT=results/v17_CrossS_m9p25_floor0p50_rate4p5_FAC0p25_LeffFeed_avalancheShift
LOG="$ROOT/live_sweep_status.log"

T=1100
RATE=4.5
DT=1e-8
TARGET_STRAIN=0.006
MAX_FREE_DX=0.5

FAC=0.25
CROSS_S=-9.25
PEIERLS_S=0.0
CROSS_FLOOR=0.50
PEIERLS_FLOOR=0.0

RHO_LIST="1e12 3e12 1e13 3e13 1e14 3e14 1e15 3e15"

rm -rf "$ROOT"
mkdir -p "$ROOT"

{
echo "============================================================"
echo "v17 one-decade higher-rate avalanche-shift test"
echo "ROOT=$ROOT"
echo "DRIVER=$DRIVER"
echo "T=$T K"
echo "RATE=$RATE s^-1"
echo "DT=$DT s"
echo "CrossS=$CROSS_S"
echo "cross_floor_frac=$CROSS_FLOOR"
echo "FAC=$FAC"
echo "tau_local_cap_mode=none"
echo "tau_local_length_mode=feed"
echo "RHO_LIST=$RHO_LIST"
echo "============================================================"
} | tee -a "$LOG"

for RHO in $RHO_LIST; do
  OUTDIR="$ROOT/bs_on/T${T}_rho${RHO}"
  mkdir -p "$OUTDIR"

  echo | tee -a "$LOG"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] START rho=$RHO" | tee -a "$LOG"
  echo "OUTDIR=$OUTDIR" | tee -a "$LOG"

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
    --expfit-cross-entropy-kB "$CROSS_S" \
    --expfit-cross-floor-frac "$CROSS_FLOOR" \
    --expfit-peierls-floor-frac "$PEIERLS_FLOOR" \
    --tau-local-cap-mode none \
    --tau-local-length-mode feed \
    --tau-local-L-eff-reduced 1.0

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE rho=$RHO" | tee -a "$LOG"

  echo "--- run_summary.txt ---" | tee -a "$LOG"
  cat "$OUTDIR/run_summary.txt" | tee -a "$LOG"

  echo "--- burst diagnostic for rho=$RHO ---" | tee -a "$LOG"
  python3 analyze_depin_burst_statistics.py \
    --run-dir "$OUTDIR" \
    --cluster-gap-steps 1 \
    --active-plastic-ratio 1.0 \
    --stress-drop-threshold-MPa 0.0 \
    --n-boot 200 \
    --show-table | tee -a "$LOG"

  echo "--- current sweep table ---" | tee -a "$LOG"
  python3 analyze_v6_results.py \
    --root "$ROOT" \
    --show-table | tee -a "$LOG"

  python3 analyze_depin_burst_statistics.py \
    --root "$ROOT" \
    --cluster-gap-steps 1 \
    --active-plastic-ratio 1.0 \
    --stress-drop-threshold-MPa 0.0 \
    --n-boot 200 \
    --show-table | tee -a "$LOG"
done

echo | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "FINAL COMBINED TABLE" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

python3 analyze_v6_results.py \
  --root "$ROOT" \
  --show-table | tee -a "$LOG"

python3 analyze_depin_burst_statistics.py \
  --root "$ROOT" \
  --cluster-gap-steps 1 \
  --active-plastic-ratio 1.0 \
  --stress-drop-threshold-MPa 0.0 \
  --n-boot 500 \
  --show-table | tee -a "$LOG"

echo "DONE: $ROOT" | tee -a "$LOG"
