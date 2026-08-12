#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER=clean_arrhenius_taylor_ddd_v17.py
ROOT=results/v17_CrossS_m9p25_floor0p50_rate0p45_FAC0p25_LeffFeed_rhoSweep
LOG="$ROOT/live_sweep_status.log"

T=1100
RATE=0.45
DT=1e-7
TARGET_STRAIN=0.006
MAX_FREE_DX=0.5
FAC=0.25
CROSS_S=-9.25
PEIERLS_S=0.0
CROSS_FLOOR=0.50
PEIERLS_FLOOR=0.0
TAU_LOCAL_CAP_MODE=none
TAU_LOCAL_LENGTH_MODE=feed
TAU_LOCAL_L_EFF_MIN=1.0

RHO_LIST="1e13 1e14 3e14 1e15 3e15 1e16 3e16"

rm -rf "$ROOT"
mkdir -p "$ROOT"

{
  echo "============================================================"
  echo "v17 uncapped physical-L_eff rho sweep"
  echo "ROOT=$ROOT"
  echo "DRIVER=$DRIVER"
  echo "T=$T K"
  echo "RATE=$RATE s^-1"
  echo "DT=$DT s"
  echo "CrossS=$CROSS_S"
  echo "cross_floor_frac=$CROSS_FLOOR"
  echo "FAC=$FAC"
  echo "tau_local_cap_mode=$TAU_LOCAL_CAP_MODE"
  echo "tau_local_length_mode=$TAU_LOCAL_LENGTH_MODE"
  echo "tau_local_L_eff_min=$TAU_LOCAL_L_EFF_MIN b"
  echo "RHO_LIST=$RHO_LIST"
  echo "============================================================"
} | tee -a "$LOG"

for RHO in $RHO_LIST; do
  OUTDIR="$ROOT/bs_on/T${T}_rho${RHO}"
  mkdir -p "$OUTDIR"

  {
    echo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START rho=$RHO"
    echo "OUTDIR=$OUTDIR"
  } | tee -a "$LOG"

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
    --tau-local-cap-mode "$TAU_LOCAL_CAP_MODE" \
    --tau-local-length-mode "$TAU_LOCAL_LENGTH_MODE" \
    --tau-local-L-eff-reduced "$TAU_LOCAL_L_EFF_MIN"

  {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE rho=$RHO"
    echo "--- run_summary.txt ---"
    cat "$OUTDIR/run_summary.txt"
    echo "--- cap/avalanche diagnostic for rho=$RHO ---"
  } | tee -a "$LOG"

  if [[ -f analyze_depin_burst_statistics.py ]]; then
    python3 analyze_depin_burst_statistics.py --run-dir "$OUTDIR" --show-table | tee -a "$LOG" || true
  fi

  {
    echo "--- current sweep table ---"
  } | tee -a "$LOG"

  python3 analyze_v6_results.py --root "$ROOT" --show-table | tee -a "$LOG" || true
  if [[ -f analyze_depin_burst_statistics.py ]]; then
    python3 analyze_depin_burst_statistics.py --root "$ROOT" --show-table | tee -a "$LOG" || true
  fi
done

{
  echo
  echo "============================================================"
  echo "FINAL COMBINED TABLE"
  echo "============================================================"
} | tee -a "$LOG"

python3 analyze_v6_results.py --root "$ROOT" --show-table | tee -a "$LOG" || true
if [[ -f analyze_depin_burst_statistics.py ]]; then
  python3 analyze_depin_burst_statistics.py --root "$ROOT" --show-table | tee -a "$LOG" || true
fi

python3 - <<'PY' | tee -a "$LOG"
from pathlib import Path
import pandas as pd
ROOT = Path("results/v17_CrossS_m9p25_floor0p50_rate0p45_FAC0p25_LeffFeed_rhoSweep")
print()
print("V17 SOFTENING / CAP / BURST CHECK")
print("=================================")
s1 = ROOT / "analysis" / "summary_v6.csv"
if s1.exists():
    df = pd.read_csv(s1).sort_values("rho_m2")
    cols = [
        "rho_m2", "tau_tail_abs_median_MPa", "epsp_over_epstotal_final",
        "n_crossed_total_final", "n_live_pins_tail_median",
        "tau_local_median_tail_MPa", "frac_tau_local_capped_tail", "flags",
    ]
    cols = [c for c in cols if c in df.columns]
    print(df[cols].to_string(index=False))
    tau = df["tau_tail_abs_median_MPa"].to_numpy()
    rho = df["rho_m2"].to_numpy()
    if len(tau):
        imax = int(tau.argmax())
        print(f"Peak tau_tail_abs_median_MPa={tau[imax]:.6g} at rho={rho[imax]:.6e}")
        if imax < len(tau)-1:
            print(f"Drop from peak to highest rho={tau[imax]-tau[-1]:.6g} MPa")
else:
    print(f"Missing {s1}")

s2 = ROOT / "depin_burst_summary_all.csv"
if s2.exists():
    av = pd.read_csv(s2).sort_values("rho_m2")
    keep = [
        "rho_m2", "total_depin", "largest_burst_depin", "largest_burst_fraction",
        "p90_burst_depin", "p99_burst_depin", "max_stress_drop_negsum_MPa",
        "frac_depin_events_over_tau_cap_reference", "frac_depin_events_clipped", "null_largest_fraction_p_ge_obs",
        "null_fano_p_ge_obs", "burst_distribution_suggestive", "interpretation",
    ]
    keep = [c for c in keep if c in av.columns]
    print()
    print("Burst statistics:")
    print(av[keep].to_string(index=False))
PY

echo "DONE: $ROOT" | tee -a "$LOG"
