#!/usr/bin/env bash
set -euo pipefail

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER=clean_arrhenius_taylor_ddd_v15.py
ROOT=results/v15_CrossS_m9p25_floor0p50_rate0p45_FAC0p25_rhoSweep
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

RHO_LIST="1e13 1e14 3e14 1e15 3e15 1e16 3e16"

rm -rf "$ROOT"
mkdir -p "$ROOT"

echo "============================================================" | tee -a "$LOG"
echo "v15 FAC=0.25 rho sweep with avalanche/cap diagnostics" | tee -a "$LOG"
echo "ROOT=$ROOT" | tee -a "$LOG"
echo "T=$T K" | tee -a "$LOG"
echo "RATE=$RATE s^-1" | tee -a "$LOG"
echo "DT=$DT s" | tee -a "$LOG"
echo "CrossS=$CROSS_S" | tee -a "$LOG"
echo "cross_floor_frac=$CROSS_FLOOR" | tee -a "$LOG"
echo "FAC=$FAC" | tee -a "$LOG"
echo "RHO_LIST=$RHO_LIST" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

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
    --expfit-peierls-floor-frac "$PEIERLS_FLOOR"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE rho=$RHO" | tee -a "$LOG"

  echo "--- run_summary.txt ---" | tee -a "$LOG"
  cat "$OUTDIR/run_summary.txt" | tee -a "$LOG"

  echo "--- cap/avalanche diagnostic for rho=$RHO ---" | tee -a "$LOG"
  python3 analyze_cap_avalanche_diagnostics.py \
    --run-dir "$OUTDIR" \
    --show-table | tee -a "$LOG"

  echo "--- current sweep table ---" | tee -a "$LOG"
  python3 analyze_v6_results.py \
    --root "$ROOT" \
    --show-table | tee -a "$LOG"

  python3 analyze_cap_avalanche_diagnostics.py \
    --root "$ROOT" \
    --show-table | tee -a "$LOG"

done

echo | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "FINAL COMBINED TABLE" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

python3 analyze_v6_results.py \
  --root "$ROOT" \
  --show-table | tee -a "$LOG"

python3 analyze_cap_avalanche_diagnostics.py \
  --root "$ROOT" \
  --show-table | tee -a "$LOG"

python3 - <<'PY' | tee -a "$LOG"
from pathlib import Path
import pandas as pd

ROOT = Path("results/v15_CrossS_m9p25_floor0p50_rate0p45_FAC0p25_rhoSweep")
s1 = ROOT / "analysis" / "summary_v6.csv"
s2 = ROOT / "cap_avalanche_summary_all.csv"

print()
print("SOFTENING / AVALANCHE CHECK")
print("===========================")

if not s1.exists():
    print(f"Missing {s1}")
    raise SystemExit(0)

df = pd.read_csv(s1)
df = df.sort_values("rho_m2")

cols = [
    "rho_m2",
    "tau_tail_abs_median_MPa",
    "epsp_over_epstotal_final",
    "n_crossed_total_final",
    "n_live_pins_tail_median",
    "tau_local_median_tail_MPa",
    "frac_tau_local_capped_tail",
    "flags",
]
cols = [c for c in cols if c in df.columns]
print(df[cols].to_string(index=False))

tau = df["tau_tail_abs_median_MPa"].to_numpy()
rho = df["rho_m2"].to_numpy()
if len(tau) >= 2:
    imax = int(tau.argmax())
    print()
    print(f"Peak tau_tail_abs_median_MPa = {tau[imax]:.3f} at rho = {rho[imax]:.3e}")
    if imax < len(tau) - 1:
        drop = tau[imax] - tau[-1]
        print(f"Drop from peak to highest-rho point = {drop:.3f} MPa")
    else:
        print("No post-peak point yet; no softening conclusion.")

if s2.exists():
    av = pd.read_csv(s2)
    keep = [
        "rho_m2",
        "largest_event_avalanche_fraction",
        "largest_stress_drop_negsum_MPa",
        "event_count_window_fano",
        "top_1pct_steps_event_fraction",
        "history_tail_frac_tau_local_capped_median",
        "frac_depin_events_at_tau_cap",
        "cap_interpretation",
        "avalanche_like",
    ]
    keep = [c for c in keep if c in av.columns]
    print()
    print("Avalanche/cap table:")
    print(av[keep].sort_values("rho_m2").to_string(index=False))
PY

echo "DONE: $ROOT" | tee -a "$LOG"
