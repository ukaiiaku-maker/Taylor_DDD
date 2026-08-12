#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:?Usage: monitor_v14_run.sh <results_root> [sleep_seconds]}"
SLEEP="${2:-60}"
TARGET="${TARGET_STRAIN:-0.006}"

echo "Monitoring ROOT=$ROOT every ${SLEEP}s"
echo "Ctrl-C stops monitoring only; it does not stop the simulation."
while true; do
  echo
  echo "================ $(date '+%Y-%m-%d %H:%M:%S') ================"
  HIST=$(find "$ROOT" -name single_glider_history.csv -type f -print0 2>/dev/null | xargs -0 ls -t 2>/dev/null | head -n 1 || true)
  if [[ -z "${HIST:-}" ]]; then
    echo "No history files found yet."
    sleep "$SLEEP"
    continue
  fi
  echo "Latest history: $HIST"
  python3 - <<PY
import pandas as pd, math
from pathlib import Path
hist=Path("$HIST")
target=float("$TARGET")
try:
    df=pd.read_csv(hist, low_memory=False)
except Exception as e:
    print(f"Could not read {hist}: {e}"); raise SystemExit(0)
for c in ["step","eps_total","eps_plastic","tau_MPa","tau_before_step_MPa","tau_after_step_MPa","d_tau_step_MPa","d_eps_p","d_eps_total","d_eps_p_over_d_eps_total","n_crossed_total","n_live_pins","n_depin","tau_local_median_MPa","frac_tau_local_capped","crossing_rate_max_s","step_walltime_s"]:
    if c in df.columns:
        df[c]=pd.to_numeric(df[c], errors="coerce")
if "eps_total" not in df:
    print("No eps_total column yet."); raise SystemExit(0)
df=df[df["eps_total"].between(0, target*1.05)]
df=df.dropna(subset=["step","eps_total"])
if len(df)==0:
    print("No valid numeric rows yet."); raise SystemExit(0)
last=df.iloc[-1]
tail=df.tail(max(20, len(df)//10))
print(f"valid_rows={len(df)}")
print(f"step={last.get('step', float('nan')):.0f}")
print(f"eps_total={last.get('eps_total', float('nan')):.6g} ({100*last.get('eps_total',0)/target:.1f}% of target {target})")
print(f"eps_plastic={last.get('eps_plastic', float('nan')):.6g}, epsp/eps={last.get('eps_plastic',0)/max(last.get('eps_total',float('nan')),1e-300):.6g}")
print(f"tau_after_MPa={last.get('tau_MPa', float('nan')):.3f}, tau_before_MPa={last.get('tau_before_step_MPa', float('nan')):.3f}, d_tau_step={last.get('d_tau_step_MPa', float('nan')):.3g}")
print(f"d_eps_p/d_eps_total={last.get('d_eps_p_over_d_eps_total', float('nan')):.3g}, d_eps_p={last.get('d_eps_p', float('nan')):.3g}, d_eps_total={last.get('d_eps_total', float('nan')):.3g}")
print(f"n_crossed_total={last.get('n_crossed_total', float('nan')):.0f}, n_depin_step={last.get('n_depin', float('nan')):.0f}, n_live_pins={last.get('n_live_pins', float('nan')):.0f}")
print(f"tau_local_median_MPa={last.get('tau_local_median_MPa', float('nan')):.3f}, frac_capped={last.get('frac_tau_local_capped', float('nan')):.3f}")
if "d_tau_step_MPa" in tail:
    dec=(tail["d_tau_step_MPa"]<0).sum(); inc=(tail["d_tau_step_MPa"]>0).sum()
    print(f"tail d_tau: inc_steps={inc}, dec_steps={dec}, median={tail['d_tau_step_MPa'].median():.3g} MPa")
if "crossing_rate_max_s" in df:
    # If params json available, use dt; otherwise infer from consecutive time if present.
    dt=None
    p=hist.parent/'clean_arrhenius_params.json'
    if p.exists():
        try:
            import json
            dt=float(json.loads(p.read_text()).get('dt'))
        except Exception:
            dt=None
    if dt:
        print(f"max_Rdt_seen={(df['crossing_rate_max_s']*dt).max():.3g}, tail_Rdt_median={(tail['crossing_rate_max_s']*dt).median():.3g}")
if "step_walltime_s" in tail:
    print(f"recent step_walltime_s median={tail['step_walltime_s'].median():.4g}")
PY
  sleep "$SLEEP"
done
