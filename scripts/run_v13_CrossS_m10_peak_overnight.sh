#!/usr/bin/env bash
set -euo pipefail

# v13 overnight peak-search sweep with a modified forest-crossing barrier.
#
# Physics choice:
#   CrossS = -10 kB instead of -9 kB.
#   This increases the zero-force forest-crossing barrier by ~kT at 1100 K
#   (~0.095 eV), increasing the unloaded pin lifetime by e ~= 2.7.
#
# Goal:
#   Move the analytically predicted peak/turnover condition to rho ~5e15 m^-2
#   while keeping Peierls unchanged.  We sweep three rates around the predicted
#   ideal rate (~0.25-0.30 s^-1).
#
# This script calls the v13 driver directly, not run_v11_density_forceWork.sh,
# so the barrier arguments are guaranteed to be passed.

cd /Volumes/Data/Data/DDD/OpenDiS
source .venv-opendis/bin/activate 2>/dev/null || true

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v13.py}"
T="${T:-1100}"
TARGET_STRAIN="${TARGET_STRAIN:-0.006}"
MAX_FREE_DX="${MAX_FREE_DX:-0.5}"
FAC="${FAC:-1.0}"
CROSS_S="${CROSS_S:--10.0}"
PEIERLS_S="${PEIERLS_S:-0.0}"
OUTROOT="${OUTROOT:-results/v13_CrossS_m10_FAC1_T1100_peakOvernight}"

# Broad but not overly dense.  Keep the target region around 5e15 well sampled.
RHO_LIST="${RHO_LIST:-1e14 3e14 1e15 3e15 5e15 7e15 1e16 3e16 1e17 3e17 1e18}"

# Three rates around the analytically expected rate for rho_peak ~5e15.
RATE_LIST="${RATE_LIST:-0.15 0.30 0.60}"

mkdir -p "$OUTROOT"

echo "============================================================"
echo "v13 CrossS=$CROSS_S peak search"
echo "OUTROOT=$OUTROOT"
echo "DRIVER=$DRIVER"
echo "T=$T K"
echo "TARGET_STRAIN=$TARGET_STRAIN"
echo "FAC=$FAC"
echo "PEIERLS_S=$PEIERLS_S"
echo "CROSS_S=$CROSS_S"
echo "RHO_LIST=$RHO_LIST"
echo "RATE_LIST=$RATE_LIST"
echo "============================================================"

for RATE in $RATE_LIST; do
  RTAG="${RATE/./p}"
  ROOT="$OUTROOT/rate${RTAG}"
  mkdir -p "$ROOT"

  # Keep the strain increment per step equal to 1e-7, as in prior tests.
  DT=$(python3 - <<PY
rate = float("$RATE")
print("{:.12g}".format(1e-7 / rate))
PY
)

  echo
  echo "============================================================"
  echo "RATE=$RATE s^-1  DT=$DT  ROOT=$ROOT"
  echo "============================================================"

  for RHO in $RHO_LIST; do
    OUTDIR="$ROOT/bs_on/T${T}_rho${RHO}"
    echo "[$(date '+%H:%M:%S')] rate=$RATE rho=$RHO -> $OUTDIR"
    mkdir -p "$OUTDIR"

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

  # Analyze after each rate block, if analyzer is available.
  if [[ -f analyze_v6_results.py ]]; then
    python3 analyze_v6_results.py \
      --root "$ROOT" \
      --show-table \
      --make-plots || true
  fi
done

# Combined v13-specific post-analysis.
python3 - <<'PY'
from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt

outroot = Path("results/v13_CrossS_m10_FAC1_T1100_peakOvernight")
outdir = outroot / "combined_analysis"
outdir.mkdir(parents=True, exist_ok=True)

rows = []
for hist_path in sorted(outroot.glob("rate*/bs_on/T*_rho*/single_glider_history.csv")):
    ppath = hist_path.parent / "clean_arrhenius_params.json"
    if not ppath.exists():
        continue
    try:
        params = json.loads(ppath.read_text())
        hist = pd.read_csv(hist_path)
    except Exception as e:
        print(f"Skipping {hist_path}: {e}")
        continue
    if len(hist) < 20:
        continue

    tail = hist.tail(max(20, len(hist)//5))
    final = hist.iloc[-1]
    eps_total = float(final["eps_total"])
    eps_p = float(final["eps_plastic"])
    dt = float(params["dt"])

    rows.append({
        "run_dir": str(hist_path.parent),
        "rate_s": float(params["strain_rate"]),
        "dt_s": dt,
        "rho_m2": float(params["forest_rho_m2"]),
        "T_K": float(params["temperature_K"]),
        "CrossS_kB": float(params.get("expfit_cross_entropy_kB", float("nan"))),
        "PeierlsS_kB": float(params.get("expfit_peierls_entropy_kB", float("nan"))),
        "FAC": float(params.get("v11_cross_force_scale_factor", float("nan"))),
        "eps_total_final": eps_total,
        "eps_plastic_final": eps_p,
        "epsp_over_epstotal_final": eps_p / eps_total if eps_total else float("nan"),
        "tau_final_MPa": float(final["tau_MPa"]),
        "tau_tail_median_MPa": float(tail["tau_MPa"].median()),
        "tau_tail_abs_median_MPa": float(tail["tau_MPa"].abs().median()),
        "n_crossed_total_final": float(final["n_crossed_total"]),
        "sum_n_depin": float(hist["n_depin"].sum()),
        "n_live_pins_tail_median": float(tail["n_live_pins"].median()),
        "tau_local_median_tail_MPa": float(tail["tau_local_median_MPa"].median()),
        "tau_local_p90_tail_MPa": float(tail["tau_local_p90_MPa"].median()),
        "frac_tau_local_capped_tail": float(tail["frac_tau_local_capped"].median()),
        "max_crossing_rate_dt": float((hist["crossing_rate_max_s"] * dt).max()),
        "median_crossing_rate_dt_tail": float((tail["crossing_rate_max_s"] * dt).median()),
    })

df = pd.DataFrame(rows)
if df.empty:
    print("No histories found for combined analysis.")
    raise SystemExit(0)

df = df.sort_values(["rate_s", "rho_m2"])
csv = outdir / "v13_CrossS_m10_peakOvernight_summary.csv"
df.to_csv(csv, index=False)
print("\nCombined summary:")
cols = [
    "rate_s", "rho_m2", "tau_tail_median_MPa",
    "epsp_over_epstotal_final", "n_crossed_total_final",
    "n_live_pins_tail_median", "tau_local_median_tail_MPa",
    "max_crossing_rate_dt",
]
print(df[cols].to_string(index=False))
print(f"\nWrote {csv}")

plot_specs = [
    ("tau_tail_median_MPa", "Tail stress, MPa", "tau_vs_rho.png"),
    ("epsp_over_epstotal_final", "eps_p / eps_total", "epsp_ratio_vs_rho.png"),
    ("n_crossed_total_final", "Total depin/cross events", "crossings_vs_rho.png"),
    ("n_live_pins_tail_median", "Tail live pins", "live_pins_vs_rho.png"),
    ("tau_local_median_tail_MPa", "Tail median local stress, MPa", "tau_local_vs_rho.png"),
    ("max_crossing_rate_dt", "max(R_cross dt)", "event_resolution_vs_rho.png"),
]
for y, ylabel, fname in plot_specs:
    plt.figure()
    for rate, g in df.groupby("rate_s"):
        g = g.sort_values("rho_m2")
        plt.plot(g["rho_m2"], g[y], marker="o", label=f"{rate:g} s$^{{-1}}$")
    plt.xscale("log")
    plt.xlabel("forest density, m$^{-2}$")
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    p = outdir / fname
    plt.savefig(p, dpi=220)
    print(f"Wrote {p}")
PY

echo
echo "DONE. Combined analysis in:"
echo "  $OUTROOT/combined_analysis"
