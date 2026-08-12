#!/usr/bin/env bash
# v17 systematic unattended rate-temperature-density campaign.
#
# Purpose:
#   Systematic physical controls: strain rate x temperature x density
#   Fixed barrier/kinetic model: CrossS=-9.25, floor=0.50, FAC=0.25
#   v17 uncapped physical L_eff: tau_local_cap_mode=none, tau_local_length_mode=feed
#
# Expected use:
#   cd /Users/sdillon/OpenDis
#   source .venv-opendis/bin/activate
#   chmod +x run_v17_systematic_rate_temp_campaign.sh
#   nohup ./run_v17_systematic_rate_temp_campaign.sh > run_v17_systematic_rate_temp_campaign.nohup.log 2>&1 &
#
# This script is resumable: if a run_summary.txt already exists for a run, it is skipped.

set -u
set -o pipefail

cd "$(dirname "$0")"

DRIVER=clean_arrhenius_taylor_ddd_v17.py
ANALYZE=analyze_v6_results.py
BURST=analyze_depin_burst_statistics.py

ROOT=results/v17_systematic_T_rate_density_FAC0p25_CrossS_m9p25_floor0p50_LeffFeed
STATUS="$ROOT/campaign_status.log"

TARGET_STRAIN=0.006
MAX_FREE_DX=0.5

FAC=0.25
CROSS_S=-9.25
PEIERLS_S=0.0
CROSS_FLOOR=0.50
PEIERLS_FLOOR=0.0

T_LIST="1000 1100 1200"
RATE_LIST="0.15 0.45 1.5 4.5"
RHO_LIST="1e12 3e12 1e13 3e13 1e14 3e14 1e15 3e15 1e16 3e16"

mkdir -p "$ROOT"

log() {
  echo "$@" | tee -a "$STATUS"
}

rate_tag() {
  echo "$1" | sed 's/\./p/g; s/+//g; s/-/m/g'
}

dt_for() {
  local RATE="$1"
  local T="$2"
  local DT
  case "$RATE" in
    0.15) DT="1e-7" ;;
    0.45) DT="1e-7" ;;
    1.5)  DT="3e-8" ;;
    4.5)  DT="1e-8" ;;
    *) echo "Unknown RATE=$RATE" >&2; return 2 ;;
  esac

  if [[ "$T" == "1200" ]]; then
    if [[ "$RATE" == "0.15" || "$RATE" == "0.45" ]]; then
      DT="5e-8"
    fi
  fi
  echo "$DT"
}

run_one() {
  local T="$1"
  local RATE="$2"
  local RHO="$3"
  local DT="$4"
  local RTAG
  RTAG="$(rate_tag "$RATE")"

  local BLOCK="$ROOT/T${T}_rate${RTAG}"
  local OUTDIR="$BLOCK/bs_on/T${T}_rho${RHO}"
  mkdir -p "$OUTDIR"

  log ""
  log "[$(date '+%Y-%m-%d %H:%M:%S')] START T=$T rate=$RATE dt=$DT rho=$RHO"
  log "OUTDIR=$OUTDIR"

  if [[ -s "$OUTDIR/run_summary.txt" ]] && grep -q "Run summary" "$OUTDIR/run_summary.txt"; then
    log "SKIP existing completed run: $OUTDIR"
    return 0
  fi

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

  local RET=$?
  if [[ "$RET" -ne 0 ]]; then
    log "FAILED T=$T rate=$RATE rho=$RHO return_code=$RET"
    echo "FAILED T=$T rate=$RATE rho=$RHO return_code=$RET" >> "$ROOT/failed_runs.txt"
    return "$RET"
  fi

  log "[$(date '+%Y-%m-%d %H:%M:%S')] DONE T=$T rate=$RATE rho=$RHO"
  if [[ -s "$OUTDIR/run_summary.txt" ]]; then
    log "--- run_summary.txt ---"
    cat "$OUTDIR/run_summary.txt" | tee -a "$STATUS"
  fi

  if [[ -f "$BURST" ]]; then
    log "--- burst diagnostic for T=$T rate=$RATE rho=$RHO ---"
    python3 "$BURST" \
      --run-dir "$OUTDIR" \
      --cluster-gap-steps 1 \
      --active-plastic-ratio 1.0 \
      --stress-drop-threshold-MPa 0.0 \
      --n-boot 100 \
      --show-table | tee -a "$STATUS" || log "WARNING: burst diagnostic failed for $OUTDIR"
  fi

  return 0
}

analyze_block() {
  local T="$1"
  local RATE="$2"
  local RTAG
  RTAG="$(rate_tag "$RATE")"
  local BLOCK="$ROOT/T${T}_rate${RTAG}"

  log ""
  log "[$(date '+%Y-%m-%d %H:%M:%S')] ANALYZE BLOCK T=$T rate=$RATE"

  if [[ -f "$ANALYZE" ]]; then
    python3 "$ANALYZE" --root "$BLOCK" --show-table | tee -a "$STATUS" || log "WARNING: analyze_v6 failed for $BLOCK"
  fi

  if [[ -f "$BURST" ]]; then
    python3 "$BURST" \
      --root "$BLOCK" \
      --cluster-gap-steps 1 \
      --active-plastic-ratio 1.0 \
      --stress-drop-threshold-MPa 0.0 \
      --n-boot 200 \
      --show-table | tee -a "$STATUS" || log "WARNING: burst summary failed for $BLOCK"
  fi
}

aggregate_final() {
  log ""
  log "[$(date '+%Y-%m-%d %H:%M:%S')] FINAL AGGREGATION"

  python3 - <<'PY' | tee -a "$ROOT/campaign_status.log"
from pathlib import Path
import re
import pandas as pd
import numpy as np

ROOT = Path("results/v17_systematic_T_rate_density_FAC0p25_CrossS_m9p25_floor0p50_LeffFeed")

def parse_block(p: Path):
    m = re.search(r"T([0-9.]+)_rate([0-9pm]+)", p.name)
    if not m:
        return None, None
    T = float(m.group(1))
    rate = float(m.group(2).replace("p", ".").replace("m", "-"))
    return T, rate

summary_rows = []
burst_rows = []

for block in sorted(ROOT.glob("T*_rate*")):
    T, rate = parse_block(block)
    if T is None:
        continue

    s = block / "analysis" / "summary_v6.csv"
    if s.exists():
        try:
            df = pd.read_csv(s)
            df["T_block_K"] = T
            df["rate_block_s"] = rate
            df["block"] = block.name
            summary_rows.append(df)
        except Exception as e:
            print(f"Could not read {s}: {e}")

    b = block / "depin_burst_summary_all.csv"
    if b.exists():
        try:
            df = pd.read_csv(b)
            df["T_block_K"] = T
            df["rate_block_s"] = rate
            df["block"] = block.name
            burst_rows.append(df)
        except Exception as e:
            print(f"Could not read {b}: {e}")

if summary_rows:
    allsum = pd.concat(summary_rows, ignore_index=True)
    allsum.to_csv(ROOT / "campaign_summary_v6_all.csv", index=False)
    print(f"Wrote {ROOT / 'campaign_summary_v6_all.csv'}")
    keep = [c for c in [
        "T_block_K", "rate_block_s", "rho_m2", "tau_tail_abs_median_MPa",
        "epsp_over_epstotal_final", "n_crossed_total_final",
        "n_live_pins_tail_median", "tau_local_median_tail_MPa",
        "frac_tau_local_capped_tail", "flags"
    ] if c in allsum.columns]
    print()
    print("Campaign mechanical summary:")
    print(allsum[keep].sort_values(["T_block_K","rate_block_s","rho_m2"]).to_string(index=False))
else:
    print("No summary_v6.csv files found.")

if burst_rows:
    allburst = pd.concat(burst_rows, ignore_index=True)
    allburst.to_csv(ROOT / "campaign_depin_burst_summary_all.csv", index=False)
    print()
    print(f"Wrote {ROOT / 'campaign_depin_burst_summary_all.csv'}")
    keep = [c for c in [
        "T_block_K", "rate_block_s", "rho_m2", "total_depin", "n_bursts",
        "largest_burst_depin", "largest_burst_fraction", "p90_burst_depin",
        "p99_burst_depin", "max_stress_drop_negsum_MPa",
        "frac_depin_events_at_tau_cap", "fit_status", "power_alpha",
        "loglik_power_minus_exp", "null_largest_fraction_p_ge_obs",
        "null_fano_p_ge_obs", "burst_distribution_suggestive", "interpretation"
    ] if c in allburst.columns]
    print()
    print("Campaign burst summary:")
    print(allburst[keep].sort_values(["T_block_K","rate_block_s","rho_m2"]).to_string(index=False))

    x = allburst.copy()
    for col in ["largest_burst_fraction", "p99_burst_depin", "max_stress_drop_negsum_MPa", "total_depin"]:
        if col not in x.columns:
            x[col] = np.nan
    x["burst_score"] = (
        np.log10(x["total_depin"].fillna(0) + 1.0)
        * x["largest_burst_fraction"].fillna(0)
        * np.log10(x["p99_burst_depin"].fillna(0) + 1.0)
    )
    x = x.sort_values("burst_score", ascending=False)
    x.to_csv(ROOT / "campaign_burst_candidates_ranked.csv", index=False)
    print()
    print(f"Wrote {ROOT / 'campaign_burst_candidates_ranked.csv'}")
    keep = [c for c in [
        "T_block_K", "rate_block_s", "rho_m2", "burst_score",
        "largest_burst_fraction", "p99_burst_depin",
        "max_stress_drop_negsum_MPa", "total_depin",
        "fit_status", "loglik_power_minus_exp", "interpretation"
    ] if c in x.columns]
    print()
    print("Top burst candidates:")
    print(x[keep].head(20).to_string(index=False))
else:
    print("No depin_burst_summary_all.csv files found.")
PY
}

log "============================================================"
log "v17 systematic unattended rate-temperature-density campaign"
log "ROOT=$ROOT"
log "DRIVER=$DRIVER"
log "TARGET_STRAIN=$TARGET_STRAIN"
log "FAC=$FAC"
log "CrossS=$CROSS_S"
log "cross_floor_frac=$CROSS_FLOOR"
log "tau_local_cap_mode=none"
log "tau_local_length_mode=feed"
log "T_LIST=$T_LIST"
log "RATE_LIST=$RATE_LIST"
log "RHO_LIST=$RHO_LIST"
log "Started: $(date '+%Y-%m-%d %H:%M:%S')"
log "============================================================"

python3 -m py_compile "$DRIVER" || { log "FATAL: driver does not compile"; exit 1; }
python3 -m py_compile "$ANALYZE" || log "WARNING: $ANALYZE does not compile"
python3 -m py_compile "$BURST" || log "WARNING: $BURST does not compile"

for T in $T_LIST; do
  for RATE in $RATE_LIST; do
    DT="$(dt_for "$RATE" "$T")"
    RTAG="$(rate_tag "$RATE")"
    log ""
    log "############################################################"
    log "BLOCK START T=$T rate=$RATE dt=$DT block=T${T}_rate${RTAG}"
    log "############################################################"

    for RHO in $RHO_LIST; do
      run_one "$T" "$RATE" "$RHO" "$DT" || true
    done

    analyze_block "$T" "$RATE"
    aggregate_final

    log "############################################################"
    log "BLOCK DONE T=$T rate=$RATE"
    log "############################################################"
  done
done

aggregate_final

date '+%Y-%m-%d %H:%M:%S' > "$ROOT/campaign_finished.ok"
log ""
log "DONE: $ROOT"
log "Finished: $(cat "$ROOT/campaign_finished.ok")"
