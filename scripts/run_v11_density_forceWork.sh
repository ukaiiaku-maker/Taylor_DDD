#!/usr/bin/env bash
# v11 density sweep for the Arrhenius Taylor/Peierls DDD harness.
# Main change relative to run_v7_density.sh:
#   - lower the effective Peierls/free-glide baseline by not assigning the same
#     large negative entropy to Peierls glide that is assigned to forest crossing;
#   - keep forest crossing slow enough to test whether pins become rate-limiting;
#   - use v8 operator splitting so a newly captured pin cannot depin in the same step.
set -euo pipefail

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v11.py}"
ROOT="${ROOT:-$(pwd)/results/v11_density_forceWork}"

LX="${LX:-3060}"; LZ="${LZ:-1530}"; NN="${NN:-96}"; NLINE="${NLINE:-4}"; SEED="${SEED:-11}"
CAP="${CAP:-16}"; REARM="${REARM:-40}"; MIN_PIN_AGE="${MIN_PIN_AGE:-0}"; MAX_FREE_DX="${MAX_FREE_DX:-0.5}"
GLIDE_JUMP="${GLIDE_JUMP:-1.0}"; SOUT_M="${SOUT_M:-1e-7}"

T="${T:-1100}"
STRAIN_RATE="${STRAIN_RATE:-10}"
DT="${DT:-1e-8}"
TARGET_STRAIN="${TARGET_STRAIN:-0.006}"

# Keep the 10x stress-scale separation, but make the Peierls branch a low-friction
# baseline.  With these defaults at T=1100 K, rate=10 s^-1, the Peierls-only
# diagnostic stress is ~1.2 MPa instead of ~43 MPa.
PEIERLS_SIGC0="${PEIERLS_SIGC0:-150}"
CROSS_SIGC0="${CROSS_SIGC0:-1500}"
PEIERLS_SCALE="${PEIERLS_SCALE:-0.02}"
CROSS_SCALE="${CROSS_SCALE:-0.40}"
PEIERLS_ENTROPY_KB="${PEIERLS_ENTROPY_KB:-0.0}"
CROSS_ENTROPY_KB="${CROSS_ENTROPY_KB:--9.0}"
FLOOR_FRAC="${FLOOR_FRAC:-0.0}"
PIN_DIAG_EVERY="${PIN_DIAG_EVERY:-100}"
BACKSTRESS_COM_PROJECTION="${BACKSTRESS_COM_PROJECTION:-external_drive}"
CROSSING_DRIVE_MODE="${CROSSING_DRIVE_MODE:-force_work}"
CROSS_FORCE_SCALE_MODE="${CROSS_FORCE_SCALE_MODE:-line_tension}"
CROSS_FORCE_SCALE_FACTOR="${CROSS_FORCE_SCALE_FACTOR:-1.0}"
CROSS_FORCE_SCALE_N="${CROSS_FORCE_SCALE_N:--1.0}"

# Start with the current v7 densities plus one intermediate high-density point.
# 1e17 is expensive and mesh-marginal; append it manually if the 3e16 run is useful.
RHO_LIST="${RHO_LIST:-1e14 1e15 3e15 1e16 3e16}"
BS_LIST="${BS_LIST:-on off}"

run_one () {
  local RHO=$1 BS=$2
  local OUT="$ROOT/bs_${BS}/T${T}_rho${RHO}"; mkdir -p "$OUT"
  echo "[$(date +%H:%M:%S)] bs=$BS T=$T rho=$RHO -> $OUT"
  python3 "$DRIVER" --outdir "$OUT" \
    --temperature-K "$T" --strain-rate "$STRAIN_RATE" \
    --target-strain "$TARGET_STRAIN" --dt "$DT" \
    --cell-lx-reduced "$LX" --cell-lz-reduced "$LZ" \
    --mobile-line-count "$NLINE" --mobile-line-nodes "$NN" \
    --forest-rho-m2 "$RHO" \
    --capture-radius-reduced "$CAP" --capture-mode swept_crossing --snap-swept-capture-to-obstacle \
    --max-free-dx-reduced "$MAX_FREE_DX" \
    --rearm-radius-reduced "$REARM" --min-pin-age-steps "$MIN_PIN_AGE" \
    --glide-jump-length-reduced "$GLIDE_JUMP" \
    --out-of-plane-spacing-mode fixed --out-of-plane-spacing-m "$SOUT_M" \
    --expfit-peierls-sigc0-MPa "$PEIERLS_SIGC0" --expfit-cross-sigc0-MPa "$CROSS_SIGC0" \
    --expfit-peierls-scale "$PEIERLS_SCALE" --expfit-cross-scale "$CROSS_SCALE" \
    --expfit-peierls-entropy-kB "$PEIERLS_ENTROPY_KB" --expfit-cross-entropy-kB "$CROSS_ENTROPY_KB" \
    --expfit-floor-frac "$FLOOR_FRAC" \
    --pin-diagnostic-every "$PIN_DIAG_EVERY" \
    --crossing-drive-mode "$CROSSING_DRIVE_MODE" \
    --cross-force-scale-mode "$CROSS_FORCE_SCALE_MODE" \
    --cross-force-scale-factor "$CROSS_FORCE_SCALE_FACTOR" \
    --cross-force-scale-N "$CROSS_FORCE_SCALE_N" \
    --backstress-mobility "$BS" --backstress-com-projection "$BACKSTRESS_COM_PROJECTION" --seed "$SEED" \
    > "$OUT/stdout.txt" 2> "$OUT/stderr.txt" \
    && echo "   ok" || { echo "   FAILED"; tail -20 "$OUT/stderr.txt"; return 1; }
}

mkdir -p "$ROOT"
echo "v11 density sweep T=$T rho=[$RHO_LIST] bs=[$BS_LIST] PeierlsS=$PEIERLS_ENTROPY_KB CrossS=$CROSS_ENTROPY_KB MAX_FREE_DX=$MAX_FREE_DX capture=swept_crossing backstress_projection=$BACKSTRESS_COM_PROJECTION crossing_drive=$CROSSING_DRIVE_MODE force_scale=$CROSS_FORCE_SCALE_MODE factor=$CROSS_FORCE_SCALE_FACTOR"
for BS in $BS_LIST; do for RHO in $RHO_LIST; do run_one "$RHO" "$BS"; done; done

echo "Done. Analyze with:"
echo "  python3 analyze_v6_results.py --root $ROOT --show-table"
echo "Also inspect clean_arrhenius_params.json for v8_peierls_only_baseline_tau_MPa, v10_backstress_com_projection, and v11_cross_force_scale_source."
