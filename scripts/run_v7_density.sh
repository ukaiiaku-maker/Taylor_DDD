#!/usr/bin/env bash
# Shorter v7 sweep: ONE temperature, a series of densities, backstress on + off.
#
# v7 fixes the obstacle burn-through: a crossed (line,obstacle) pair is only BLOCKED
# until the lead node glides --rearm-radius-reduced away in x, then it can re-pin on the
# periodic wrap.  --min-pin-age-steps lets a fresh pin load (bow) before it can depin.
#
# IMPORTANT (from validation): the forest only gates the bulk once the pinned fraction is
# appreciable, which happens at HIGH density.  At 1e14 the forest is too sparse to move the
# flow stress; the Taylor trend should appear toward 1e16+.  So this sweep targets the high
# end and runs to flow.  The 1e16 case is the slow one (~2900 obstacles); 1e17 is slower
# still and mesh-marginal -- include it only if you have the time.
set -euo pipefail

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v7.py}"
ROOT="${ROOT:-$(pwd)/results/v7_density}"

LX="${LX:-3060}"; LZ="${LZ:-1530}"; NN="${NN:-96}"; NLINE="${NLINE:-4}"; SEED="${SEED:-11}"
CAP="${CAP:-16}"; REARM="${REARM:-40}"; MIN_PIN_AGE="${MIN_PIN_AGE:-5}"
GLIDE_JUMP="${GLIDE_JUMP:-1.0}"; SOUT_M="${SOUT_M:-1e-7}"

T="${T:-1100}"
STRAIN_RATE="${STRAIN_RATE:-10}"
DT="${DT:-1e-8}"
TARGET_STRAIN="${TARGET_STRAIN:-0.006}"

PEIERLS_SIGC0="${PEIERLS_SIGC0:-150}"
CROSS_SIGC0="${CROSS_SIGC0:-1500}"
PEIERLS_SCALE="${PEIERLS_SCALE:-0.02}"
CROSS_SCALE="${CROSS_SCALE:-0.40}"
PEIERLS_ENTROPY_KB="${PEIERLS_ENTROPY_KB:--9.0}"
CROSS_ENTROPY_KB="${CROSS_ENTROPY_KB:--9.0}"
FLOOR_FRAC="${FLOOR_FRAC:-0.0}"
PIN_DIAG_EVERY="${PIN_DIAG_EVERY:-50}"

RHO_LIST="${RHO_LIST:-1e14 1e15 1e16}"
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
    --capture-radius-reduced "$CAP" --max-free-dx-reduced "$CAP" \
    --rearm-radius-reduced "$REARM" --min-pin-age-steps "$MIN_PIN_AGE" \
    --glide-jump-length-reduced "$GLIDE_JUMP" \
    --out-of-plane-spacing-mode fixed --out-of-plane-spacing-m "$SOUT_M" \
    --expfit-peierls-sigc0-MPa "$PEIERLS_SIGC0" --expfit-cross-sigc0-MPa "$CROSS_SIGC0" \
    --expfit-peierls-scale "$PEIERLS_SCALE" --expfit-cross-scale "$CROSS_SCALE" \
    --expfit-peierls-entropy-kB "$PEIERLS_ENTROPY_KB" --expfit-cross-entropy-kB "$CROSS_ENTROPY_KB" \
    --expfit-floor-frac "$FLOOR_FRAC" \
    --pin-diagnostic-every "$PIN_DIAG_EVERY" \
    --backstress-mobility "$BS" --seed "$SEED" \
    > "$OUT/stdout.txt" 2> "$OUT/stderr.txt" \
    && echo "   ok" || { echo "   FAILED"; tail -5 "$OUT/stderr.txt"; }
}

mkdir -p "$ROOT"
echo "v7 density sweep  T=$T  rho=[$RHO_LIST]  bs=[$BS_LIST]  rearm=$REARM min_pin_age=$MIN_PIN_AGE cross_scale=$CROSS_SCALE"
for BS in $BS_LIST; do for RHO in $RHO_LIST; do run_one "$RHO" "$BS"; done; done

echo "Done. Analyze with:"
echo "  python3 analyze_v6_results.py --root $ROOT --show-table"
echo "Pass checks: n_live_pins_tail_median > 0 at high rho; tau rises with rho; bs_on > bs_off."
