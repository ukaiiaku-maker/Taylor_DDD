#!/usr/bin/env bash
# OVERNIGHT sweep for clean_arrhenius_taylor_ddd_v6.py
#
# KEY DESIGN DECISIONS (different from the v5 scaled-cell script):
#  * FIXED cell.  Obstacle count grows with density (nobs = rho * cell_area), which
#    is the only way the forest areal density -- and therefore the Taylor flow stress
#    -- can change with rho.  The scaled (1/sqrt(rho)) cell is self-similar and gives
#    density-INDEPENDENT results by construction; use it only for convergence checks.
#  * DECOUPLED barriers.  Peierls (lattice friction) sigc0 is LOW so the bulk flows at
#    low stress; crossing (forest) sigc0 is ~10x higher so pins hold and load up.
#  * Run to FLOW.  target strain is large enough that eps_p/eps_total -> ~0.9 and the
#    stress plateaus; only the plateau is a meaningful flow stress.
#  * strain rate 10 (not 1): with the low Peierls sigc0 the bulk yields easily, so a
#    higher rate is no longer "overdriven" and costs ~10x fewer steps.
#
# Resolution note: with NN=128, dz ~ 12 b, so obstacle spacing stays > dz up to
# ~1e16-1e17.  Keep rho <= 1e16 for physics; 1e17 is borderline, 1e18 is unresolved.
set -euo pipefail

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v6.py}"
ROOT="${ROOT:-$(pwd)/results/v6_overnight}"

# fixed cell + discretization
LX="${LX:-3060}"; LZ="${LZ:-1530}"; NN="${NN:-128}"; NLINE="${NLINE:-4}"; SEED="${SEED:-11}"
CAP="${CAP:-16}"            # capture radius ~ node spacing (fixed cell)
GLIDE_JUMP="${GLIDE_JUMP:-1.0}"
SOUT_M="${SOUT_M:-1e-7}"   # fixed out-of-plane spacing (clean for a fixed cell)

# loading
STRAIN_RATE="${STRAIN_RATE:-10}"
DT="${DT:-1e-8}"
TARGET_STRAIN="${TARGET_STRAIN:-0.012}"

# decoupled barriers (factor ~10 in stress scale)
PEIERLS_SIGC0="${PEIERLS_SIGC0:-150}"
CROSS_SIGC0="${CROSS_SIGC0:-1500}"
PEIERLS_SCALE="${PEIERLS_SCALE:-0.02}"
CROSS_SCALE="${CROSS_SCALE:-0.10}"
PEIERLS_ENTROPY_KB="${PEIERLS_ENTROPY_KB:--9.0}"
CROSS_ENTROPY_KB="${CROSS_ENTROPY_KB:--9.0}"
FLOOR_FRAC="${FLOOR_FRAC:-0.0}"

PIN_DIAG_EVERY="${PIN_DIAG_EVERY:-100}"

T_LIST="${T_LIST:-900 1100 1300}"
RHO_LIST="${RHO_LIST:-1e13 1e14 1e15 1e16}"
# backstress sweep: 'on' is the full grid; 'off' is the A/B control (kept smaller below)
BS_LIST="${BS_LIST:-on}"

run_one () {
  local T=$1 RHO=$2 BS=$3
  local OUT="$ROOT/bs_${BS}/T${T}_rho${RHO}"; mkdir -p "$OUT"
  echo "[$(date +%H:%M:%S)] bs=$BS T=$T rho=$RHO -> $OUT"
  python3 "$DRIVER" --outdir "$OUT" \
    --temperature-K "$T" --strain-rate "$STRAIN_RATE" \
    --target-strain "$TARGET_STRAIN" --dt "$DT" \
    --cell-lx-reduced "$LX" --cell-lz-reduced "$LZ" \
    --mobile-line-count "$NLINE" --mobile-line-nodes "$NN" \
    --forest-rho-m2 "$RHO" \
    --capture-radius-reduced "$CAP" --max-free-dx-reduced "$CAP" \
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
echo "=================================================================="
echo "v6 OVERNIGHT SWEEP  fixed cell ${LX}x${LZ}  NN=$NN"
echo "  peierls_sigc0=$PEIERLS_SIGC0  cross_sigc0=$CROSS_SIGC0 (ratio $(python3 -c "print($CROSS_SIGC0/$PEIERLS_SIGC0)"))"
echo "  rate=$STRAIN_RATE dt=$DT target=$TARGET_STRAIN  T=[$T_LIST]  rho=[$RHO_LIST]"
echo "=================================================================="

# 1) main grid: backstress ON, all T x rho
for T in $T_LIST; do for RHO in $RHO_LIST; do run_one "$T" "$RHO" on; done; done

# 2) A/B control: backstress OFF at one temperature, all rho
for RHO in $RHO_LIST; do run_one 1100 "$RHO" off; done

echo "=================================================================="
echo "ALL DONE. Summaries:"
echo "  python3 summarize_clean_arrhenius_taylor_v2.py --root $ROOT/bs_on  --tail-fraction 0.5"
echo "  python3 summarize_clean_arrhenius_taylor_v2.py --root $ROOT/bs_off --tail-fraction 0.5"
echo "  python3 diagnose_pin_amplification_v2.py       --root $ROOT/bs_on"
echo "=================================================================="
