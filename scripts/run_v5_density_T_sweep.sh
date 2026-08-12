#!/usr/bin/env bash
# Density x temperature sweep for clean_arrhenius_taylor_ddd_v5.py
#
# Cell scales as 1/sqrt(rho) so N_obs ~ MX*MZ and nodes/spacing stay constant.
#
# FLOW-STRESS TUNING (read this):
#   The flow stress is set by the barrier's CHARACTERISTIC STRESS, --expfit-sigc0-MPa,
#   NOT by the barrier height (--expfit-*-scale).  At this geometry/rate the Peierls
#   friction is ~0.11*sigc0, so sigc0=120 -> ~13 MPa friction floor; the forest/Taylor
#   contribution sits on top of that.  Lower SIGC0 to lower the whole stress scale.
#   sigc0 is SHARED by the Peierls (glide) and crossing (depin) branches.
#
# STRAIN RATE / TIME STEP:
#   STRAIN_RATE=1 keeps the system from being "overdriven" (at high rate the Peierls
#   branch cannot supply the required glide rate and the cell never yields).  Because
#   dt must resolve the kinetics, low rate means many steps; DT=1e-8 is a tractable
#   compromise (coarsens depin-event timing ~10x vs 1e-9).  For a final quantitative
#   number, re-run a few cases at DT=1e-9 to check dt-convergence.
set -euo pipefail

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v5.py}"
ROOT="${ROOT:-$(pwd)/results/v5_density_T_sweep}"
B="${B:-2.48e-10}"
MX="${MX:-24}"; MZ="${MZ:-12}"; NN="${NN:-128}"; NLINE="${NLINE:-4}"; SEED="${SEED:-11}"

STRAIN_RATE="${STRAIN_RATE:-1}"
DT="${DT:-1e-8}"
TARGET_STRAIN="${TARGET_STRAIN:-0.0002}"
BACKSTRESS="${BACKSTRESS:-on}"
GLIDE_JUMP="${GLIDE_JUMP:-1.0}"

# barrier knobs (tunable)
SIGC0="${SIGC0:-120}"
CROSS_SCALE="${CROSS_SCALE:-0.10}"
PEIERLS_SCALE="${PEIERLS_SCALE:-0.02}"
CROSS_ENTROPY_KB="${CROSS_ENTROPY_KB:--9.0}"
PEIERLS_ENTROPY_KB="${PEIERLS_ENTROPY_KB:--9.0}"
FLOOR_FRAC="${FLOOR_FRAC:-0.0}"

# out-of-plane spacing (plastic-strain normalization): fixed is cleaner for a
# density sweep; set SOUT_MODE=forest_spacing to scale it with 1/sqrt(rho).
SOUT_MODE="${SOUT_MODE:-fixed}"
SOUT_M="${SOUT_M:-1e-7}"

T_LIST="${T_LIST:-900 1100 1300}"
RHO_LIST="${RHO_LIST:-1e11 1e12 1e13 1e14 1e15 1e16 1e17 1e18}"

mkdir -p "$ROOT"
echo "sigc0=$SIGC0  rate=$STRAIN_RATE  dt=$DT  target=$TARGET_STRAIN  backstress=$BACKSTRESS  s_out=$SOUT_MODE/$SOUT_M"
for T in $T_LIST; do
  for RHO in $RHO_LIST; do
    read LX LZ CAP < <(python3 - "$RHO" "$B" "$MX" "$MZ" "$NN" <<'PY'
import sys, math
rho, b, mx, mz, nn = (float(a) for a in sys.argv[1:6])
sp = (1.0/math.sqrt(rho))/b
print(f"{mx*sp:.6g} {mz*sp:.6g} {2.0*(mz*sp/nn):.6g}")
PY
)
    OUT="$ROOT/T${T}_rho${RHO}"; mkdir -p "$OUT"
    echo "T=$T rho=$RHO  Lx=$LX Lz=$LZ cap=$CAP"
    python3 "$DRIVER" --outdir "$OUT" \
      --temperature-K "$T" --strain-rate "$STRAIN_RATE" \
      --target-strain "$TARGET_STRAIN" --dt "$DT" \
      --cell-lx-reduced "$LX" --cell-lz-reduced "$LZ" \
      --mobile-line-count "$NLINE" --mobile-line-nodes "$NN" \
      --forest-rho-m2 "$RHO" \
      --capture-radius-reduced "$CAP" --max-free-dx-reduced "$CAP" \
      --glide-jump-length-reduced "$GLIDE_JUMP" \
      --out-of-plane-spacing-mode "$SOUT_MODE" --out-of-plane-spacing-m "$SOUT_M" \
      --expfit-sigc0-MPa "$SIGC0" \
      --expfit-cross-scale "$CROSS_SCALE" --expfit-peierls-scale "$PEIERLS_SCALE" \
      --expfit-cross-entropy-kB "$CROSS_ENTROPY_KB" --expfit-peierls-entropy-kB "$PEIERLS_ENTROPY_KB" \
      --expfit-floor-frac "$FLOOR_FRAC" \
      --backstress-mobility "$BACKSTRESS" --seed "$SEED" \
      > "$OUT/stdout.txt" 2> "$OUT/stderr.txt" \
      && echo "   ok" || { echo "   FAILED"; tail -5 "$OUT/stderr.txt"; }
  done
done
echo "sweep done -> $ROOT"
