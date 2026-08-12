#!/usr/bin/env bash
# Density x temperature sweep for clean_arrhenius_taylor_ddd_v4.py
#
# Geometry consistency (NOT automatic in the driver):
#   The driver fixes the cell and scales only the obstacle COUNT with density.
#   Here we instead scale the CELL as 1/sqrt(rho) so that at EVERY density the
#   cell spans the same number of obstacle spacings.  Consequences:
#     - obstacle count is constant  (N_obs ~ MX*MZ ~ 288)
#     - nodes-per-spacing is constant (~ NN/MZ ~ 21)
#     - capture radius and max glide step scale with the spacing, so obstacles
#       are neither missed (capture >= node spacing) nor skipped (step <= capture)
#   Density then enters through spacing/geometry, exactly as intended.
set -euo pipefail

DRIVER="${DRIVER:-clean_arrhenius_taylor_ddd_v4.py}"
ROOT="${ROOT:-$(pwd)/results/v4_density_T_sweep}"
B="${B:-2.48e-10}"          # Burgers vector (m)
MX="${MX:-24}"              # cell width  in obstacle spacings
MZ="${MZ:-12}"              # cell height in obstacle spacings  (MX*MZ ~ N_obstacles)
NN="${NN:-256}"             # nodes per line (~ NN/MZ nodes per spacing)
NLINE="${NLINE:-4}"
SEED="${SEED:-11}"
TARGET_STRAIN="${TARGET_STRAIN:-0.00025}"
DT="${DT:-1e-9}"
STRAIN_RATE="${STRAIN_RATE:-100}"
BACKSTRESS="${BACKSTRESS:-on}"   # 'on' = v4 force balance; 'off' = v3 baseline for A/B

T_LIST="${T_LIST:-900 1100 1300}"
RHO_LIST="${RHO_LIST:-1e12 1e13 1e14 1e15 1e16 1e17 1e18}"

mkdir -p "$ROOT"
for T in $T_LIST; do
  for RHO in $RHO_LIST; do
    read LX LZ CAP MAXDX < <(python3 - "$RHO" "$B" "$MX" "$MZ" "$NN" <<'PY'
import sys, math
rho, b, mx, mz, nn = (float(a) for a in sys.argv[1:6])
spacing_red = (1.0/math.sqrt(rho))/b      # obstacle spacing in units of b
lx = mx*spacing_red
lz = mz*spacing_red
cap = 2.0*(lz/nn)                          # capture radius ~ 2 node spacings
maxdx = cap                                # CFL: never step past an obstacle
print(f"{lx:.6g} {lz:.6g} {cap:.6g} {maxdx:.6g}")
PY
)
    OUT="$ROOT/T${T}_rho${RHO}"
    mkdir -p "$OUT"
    echo "T=$T rho=$RHO  Lx=$LX Lz=$LZ cap=$CAP maxdx=$MAXDX"
    python3 "$DRIVER" --outdir "$OUT" \
      --temperature-K "$T" --strain-rate "$STRAIN_RATE" \
      --target-strain "$TARGET_STRAIN" --dt "$DT" \
      --cell-lx-reduced "$LX" --cell-lz-reduced "$LZ" \
      --mobile-line-count "$NLINE" --mobile-line-nodes "$NN" \
      --forest-rho-m2 "$RHO" \
      --capture-radius-reduced "$CAP" --max-free-dx-reduced "$MAXDX" \
      --out-of-plane-spacing-mode forest_spacing \
      --backstress-mobility "$BACKSTRESS" \
      --pin-diagnostic-every 0 \
      --seed "$SEED" \
      > "$OUT/stdout.txt" 2> "$OUT/stderr.txt" \
      && echo "   ok" || { echo "   FAILED"; tail -5 "$OUT/stderr.txt"; }
  done
done
echo "sweep done -> $ROOT"
