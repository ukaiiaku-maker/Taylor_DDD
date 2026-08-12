#!/usr/bin/env python3
"""
analytic_cross_entropy_target.py

Small analytic selector for the v13 OpenDiS Arrhenius-Taylor sweep.

Purpose:
  Choose a modified forest-crossing entropy so that the Zener-Hollomon-like
  peak density and the pin-survival/loading timescale overlap near a target
  density.

This is a deliberately lightweight scaling calculation calibrated to the
current v13 observations:
  - T = 1100 K
  - current crossing entropy CrossS = -9 kB
  - current approximate analytical peak:
        rho_peak ~= 5.28e15 m^-2 at strain_rate = 0.1 s^-1
  - current pin-survival threshold at rho = 1e16:
        strain_rate_survival ~= 0.52 s^-1
"""

import math
import argparse
from pathlib import Path
import pandas as pd

KB_EV = 8.617333262145e-5

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=1100.0)
    ap.add_argument("--target-rho", type=float, default=5e15)
    ap.add_argument("--current-cross-entropy", type=float, default=-9.0)
    ap.add_argument("--candidate-cross-entropies", default="-9.0 -9.5 -10.0 -10.5 -11.0")
    ap.add_argument("--outdir", default="results/v13_crossEntropy_analytic")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    kT = KB_EV * args.T

    # Empirical current analytic scaling from previous scan:
    # rho_peak_current(rate) ~= 5.284764e15 * sqrt(rate/0.1)
    rho_ref = 5.284764e15
    rate_ref = 0.1

    # Current pin-survival threshold:
    # rate_survival_current(rho) ~= 0.52*sqrt(1e16/rho)
    surv_pref = 0.52

    rows = []
    for S in [float(x) for x in args.candidate_cross_entropies.split()]:
        # More negative entropy increases DeltaG by -kT*(S - S0).
        # The kinetic slowdown factor is exp(deltaG/kT) = exp(-(S-S0)).
        deltaG = -kT * (S - args.current_cross_entropy)
        slowdown = math.exp(deltaG / kT)

        # New ZH peak rate needed to put rho_peak at target rho.
        # New barrier means the same rho_peak occurs at current-rate * slowdown.
        rate_peak_target = rate_ref * (args.target_rho / rho_ref)**2 * slowdown

        # New survival threshold lowers by slowdown because zero-force lifetime increases.
        rate_survival_target = surv_pref * math.sqrt(1e16 / args.target_rho) / slowdown

        rows.append({
            "CrossS_kB": S,
            "deltaG0_eV_vs_current": deltaG,
            "slowdown_factor_R0": slowdown,
            "zero_force_lifetime_factor": slowdown,
            "rate_for_ZH_peak_at_target_rho_s^-1": rate_peak_target,
            "rate_for_pin_survival_at_target_rho_s^-1": rate_survival_target,
            "ratio_peak_rate_over_survival_rate": rate_peak_target / rate_survival_target,
        })

    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    out = outdir / "cross_entropy_target_table.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    print("\nRecommended starting point:")
    print("  CrossS = -10 kB")
    print("  It adds about +0.095 eV to the zero-force crossing barrier at 1100 K,")
    print("  moving the predicted peak near rho ~5e15 to rate ~0.25-0.30 s^-1,")
    print("  while also moving the pin-survival threshold into the same rate window.")
    print("  Suggested rates: 0.15, 0.30, 0.60 s^-1.")

if __name__ == "__main__":
    main()
