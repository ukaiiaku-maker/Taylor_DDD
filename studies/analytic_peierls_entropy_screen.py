#!/usr/bin/env python3
import numpy as np
import pandas as pd

kB_eV = 8.617333262145e-5

T0 = 1100.0
G00_eV = 1.908192
gT = 1.241743865563325
sigc0_MPa = 1497.042242375928
sT = 0.10850578873777168
a = 2.2056211004282904
n = 2.5207319790155385

def G0(T, scale):
    return scale * G00_eV * np.exp(-gT*(T-T0)/T0)

def sigc(T):
    return sigc0_MPa * np.exp(-sT*(T-T0)/T0)

def H_expfit(tau_MPa, T, scale):
    r = max(tau_MPa, 0.0) / max(sigc(T), 1e-300)
    return G0(T, scale) * np.exp(-a*r**n)

def G_total(tau_MPa, T, scale, S_kB):
    return max(0.0, H_expfit(tau_MPa, T, scale) - kB_eV*T*S_kB)

def main():
    rows = []
    for scale in [0.002,0.005,0.01,0.015,0.02,0.03,0.04,0.05]:
        for S in np.linspace(-12, -5, 29):
            G1000_0 = G_total(0, 1000, scale, S)
            G1100_0 = G_total(0, 1100, scale, S)
            G1200_0 = G_total(0, 1200, scale, S)
            rows.append({
                "scale": scale,
                "S_kB": S,
                "G0_1000_eV": G1000_0,
                "G0_1100_eV": G1100_0,
                "G0_1200_eV": G1200_0,
                "G0_1000_over_kT": G1000_0/(kB_eV*1000),
                "G0_1100_over_kT": G1100_0/(kB_eV*1100),
                "G0_1200_over_kT": G1200_0/(kB_eV*1200),
            })

    df = pd.DataFrame(rows)

    # Useful window: not zero-stress saturated at 1100/1200,
    # but not huge at 1000.
    cand = df[
        (df["G0_1000_over_kT"] > 5) &
        (df["G0_1000_over_kT"] < 15) &
        (df["G0_1100_over_kT"] > 5) &
        (df["G0_1200_over_kT"] > 5)
    ].copy()

    print("\nCandidates with finite zero-stress barriers through 1200 K:")
    print(cand.sort_values(["scale","S_kB"]).to_string(index=False))

    df.to_csv("results/analytic_peierls_entropy_screen.csv", index=False)
    cand.to_csv("results/analytic_peierls_entropy_candidates.csv", index=False)
    print("\nwrote results/analytic_peierls_entropy_screen.csv")
    print("wrote results/analytic_peierls_entropy_candidates.csv")

if __name__ == "__main__":
    main()
