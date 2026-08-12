#!/usr/bin/env python3
"""
Analytic Arrhenius-Taylor turnover triage for the v13 OpenDiS reduced DDD model.

Purpose
-------
Estimate where the peak / high-density turnover should occur before running
expensive DDD sweeps.  This is an analytical triage model, not a replacement for
DDD.  It uses the canonical Taylor Arrhenius prefactor and calibrates the
unknown effective stress-concentration factor from the v13 rho=3e16, T=1100 K,
strain_rate=10 s^-1 cross-force-factor diagnostic.

Assumed relation
----------------
  X = 1 / sqrt(2 rho)
  epsdot = nu0 * (b/X)^m * exp[-DeltaG(tau_loc,T)/(kT)]
  tau_loc = A(FAC) * tau_macro * X/b
  A(FAC) = A0 * FAC^(-p)

The EXP barrier is
  DeltaG = G0*exp[-a*(tau_loc/sigc)^n] - kT*S_kB
with floor_frac supported but defaulting to zero, matching the current runs.

Default calibration data are from the v13 actual-swept rho=3e16 runs:
  FAC=1.0:  tau_tail ~= 252.1117 MPa
  FAC=0.5:  tau_tail ~= 129.5741 MPa
  FAC=0.25: tau_tail ~= 85.0404 MPa
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

KB_EV = 8.617333262145e-5


def parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.replace(",", " ").split() if x.strip()]


def barrier_params(T: float, args):
    G0 = args.cross_scale * args.G00_eV * math.exp(-args.gT * (T - args.T0_K) / args.T0_K)
    Gfloor = args.floor_frac * G0
    sigc = args.sigc0_MPa * math.exp(-args.sT * (T - args.T0_K) / args.T0_K)
    entropy_term = -KB_EV * T * args.cross_entropy_kB
    return G0, Gfloor, sigc, entropy_term


def target_barrier_eV(rho: float, rate: float, T: float, args) -> float:
    X = 1.0 / math.sqrt(args.rho_prefactor * rho)
    bx = args.b_m / X
    pref = args.nu0_s * (bx ** args.site_prefactor_power)
    return KB_EV * T * math.log(pref / rate)


def tau_local_from_rate_MPa(rho: float, rate: float, T: float, args) -> float:
    G0, Gfloor, sigc, entropy_term = barrier_params(T, args)
    Gtar = target_barrier_eV(rho, rate, T, args)

    # Barrier range for nonnegative local stress.
    Gmax = entropy_term + G0
    Gmin = entropy_term + Gfloor

    if Gtar >= Gmax:
        return 0.0  # zero local stress already gives <= requested rate
    if Gtar <= Gmin:
        return float("nan")  # requested rate exceeds floor-limited maximum

    y = (Gtar - entropy_term - Gfloor) / max(G0 - Gfloor, 1e-300)
    y = min(max(y, 1e-300), 1.0)
    return sigc * ((-math.log(y) / args.a) ** (1.0 / args.n))


def calibrate_Afac(args):
    facs = np.asarray(parse_float_list(args.cal_facs), dtype=float)
    taus = np.asarray(parse_float_list(args.cal_taus_MPa), dtype=float)
    if len(facs) != len(taus):
        raise ValueError("--cal-facs and --cal-taus-MPa must have the same length")

    tau_loc = tau_local_from_rate_MPa(args.cal_rho_m2, args.cal_rate_s, args.temperature_K, args)
    X = 1.0 / math.sqrt(args.rho_prefactor * args.cal_rho_m2)
    bx = args.b_m / X
    A = tau_loc * bx / taus

    # ln A = ln A0 - p ln FAC
    coeff = np.polyfit(np.log(facs), np.log(A), deg=1)
    p = -float(coeff[0])
    A0 = float(math.exp(coeff[1]))
    pred = A0 * facs ** (-p)
    return A0, p, facs, taus, A, pred, tau_loc


def macro_tau_MPa(rho: float, rate: float, fac: float, T: float, A0: float, p: float, args) -> float:
    tau_loc = tau_local_from_rate_MPa(rho, rate, T, args)
    if not np.isfinite(tau_loc):
        return float("nan")
    X = 1.0 / math.sqrt(args.rho_prefactor * rho)
    bx = args.b_m / X
    Afac = A0 * fac ** (-p)
    return tau_loc * bx / Afac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="analytic_turnover_v13")
    ap.add_argument("--temperature-K", type=float, default=1100.0)
    ap.add_argument("--rates", default="0.01 0.1 1 10 100")
    ap.add_argument("--facs", default="1.0 0.5 0.25")
    ap.add_argument("--rho-min", type=float, default=1e13)
    ap.add_argument("--rho-max", type=float, default=1e19)
    ap.add_argument("--rho-n", type=int, default=900)

    # Current EXP barrier defaults.
    ap.add_argument("--b-m", type=float, default=2.48e-10)
    ap.add_argument("--nu0-s", type=float, default=1e12)
    ap.add_argument("--T0-K", type=float, default=1100.0)
    ap.add_argument("--G00-eV", type=float, default=1.908192)
    ap.add_argument("--gT", type=float, default=1.241743865563325)
    ap.add_argument("--sigc0-MPa", type=float, default=1497.042242375928)
    ap.add_argument("--sT", type=float, default=0.10850578873777168)
    ap.add_argument("--a", type=float, default=2.2056211004282904)
    ap.add_argument("--n", type=float, default=2.5207319790155385)
    ap.add_argument("--floor-frac", type=float, default=0.0)
    ap.add_argument("--cross-scale", type=float, default=0.40)
    ap.add_argument("--cross-entropy-kB", type=float, default=-9.0)

    # Taylor analytic geometry choices.
    ap.add_argument("--site-prefactor-power", type=float, default=4.0,
                    help="m in (b/X)^m; canonical Taylor value used previously was m=4")
    ap.add_argument("--rho-prefactor", type=float, default=2.0,
                    help="X=1/sqrt(rho_prefactor*rho); use 2 for X≈1/sqrt(2rho)")

    # v13 calibration point.
    ap.add_argument("--cal-rho-m2", type=float, default=3e16)
    ap.add_argument("--cal-rate-s", type=float, default=10.0)
    ap.add_argument("--cal-facs", default="1.0 0.5 0.25")
    ap.add_argument("--cal-taus-MPa", default="252.111671 129.574125 85.040429")

    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    A0, p, cal_facs, cal_taus, cal_A, cal_A_pred, cal_tau_loc = calibrate_Afac(args)

    rates = parse_float_list(args.rates)
    facs = parse_float_list(args.facs)
    rhos = np.logspace(math.log10(args.rho_min), math.log10(args.rho_max), args.rho_n)

    rows = []
    peak_rows = []
    for fac in facs:
        for rate in rates:
            taus = np.asarray([macro_tau_MPa(rho, rate, fac, args.temperature_K, A0, p, args) for rho in rhos])
            finite = np.isfinite(taus)
            for rho, tau in zip(rhos, taus):
                rows.append({"FAC": fac, "strain_rate_s": rate, "rho_m2": rho, "tau_MPa": tau})
            if np.any(finite):
                idxs = np.where(finite)[0]
                imax = idxs[np.nanargmax(taus[finite])]
                peak_rows.append({
                    "FAC": fac,
                    "strain_rate_s": rate,
                    "rho_peak_m2": rhos[imax],
                    "tau_peak_MPa": taus[imax],
                    "rho_first_softening_m2": rhos[min(imax + 1, len(rhos)-1)],
                    "tau_at_3e16_MPa": np.interp(3e16, rhos[finite], taus[finite]) if np.any(finite) else np.nan,
                })
            else:
                peak_rows.append({"FAC": fac, "strain_rate_s": rate, "rho_peak_m2": np.nan, "tau_peak_MPa": np.nan})

    df = pd.DataFrame(rows)
    pk = pd.DataFrame(peak_rows)
    df.to_csv(outdir / "analytic_stress_vs_rho.csv", index=False)
    pk.to_csv(outdir / "analytic_peak_summary.csv", index=False)

    with open(outdir / "calibration.txt", "w") as f:
        f.write(f"A0={A0:.8g}\n")
        f.write(f"p={p:.8g}\n")
        f.write(f"cal_tau_local_MPa={cal_tau_loc:.8g}\n")
        f.write("FAC, tau_measured_MPa, A_inferred, A_fit\n")
        for a_fac, tau, A, Ap in zip(cal_facs, cal_taus, cal_A, cal_A_pred):
            f.write(f"{a_fac:g}, {tau:.8g}, {A:.8g}, {Ap:.8g}\n")

    # Plot stress vs rho, one figure per FAC to keep readability.
    for fac in facs:
        plt.figure(figsize=(7.0, 5.0))
        for rate in rates:
            sub = df[(df["FAC"] == fac) & (df["strain_rate_s"] == rate)]
            plt.plot(sub["rho_m2"], sub["tau_MPa"], label=f"rate={rate:g} s$^{{-1}}$")
        plt.xscale("log")
        plt.xlabel(r"forest density $\rho$ (m$^{-2}$)")
        plt.ylabel(r"analytic $\tau$ (MPa)")
        plt.title(f"Analytic Arrhenius-Taylor turnover, FAC={fac:g}")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"analytic_stress_vs_rho_FAC_{str(fac).replace('.', 'p')}.png", dpi=200)
        plt.close()

    plt.figure(figsize=(7.0, 5.0))
    for fac in facs:
        sub = pk[pk["FAC"] == fac].sort_values("strain_rate_s")
        plt.plot(sub["strain_rate_s"], sub["rho_peak_m2"], marker="o", label=f"FAC={fac:g}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"strain rate $\dot\epsilon$ (s$^{-1}$)")
    plt.ylabel(r"predicted peak density $\rho_\mathrm{peak}$ (m$^{-2}$)")
    plt.title("Predicted density for peak / onset of turnover")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "analytic_peak_density_vs_rate.png", dpi=200)
    plt.close()

    print("Calibration:")
    print(f"  A(FAC)=A0*FAC^(-p), A0={A0:.4g}, p={p:.3g}")
    print("Peak summary:")
    print(pk.to_string(index=False))
    print(f"Wrote: {outdir}")


if __name__ == "__main__":
    main()
