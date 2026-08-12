#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

kB_eV = 8.617333262145e-5

# EXP posterior median parameters.
T0 = 1100.0
G00_eV_raw = 1.908192
gT = 1.241743865563325
sigc0_MPa = 1497.042242375928
sT = 0.10850578873777168
a = 2.2056211004282904
n = 2.5207319790155385
b_m = 2.48e-10

def G0_eV(T, scale=1.0):
    return scale * G00_eV_raw * np.exp(-gT*(T-T0)/T0)

def sigc_MPa(T):
    return sigc0_MPa * np.exp(-sT*(T-T0)/T0)

def G_exp_eV(tau_local_MPa, T, scale=1.0, floor_frac=0.0):
    tau = np.maximum(np.asarray(tau_local_MPa, dtype=float), 0.0)
    G0 = G0_eV(T, scale)
    Gf = floor_frac * G0
    r = tau / max(sigc_MPa(T), 1e-300)
    return Gf + (G0 - Gf) * np.exp(-a*r**n)

def vstar_eV_per_MPa(tau_local_MPa, T, scale=1.0, floor_frac=0.0):
    # vstar = -dG/dtau, units eV/MPa.
    tau = np.maximum(np.asarray(tau_local_MPa, dtype=float), 1e-300)
    G0 = G0_eV(T, scale)
    Gf = floor_frac * G0
    sc = max(sigc_MPa(T), 1e-300)
    r = tau/sc
    return (G0 - Gf) * np.exp(-a*r**n) * a*n*r**(n-1.0) / sc

def invert_G_to_tau_local(Greq_eV, T, scale=1.0, floor_frac=0.0):
    G0 = G0_eV(T, scale)
    Gf = floor_frac * G0
    if Greq_eV > G0:
        return np.nan, "Greq_above_G0"
    if Greq_eV <= Gf:
        return np.inf, "Greq_below_floor"
    y = (Greq_eV - Gf) / max(G0 - Gf, 1e-300)
    r = (-np.log(np.clip(y, 1e-300, 1.0)) / a)**(1.0/n)
    return sigc_MPa(T) * r, "ok"

def analytic_curve(T, edot, rho_grid, scale, floor_frac, p, nu0):
    rows = []
    for rho in rho_grid:
        X_over_b = 1.0/(b_m*np.sqrt(rho))
        pref = nu0 * (1.0/X_over_b)**p   # (b/X)^p
        Greq = kB_eV*T*np.log(max(pref/edot, 1e-300))

        tau_local, status = invert_G_to_tau_local(Greq, T, scale, floor_frac)
        if np.isfinite(tau_local):
            tau_app = tau_local / X_over_b
            Gval = G_exp_eV(tau_local, T, scale, floor_frac)
            vloc = vstar_eV_per_MPa(tau_local, T, scale, floor_frac)
            peak_metric = tau_local * vloc / (kB_eV*T)
        else:
            tau_app = np.nan
            Gval = np.nan
            vloc = np.nan
            peak_metric = np.nan

        rows.append({
            "T_K": T,
            "edot_s": edot,
            "rho_m2": rho,
            "barrier_scale": scale,
            "floor_frac": floor_frac,
            "taylor_power": p,
            "X_over_b": X_over_b,
            "prefactor_s": pref,
            "Greq_eV": Greq,
            "tau_local_MPa": tau_local,
            "tau_app_MPa": tau_app,
            "G_eV": Gval,
            "vstar_eV_per_MPa": vloc,
            "peak_metric_tau_v_over_kT": peak_metric,
            "solve_status": status,
        })
    return pd.DataFrame(rows)

def find_peak(df):
    d = df[np.isfinite(df["tau_app_MPa"]) & (df["tau_app_MPa"] > 0)].copy()
    if len(d) == 0:
        return None
    i = d["tau_app_MPa"].idxmax()
    return d.loc[i].to_dict()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results/analytic_peak_vs_barrier_scale")
    ap.add_argument("--T-list", default="1100,1200,1300")
    ap.add_argument("--edot", type=float, default=100.0)
    ap.add_argument("--rho-min", type=float, default=1e11)
    ap.add_argument("--rho-max", type=float, default=1e18)
    ap.add_argument("--n-rho", type=int, default=500)
    ap.add_argument("--scales", default="0.25,0.35,0.5,0.7,1.0,1.4,2.0")
    ap.add_argument("--floor-frac", type=float, default=0.0)
    ap.add_argument("--taylor-power", type=float, default=4.0)
    ap.add_argument("--nu0", type=float, default=1e12)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    T_list = [float(x) for x in args.T_list.split(",")]
    scales = [float(x) for x in args.scales.split(",")]
    rho_grid = np.logspace(np.log10(args.rho_min), np.log10(args.rho_max), args.n_rho)

    all_curves = []
    peaks = []

    for T in T_list:
        for scale in scales:
            df = analytic_curve(
                T=T,
                edot=args.edot,
                rho_grid=rho_grid,
                scale=scale,
                floor_frac=args.floor_frac,
                p=args.taylor_power,
                nu0=args.nu0,
            )
            all_curves.append(df)
            pk = find_peak(df)
            if pk is not None:
                peaks.append(pk)

    curves = pd.concat(all_curves, ignore_index=True)
    peak_df = pd.DataFrame(peaks)

    curves.to_csv(outdir/"analytic_density_strength_curves.csv", index=False)
    peak_df.to_csv(outdir/"analytic_peak_table.csv", index=False)

    print("\nAnalytical peak table:")
    if len(peak_df):
        cols = [
            "T_K", "edot_s", "barrier_scale", "rho_m2",
            "tau_app_MPa", "tau_local_MPa",
            "G_eV", "peak_metric_tau_v_over_kT",
            "solve_status",
        ]
        print(peak_df[cols].sort_values(["T_K","barrier_scale"]).to_string(index=False))
    else:
        print("No finite peaks found.")

    # Plot density-strength curves.
    for T in T_list:
        fig, ax = plt.subplots(figsize=(8.0, 5.6))
        for scale in scales:
            d = curves[(curves["T_K"] == T) & (curves["barrier_scale"] == scale)]
            good = np.isfinite(d["tau_app_MPa"]) & (d["tau_app_MPa"] > 0)
            if good.sum() == 0:
                continue
            ax.plot(d.loc[good, "rho_m2"], d.loc[good, "tau_app_MPa"],
                    linewidth=2, label=f"scale={scale:g}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"density $\rho$ (m$^{-2}$)")
        ax.set_ylabel("analytical applied stress (MPa)")
        ax.set_title(f"Analytical Arrhenius-Taylor curve, T={T:g} K, rate={args.edot:g} s$^{{-1}}$")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir/f"analytic_density_strength_T{int(T)}.png", dpi=300)
        plt.close(fig)

    # Plot peak density vs scale.
    if len(peak_df):
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        for T, g in peak_df.groupby("T_K"):
            g = g.sort_values("barrier_scale")
            ax.plot(g["barrier_scale"], g["rho_m2"], marker="o", linewidth=2, label=f"T={T:g} K")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("crossing barrier scale")
        ax.set_ylabel(r"peak density $\rho_{\rm peak}$ (m$^{-2}$)")
        ax.axhspan(1e15, 3e16, alpha=0.15, label="target DDD window")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir/"peak_density_vs_barrier_scale.png", dpi=300)
        plt.close(fig)

    print(f"\nwrote: {outdir/'analytic_peak_table.csv'}")
    print(f"plots written to: {outdir}")

if __name__ == "__main__":
    main()
