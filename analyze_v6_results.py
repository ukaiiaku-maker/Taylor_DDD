#!/usr/bin/env python3
"""
Analyze/plot clean_arrhenius_taylor_ddd_v6.py overnight sweep results.

Expected layout from run_v6_overnight.sh:
    results/v6_overnight/bs_on/T900_rho1e13/...
    results/v6_overnight/bs_off/T1100_rho1e13/...

Outputs:
    summary_v6.csv
    01_flow_stress_vs_rho.png
    02_epsp_ratio_vs_rho.png
    03_tau_local_vs_rho.png
    04_pin_cap_fraction_vs_rho.png
    05_backstress_ab_T1100.png
    06_representative_stress_strain.png
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_float_from_name(name: str, key: str) -> float:
    m = re.search(rf"{key}([0-9.eE+-]+)", name)
    return float(m.group(1)) if m else np.nan


def q(x, p):
    x = pd.Series(x).replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.quantile(p)) if len(x) else np.nan


def med(x):
    x = pd.Series(x).replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.median()) if len(x) else np.nan


def safe_last(df: pd.DataFrame, col: str, default=np.nan):
    return float(df[col].iloc[-1]) if col in df and len(df) else default


def load_params(run_dir: Path) -> dict:
    p = run_dir / "clean_arrhenius_params.json"
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}


def summarize_run(run_dir: Path, tail_fraction: float) -> dict | None:
    hist_file = run_dir / "single_glider_history.csv"
    if not hist_file.exists():
        return None

    try:
        h = pd.read_csv(hist_file)
    except Exception as e:
        print(f"Could not read {hist_file}: {e}")
        return None
    if h.empty:
        return None

    p = load_params(run_dir)

    # Prefer JSON parameters; fall back to directory names.
    T = float(p.get("temperature_K", parse_float_from_name(run_dir.name, "T")))
    rho = float(p.get("forest_rho_m2", parse_float_from_name(run_dir.name, "rho")))
    bs = str(p.get("backstress_mobility", "unknown"))
    for part in run_dir.parts[::-1]:
        if part in ("bs_on", "bs_off"):
            bs = part.replace("bs_", "")

    n_tail = max(5, int(np.ceil(len(h) * tail_fraction)))
    ht = h.tail(n_tail)

    tau_col = "tau_MPa" if "tau_MPa" in h else "sigma_MPa"
    eps_total = safe_last(h, "eps_total")
    eps_p = safe_last(h, "eps_plastic")
    eps_ratio = eps_p / eps_total if np.isfinite(eps_total) and abs(eps_total) > 0 else np.nan

    row = {
        "run_dir": str(run_dir),
        "backstress": bs,
        "T_K": T,
        "rho_m2": rho,
        "finished": (run_dir / "run.finished").exists(),
        "n_steps": int(safe_last(h, "step", len(h))),
        "n_rows": len(h),
        "eps_total_final": eps_total,
        "eps_plastic_final": eps_p,
        "epsp_over_epstotal_final": eps_ratio,

        # Signed stress tells you whether detailed-balance correction is over/under-shooting.
        "tau_tail_median_MPa": med(ht[tau_col]),
        "tau_tail_p10_MPa": q(ht[tau_col], 0.10),
        "tau_tail_p90_MPa": q(ht[tau_col], 0.90),
        "tau_tail_abs_median_MPa": med(np.abs(ht[tau_col])),
        "tau_final_MPa": safe_last(h, tau_col),
        "tau_max_abs_MPa": float(np.nanmax(np.abs(h[tau_col]))),

        "n_obstacles": safe_last(h, "n_obstacles_active"),
        "n_crossed_total_final": safe_last(h, "n_crossed_total", 0),
        "n_capture_total": float(h["n_capture"].sum()) if "n_capture" in h else np.nan,
        "n_depin_total": float(h["n_depin"].sum()) if "n_depin" in h else np.nan,
        "n_pinned_tail_median": med(ht["n_pinned_nodes"]) if "n_pinned_nodes" in h else np.nan,
        "n_live_pins_tail_median": med(ht["n_live_pins"]) if "n_live_pins" in h else np.nan,

        "tau_local_median_tail_MPa": med(ht["tau_local_median_MPa"]) if "tau_local_median_MPa" in h else np.nan,
        "tau_local_p90_tail_MPa": med(ht["tau_local_p90_MPa"]) if "tau_local_p90_MPa" in h else np.nan,
        "tau_local_max_tail_MPa": med(ht["tau_local_max_MPa"]) if "tau_local_max_MPa" in h else np.nan,
        "tau_local_uncapped_max_tail_MPa": med(ht["tau_local_uncapped_max_MPa"]) if "tau_local_uncapped_max_MPa" in h else np.nan,
        "frac_tau_local_capped_tail": med(ht["frac_tau_local_capped"]) if "frac_tau_local_capped" in h else np.nan,
        "phi_median_tail": med(ht["phi_median"]) if "phi_median" in h else np.nan,
        "phi_p90_tail": med(ht["phi_p90"]) if "phi_p90" in h else np.nan,
        "F_line_tension_median_tail_N": med(ht["F_line_tension_median_N"]) if "F_line_tension_median_N" in h else np.nan,
        "crossing_rate_max_tail_s": med(ht["crossing_rate_max_s"]) if "crossing_rate_max_s" in h else np.nan,
        "crossing_expected_events_tail": med(ht["crossing_expected_events_step"]) if "crossing_expected_events_step" in h else np.nan,
        "walltime_total_s": float(h["step_walltime_s"].sum()) if "step_walltime_s" in h else np.nan,
    }

    # Compact status flag for the main failure modes.
    flags = []
    if row["epsp_over_epstotal_final"] < 0.5:
        flags.append("not_flowing")
    if abs(row["tau_tail_median_MPa"]) < 1e-6 and row["epsp_over_epstotal_final"] > 1.2:
        flags.append("zero_stress_overshoot")
    if row["frac_tau_local_capped_tail"] > 0.05:
        flags.append("tau_local_cap_active")
    if row["n_live_pins_tail_median"] == 0 and row["n_obstacles"] > 0:
        flags.append("no_live_pins_tail")
    row["flags"] = ";".join(flags)
    return row


def collect(root: Path, tail_fraction: float) -> pd.DataFrame:
    rows = []
    for hist in sorted(root.rglob("single_glider_history.csv")):
        row = summarize_run(hist.parent, tail_fraction)
        if row:
            rows.append(row)
    if not rows:
        raise SystemExit(f"No single_glider_history.csv files found under {root}")
    return pd.DataFrame(rows).sort_values(["backstress", "T_K", "rho_m2"])


def plot_lines(df: pd.DataFrame, y: str, ylabel: str, title: str, outfile: Path):
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for (bs, T), g in df.groupby(["backstress", "T_K"], dropna=False):
        g = g.sort_values("rho_m2")
        ax.plot(g["rho_m2"], g[y], marker="o", label=f"bs={bs}, T={T:g} K")
    ax.set_xscale("log")
    ax.set_xlabel(r"Forest density $\rho$ (m$^{-2}$)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=220)
    plt.close(fig)


def plot_ab(df: pd.DataFrame, outfile: Path, T_ab: float = 1100.0):
    g = df[np.isclose(df["T_K"], T_ab)].copy()
    if g.empty or g["backstress"].nunique() < 2:
        return
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for bs, gg in g.groupby("backstress"):
        gg = gg.sort_values("rho_m2")
        ax.plot(gg["rho_m2"], gg["tau_tail_abs_median_MPa"], marker="o", label=f"bs={bs}")
    ax.set_xscale("log")
    ax.set_xlabel(r"Forest density $\rho$ (m$^{-2}$)")
    ax.set_ylabel("Tail median |tau| (MPa)")
    ax.set_title(f"Backstress A/B control at {T_ab:g} K")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=220)
    plt.close(fig)


def plot_representative_curves(df: pd.DataFrame, outfile: Path, max_curves: int = 12):
    d = df[df["backstress"].eq("on")].copy()
    if d.empty:
        d = df.copy()
    d = d.sort_values(["T_K", "rho_m2"])
    if len(d) > max_curves:
        idx = np.linspace(0, len(d) - 1, max_curves).round().astype(int)
        d = d.iloc[idx]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for _, r in d.iterrows():
        h = pd.read_csv(Path(r["run_dir"]) / "single_glider_history.csv")
        tau_col = "tau_MPa" if "tau_MPa" in h else "sigma_MPa"
        label = f"T{r['T_K']:g}, rho={r['rho_m2']:.0e}, bs={r['backstress']}"
        ax.plot(h["eps_total"], h[tau_col], lw=1.2, label=label)
    ax.set_xlabel("Total strain")
    ax.set_ylabel("Signed tau (MPa)")
    ax.set_title("Representative stress-strain histories")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outfile, dpi=220)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/v6_overnight",
                    help="Root result directory to scan.")
    ap.add_argument("--outdir", default=None,
                    help="Directory for summary CSV and plots. Default: <root>/analysis")
    ap.add_argument("--tail-fraction", type=float, default=0.5,
                    help="Fraction of final history rows used for plateau/tail stats.")
    ap.add_argument("--show-table", action="store_true",
                    help="Print a compact summary table.")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else root / "analysis"
    outdir.mkdir(parents=True, exist_ok=True)

    df = collect(root, args.tail_fraction)
    df.to_csv(outdir / "summary_v6.csv", index=False)

    plot_lines(df, "tau_tail_abs_median_MPa", "Tail median |tau| (MPa)",
               "Flow stress magnitude vs forest density",
               outdir / "01_flow_stress_vs_rho.png")

    plot_lines(df, "epsp_over_epstotal_final", r"Final $\epsilon_p/\epsilon_{total}$",
               "Plastic strain fraction vs forest density",
               outdir / "02_epsp_ratio_vs_rho.png")

    plot_lines(df, "tau_local_median_tail_MPa", "Tail median local pin stress (MPa)",
               "Pin local stress vs forest density",
               outdir / "03_tau_local_vs_rho.png")

    plot_lines(df, "frac_tau_local_capped_tail", "Tail fraction of live pins capped",
               "Tau-local cap activity vs forest density",
               outdir / "04_pin_cap_fraction_vs_rho.png")

    plot_ab(df, outdir / "05_backstress_ab_T1100.png", T_ab=1100.0)
    plot_representative_curves(df, outdir / "06_representative_stress_strain.png")

    compact_cols = [
        "backstress", "T_K", "rho_m2", "tau_tail_abs_median_MPa",
        "tau_tail_median_MPa", "epsp_over_epstotal_final",
        "n_crossed_total_final", "n_live_pins_tail_median",
        "tau_local_median_tail_MPa", "frac_tau_local_capped_tail", "flags"
    ]
    if args.show_table:
        with pd.option_context("display.max_rows", 200, "display.width", 180):
            print(df[compact_cols].to_string(index=False))

    print(f"Wrote: {outdir / 'summary_v6.csv'}")
    print(f"Wrote plots in: {outdir}")


if __name__ == "__main__":
    main()
