#!/usr/bin/env python3
"""
audit_ddd_event_strain_and_state.py

Audit whether forest crossing/depin events are coupled consistently to swept-area
plastic strain in clean_arrhenius_taylor_ddd_v12-style runs, and optionally plot
final line/obstacle states when final node files are present.

Usage examples:
  python3 audit_ddd_event_strain_and_state.py \
    results/v12_crossForceFactor_1p0_T1100 \
    results/v12_crossForceFactor_0p5_T1100 \
    results/v12_crossForceFactor_0p25_T1100 \
    --out audit_cross_force_factor --plot-final-states

  python3 audit_ddd_event_strain_and_state.py results/v12_crossForceFactor_0p25_T1100 \
    --plot-final-states --max-state-plots 20
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def safe_float(x: Any, default: float = float("nan")) -> float:
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def parse_from_path(run_dir: Path) -> Dict[str, Any]:
    s = str(run_dir)
    out: Dict[str, Any] = {}
    m = re.search(r"T([0-9]+(?:p[0-9]+)?)", s)
    if m:
        out["T_K_path"] = safe_float(m.group(1).replace("p", "."))
    m = re.search(r"rho([0-9]+(?:p[0-9]+)?e[+-]?[0-9]+)", s)
    if m:
        out["rho_m2_path"] = safe_float(m.group(1).replace("p", "."))
    m = re.search(r"bs_([A-Za-z0-9]+)", s)
    if m:
        out["backstress_path"] = m.group(1)
    return out


def find_run_dirs(roots: Iterable[Path]) -> List[Path]:
    dirs = []
    for root in roots:
        if root.is_file() and root.name == "single_glider_history.csv":
            dirs.append(root.parent)
        elif (root / "single_glider_history.csv").exists():
            dirs.append(root)
        else:
            dirs.extend(p.parent for p in root.rglob("single_glider_history.csv"))
    # de-duplicate preserving order
    seen = set()
    unique = []
    for d in dirs:
        rp = str(d.resolve())
        if rp not in seen:
            seen.add(rp)
            unique.append(d)
    return unique


def expected_mean_node_dx_from_epsp(epsp: float, params: Dict[str, Any]) -> float:
    """Mean accumulated reduced x displacement per node inferred from eps_p.

    Code relation: d_eps = sum(dx_nodes)/(Lx*nn) * b/s_out.
    Therefore mean dx per node over nline*nn nodes is:
        eps * Lx * nn * s_out/b / (nline*nn) = eps*Lx*s_out/(b*nline)
    """
    Lx = safe_float(params.get("cell_lx_reduced"), 3060.0)
    b = safe_float(params.get("b_m"), 2.48e-10)
    nline = safe_float(params.get("mobile_line_count"), 4.0)
    s_out = safe_float(params.get("out_of_plane_spacing_m_used"), params.get("out_of_plane_spacing_m", 1e-6))
    if not np.isfinite(epsp) or b <= 0 or nline <= 0:
        return float("nan")
    return epsp * Lx * s_out / (b * nline)


def summarize_run(run_dir: Path) -> Dict[str, Any]:
    hist_path = run_dir / "single_glider_history.csv"
    ev_path = run_dir / "single_glider_crossing_events.csv"
    params_path = run_dir / "clean_arrhenius_params.json"
    params = load_json(params_path)
    path_info = parse_from_path(run_dir)

    h = pd.read_csv(hist_path)
    if len(h) == 0:
        raise ValueError(f"empty history: {hist_path}")
    last = h.iloc[-1]

    epsp_final = safe_float(last.get("eps_plastic"))
    epst_final = safe_float(last.get("eps_total"))
    n_cross_final = safe_float(last.get("n_crossed_total"), 0.0)
    n_depin_sum = safe_float(h.get("n_depin", pd.Series(dtype=float)).sum(), 0.0)
    n_capture_sum = safe_float(h.get("n_capture", pd.Series(dtype=float)).sum(), 0.0)
    n_candidate_sum = safe_float(h.get("n_candidate_tests", pd.Series(dtype=float)).sum(), 0.0)
    d_eps_pos = safe_float(h.get("d_eps_p", pd.Series(dtype=float)).clip(lower=0).sum(), float("nan")) if "d_eps_p" in h else float("nan")
    d_eps_abs = safe_float(h.get("d_eps_p", pd.Series(dtype=float)).abs().sum(), float("nan")) if "d_eps_p" in h else float("nan")

    depin_per_epsp = n_depin_sum / epsp_final if epsp_final and epsp_final > 0 else float("nan")
    capture_per_epsp = n_capture_sum / epsp_final if epsp_final and epsp_final > 0 else float("nan")
    epsp_per_depin = epsp_final / n_depin_sum if n_depin_sum > 0 else float("nan")
    epsp_per_capture = epsp_final / n_capture_sum if n_capture_sum > 0 else float("nan")

    summary: Dict[str, Any] = {
        "run_dir": str(run_dir),
        "root_name": run_dir.parents[2].name if len(run_dir.parents) >= 3 else run_dir.parent.name,
        "T_K": safe_float(params.get("temperature_K"), path_info.get("T_K_path", float("nan"))),
        "rho_m2": safe_float(params.get("forest_rho_m2"), path_info.get("rho_m2_path", float("nan"))),
        "backstress": params.get("backstress_mobility", path_info.get("backstress_path", "")),
        "cross_force_scale_factor": safe_float(params.get("cross_force_scale_factor", params.get("v11_cross_force_scale_factor")), float("nan")),
        "cross_force_scale_N": safe_float(params.get("v11_cross_force_scale_N"), float("nan")),
        "crossing_drive_mode": params.get("crossing_drive_mode", params.get("v11_crossing_drive_mode", "")),
        "plastic_strain_source": params.get("plastic_strain_source", params.get("v3_plastic_strain_source", "")),
        "capture_mode": params.get("capture_mode", ""),
        "rearm_radius_reduced": safe_float(params.get("v7_rearm_radius_reduced", params.get("rearm_radius_reduced")), float("nan")),
        "capture_radius_reduced": safe_float(params.get("capture_radius_reduced"), float("nan")),
        "Lx_reduced": safe_float(params.get("cell_lx_reduced"), float("nan")),
        "Lz_reduced": safe_float(params.get("cell_lz_reduced"), float("nan")),
        "nline": safe_float(params.get("mobile_line_count"), float("nan")),
        "nnodes_per_line": safe_float(params.get("mobile_line_nodes"), float("nan")),
        "nsteps": int(safe_float(last.get("step"), 0)),
        "eps_total_final": epst_final,
        "eps_plastic_final": epsp_final,
        "epsp_over_epstotal_final": epsp_final / epst_final if epst_final and epst_final > 0 else float("nan"),
        "tau_final_MPa": safe_float(last.get("tau_MPa", last.get("sigma_MPa"))),
        "tau_tail_median_MPa": safe_float(h["tau_MPa"].tail(max(1, len(h)//5)).median()) if "tau_MPa" in h else float("nan"),
        "n_crossed_total_final": n_cross_final,
        "n_depin_sum_history": n_depin_sum,
        "n_capture_sum_history": n_capture_sum,
        "n_candidate_tests_sum": n_candidate_sum,
        "depin_per_capture": n_depin_sum / n_capture_sum if n_capture_sum > 0 else float("nan"),
        "capture_per_candidate": n_capture_sum / n_candidate_sum if n_candidate_sum > 0 else float("nan"),
        "depin_per_epsp": depin_per_epsp,
        "capture_per_epsp": capture_per_epsp,
        "epsp_per_depin": epsp_per_depin,
        "epsp_per_capture": epsp_per_capture,
        "positive_d_eps_p_sum": d_eps_pos,
        "abs_d_eps_p_sum": d_eps_abs,
        "mean_node_dx_from_epsp_reduced": expected_mean_node_dx_from_epsp(epsp_final, params),
        "mean_node_cell_wraps_from_epsp": expected_mean_node_dx_from_epsp(epsp_final, params) / safe_float(params.get("cell_lx_reduced"), 3060.0),
        "n_live_pins_tail_median": safe_float(h.get("n_live_pins", pd.Series(dtype=float)).tail(max(1, len(h)//5)).median(), float("nan")),
        "n_pinned_nodes_tail_median": safe_float(h.get("n_pinned_nodes", pd.Series(dtype=float)).tail(max(1, len(h)//5)).median(), float("nan")),
        "n_pinned_raw_lines_tail_median": safe_float(h.get("n_pinned_raw_lines", pd.Series(dtype=float)).tail(max(1, len(h)//5)).median(), float("nan")),
        "mean_dx_after_projection_tail": safe_float(h.get("mean_dx_after_projection_free_reduced", pd.Series(dtype=float)).tail(max(1, len(h)//5)).median(), float("nan")),
        "dx_app_tail": safe_float(h.get("dx_app_reduced", pd.Series(dtype=float)).tail(max(1, len(h)//5)).median(), float("nan")),
    }

    if ev_path.exists() and ev_path.stat().st_size > 0:
        try:
            ev = pd.read_csv(ev_path)
            dep = ev[ev.get("event", "") == "depin_cross"] if "event" in ev else pd.DataFrame()
            cap = ev[ev.get("event", "") == "capture_pin"] if "event" in ev else pd.DataFrame()
            summary["event_rows"] = len(ev)
            summary["depin_event_rows"] = len(dep)
            summary["capture_event_rows"] = len(cap)
            if len(dep) and {"line_id", "obs_id"}.issubset(dep.columns):
                gp = dep.groupby(["line_id", "obs_id"]).size()
                summary["unique_depin_line_obs_pairs"] = int(len(gp))
                summary["mean_depin_per_unique_pair"] = float(gp.mean())
                summary["median_depin_per_unique_pair"] = float(gp.median())
                summary["max_depin_per_unique_pair"] = int(gp.max())
                summary["frac_depin_in_repeated_pairs"] = float(gp[gp > 1].sum() / gp.sum()) if gp.sum() > 0 else float("nan")
            if len(cap) and {"line_id", "obs_id"}.issubset(cap.columns):
                gp = cap.groupby(["line_id", "obs_id"]).size()
                summary["unique_capture_line_obs_pairs"] = int(len(gp))
                summary["mean_capture_per_unique_pair"] = float(gp.mean())
                summary["max_capture_per_unique_pair"] = int(gp.max())
        except Exception as e:
            summary["event_read_error"] = str(e)
    return summary


def plot_scalar_time(run_dir: Path, outdir: Path) -> None:
    h = pd.read_csv(run_dir / "single_glider_history.csv")
    rho = "unknown"
    params = load_json(run_dir / "clean_arrhenius_params.json")
    if params.get("forest_rho_m2") is not None:
        rho = f"rho{params['forest_rho_m2']:.0e}"
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{run_dir.parent.parent.name}_{run_dir.parent.name}_{run_dir.name}_{rho}")

    xcol = "eps_plastic" if "eps_plastic" in h.columns else "step"
    for cols, suffix, ylabel in [
        (["n_crossed_total", "n_capture", "n_depin"], "events_vs_epsp", "events / step counts"),
        (["d_eps_p"], "depsp_vs_epsp", "d_eps_p"),
        (["n_live_pins", "n_pinned_nodes"], "pins_vs_epsp", "pins"),
        (["mean_dx_after_projection_free_reduced", "dx_app_reduced", "mean_dx_raw_free_reduced"], "dx_vs_epsp", "dx reduced"),
    ]:
        have = [c for c in cols if c in h.columns]
        if not have:
            continue
        plt.figure(figsize=(7, 4.5))
        for c in have:
            y = h[c]
            if c in ["n_capture", "n_depin"]:
                y = y.cumsum()
                label = "cum_" + c
            else:
                label = c
            plt.plot(h[xcol], y, label=label, linewidth=1.2)
        plt.xlabel(xcol)
        plt.ylabel(ylabel)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / f"{tag}_{suffix}.png", dpi=200)
        plt.close()


def plot_final_state(run_dir: Path, outdir: Path, max_obstacles: int = 5000) -> Optional[Path]:
    nodes_path = run_dir / "single_glider_final_nodes.csv"
    obs_path = run_dir / "fixed_forest_obstacles.csv"
    if not nodes_path.exists():
        return None
    nodes = pd.read_csv(nodes_path)
    obs = pd.read_csv(obs_path) if obs_path.exists() else pd.DataFrame()
    params = load_json(run_dir / "clean_arrhenius_params.json")
    Lx = safe_float(params.get("cell_lx_reduced"), nodes["x_reduced"].max() if "x_reduced" in nodes else 1)
    rho = safe_float(params.get("forest_rho_m2"), float("nan"))
    fac = safe_float(params.get("v11_cross_force_scale_factor", params.get("cross_force_scale_factor")), float("nan"))

    plt.figure(figsize=(9, 4.8))
    if len(obs):
        if len(obs) > max_obstacles:
            obs_plot = obs.sample(max_obstacles, random_state=1)
        else:
            obs_plot = obs
        plt.scatter(obs_plot["x_reduced"], obs_plot["z_reduced"], s=2, alpha=0.25, label="forest obstacles")

    for li, g in nodes.groupby("line_id"):
        g = g.sort_values("node_id")
        x = g["x_reduced"].to_numpy(dtype=float)
        z = g["z_reduced"].to_numpy(dtype=float)
        # Draw broken line segments to avoid spurious lines across periodic x jumps.
        start = 0
        for k in range(1, len(x)):
            dx = (x[k] - x[k-1] + 0.5 * Lx) % Lx - 0.5 * Lx
            if abs(dx) > 0.4 * Lx:
                if k - start > 1:
                    plt.plot(x[start:k], z[start:k], linewidth=0.8)
                start = k
        if len(x) - start > 1:
            plt.plot(x[start:], z[start:], linewidth=0.8)
        if "pinned" in g.columns:
            pg = g[g["pinned"].astype(int) == 1]
            if len(pg):
                plt.scatter(pg["x_reduced"], pg["z_reduced"], s=18, marker="x", label="pinned nodes" if li == 0 else None)

    plt.xlabel("x / b")
    plt.ylabel("z / b")
    plt.title(f"Final state: rho={rho:.2e}, force factor={fac:g}\n{run_dir}", fontsize=9)
    plt.tight_layout()
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"final_state_{run_dir.parent.parent.name}_{run_dir.parent.name}_{run_dir.name}")
    path = outdir / f"{tag}.png"
    plt.savefig(path, dpi=220)
    plt.close()
    return path


def make_cross_root_plots(df: pd.DataFrame, outdir: Path) -> None:
    if df.empty:
        return
    df = df.sort_values(["cross_force_scale_factor", "rho_m2"])
    for y, ylabel, fname in [
        ("depin_per_epsp", "depin events / eps_p", "depin_per_epsp_vs_rho.png"),
        ("epsp_per_depin", "eps_p / depin event", "epsp_per_depin_vs_rho.png"),
        ("n_depin_sum_history", "total depin events", "depin_total_vs_rho.png"),
        ("n_capture_sum_history", "total capture events", "capture_total_vs_rho.png"),
        ("tau_tail_median_MPa", "tail tau median (MPa)", "tau_tail_vs_rho.png"),
        ("epsp_over_epstotal_final", "eps_p / eps_total", "epsp_ratio_vs_rho.png"),
        ("mean_depin_per_unique_pair", "mean depins per unique line-obstacle pair", "repeat_depin_pairs_vs_rho.png"),
        ("frac_depin_in_repeated_pairs", "fraction depins in repeated pairs", "repeat_fraction_vs_rho.png"),
        ("mean_node_cell_wraps_from_epsp", "mean node cell wraps inferred from eps_p", "mean_wraps_vs_rho.png"),
    ]:
        if y not in df.columns or df[y].dropna().empty:
            continue
        plt.figure(figsize=(7, 4.8))
        group_col = "cross_force_scale_factor" if "cross_force_scale_factor" in df.columns else "root_name"
        for key, g in df.groupby(group_col):
            gg = g.dropna(subset=["rho_m2", y]).sort_values("rho_m2")
            if len(gg) == 0:
                continue
            plt.plot(gg["rho_m2"], gg[y], marker="o", label=str(key))
        plt.xscale("log")
        if y not in ["epsp_over_epstotal_final", "frac_depin_in_repeated_pairs"]:
            vals = df[y].dropna()
            if len(vals) and (vals > 0).all() and vals.max() / max(vals.min(), 1e-300) > 20:
                plt.yscale("log")
        plt.xlabel("forest rho (m$^{-2}$)")
        plt.ylabel(ylabel)
        plt.legend(title=group_col, fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / fname, dpi=220)
        plt.close()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("roots", nargs="+", help="Run roots or run directories containing single_glider_history.csv")
    ap.add_argument("--out", default="audit_ddd_event_strain", help="Output directory")
    ap.add_argument("--plot-final-states", action="store_true")
    ap.add_argument("--max-state-plots", type=int, default=30)
    ap.add_argument("--plot-time-series", action="store_true")
    args = ap.parse_args()

    roots = [Path(r) for r in args.roots]
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    run_dirs = find_run_dirs(roots)
    if not run_dirs:
        raise SystemExit("No single_glider_history.csv files found.")

    rows: List[Dict[str, Any]] = []
    for rd in run_dirs:
        try:
            rows.append(summarize_run(rd))
        except Exception as e:
            rows.append({"run_dir": str(rd), "error": str(e)})

    df = pd.DataFrame(rows)
    if "rho_m2" in df.columns:
        df = df.sort_values([c for c in ["cross_force_scale_factor", "T_K", "rho_m2"] if c in df.columns])
    csv_path = outdir / "audit_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path}")
    with pd.option_context("display.max_columns", 999, "display.width", 240):
        cols = [c for c in [
            "cross_force_scale_factor", "rho_m2", "tau_tail_median_MPa", "epsp_over_epstotal_final",
            "n_depin_sum_history", "n_capture_sum_history", "depin_per_epsp", "epsp_per_depin",
            "unique_depin_line_obs_pairs", "mean_depin_per_unique_pair", "max_depin_per_unique_pair",
            "frac_depin_in_repeated_pairs", "mean_node_cell_wraps_from_epsp", "run_dir"
        ] if c in df.columns]
        print(df[cols].to_string(index=False))

    make_cross_root_plots(df, outdir)

    if args.plot_time_series:
        ts_dir = outdir / "time_series"
        ts_dir.mkdir(exist_ok=True)
        for rd in run_dirs:
            try:
                plot_scalar_time(rd, ts_dir)
            except Exception as e:
                print(f"time-series plot failed for {rd}: {e}")

    if args.plot_final_states:
        st_dir = outdir / "final_states"
        st_dir.mkdir(exist_ok=True)
        n = 0
        for rd in run_dirs:
            if n >= args.max_state_plots:
                break
            try:
                p = plot_final_state(rd, st_dir)
                if p is not None:
                    print(f"Wrote {p}")
                    n += 1
            except Exception as e:
                print(f"state plot failed for {rd}: {e}")

    print(f"Wrote plots in {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
