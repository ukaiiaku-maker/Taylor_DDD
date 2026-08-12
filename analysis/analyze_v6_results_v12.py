#!/usr/bin/env python3
"""
analyze_v6_results_v12.py

Drop-in analyzer for the reduced OpenDiS Arrhenius Taylor / Peierls DDD runs.
It preserves the old v6-style summary quantities and adds v12-specific
center-of-mass projection diagnostics:

  mean_dx_raw_free_reduced
  mean_dx_after_projection_free_reduced
  dx_app_reduced
  n_com_projected_lines
  n_pinned_raw_lines

The main purpose is to distinguish:
  v12 default: pinned lines are left raw, so n_pinned_raw_lines > 0 in the tail
  compatibility mode: --project-backstress-on-pinned-lines, so pinned lines are
                      also COM projected and n_pinned_raw_lines should stay ~0

Usage:
  python3 analyze_v6_results_v12.py --root results/v12_pinnedLinesRaw_T1100_AB --show-table
  python3 analyze_v6_results_v12.py --root results/v12_pinnedLinesRaw_T1100_AB \
      results/v12_projectPinnedLines_T1100_AB --out summary_v12_AB.csv --show-table

You can also copy this over analyze_v6_results.py if you want to keep the old
command name.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


HISTORY_NAME = "single_glider_history.csv"
PARAMS_NAME = "clean_arrhenius_params.json"
FINISHED_NAME = "run.finished"


def _safe_float(x: Any, default: float = float("nan")) -> float:
    try:
        if x is None:
            return default
        y = float(x)
        return y if math.isfinite(y) else default
    except Exception:
        return default


def _safe_bool(x: Any, default: bool = False) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    if s in {"1", "true", "yes", "on"}:
        return True
    if s in {"0", "false", "no", "off"}:
        return False
    return default


def _nan_if_missing(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series(np.nan, index=df.index, dtype=float)


def _tail(df: pd.DataFrame, tail_frac: float, tail_min_rows: int) -> pd.DataFrame:
    if len(df) == 0:
        return df
    n_tail = max(int(math.ceil(len(df) * float(tail_frac))), int(tail_min_rows))
    n_tail = min(max(n_tail, 1), len(df))
    return df.tail(n_tail)


def _q(s: pd.Series, q: float) -> float:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if len(s) == 0:
        return float("nan")
    return float(s.quantile(q))


def _median(s: pd.Series) -> float:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if len(s) == 0:
        return float("nan")
    return float(s.median())


def _max(s: pd.Series) -> float:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if len(s) == 0:
        return float("nan")
    return float(s.max())


def _sum(s: pd.Series) -> float:
    s = pd.to_numeric(s, errors="coerce").dropna()
    if len(s) == 0:
        return 0.0
    return float(s.sum())


def _last(s: pd.Series, default: float = float("nan")) -> float:
    if s is None or len(s) == 0:
        return default
    ss = pd.to_numeric(s, errors="coerce").dropna()
    if len(ss) == 0:
        return default
    return float(ss.iloc[-1])


def _read_params(run_dir: Path) -> Dict[str, Any]:
    p = run_dir / PARAMS_NAME
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}


def _parse_rho_from_path(path: Path) -> float:
    m = re.search(r"rho([0-9.eE+-]+)", str(path))
    return _safe_float(m.group(1)) if m else float("nan")


def _parse_T_from_path(path: Path) -> float:
    m = re.search(r"T([0-9.]+)", str(path))
    return _safe_float(m.group(1)) if m else float("nan")


def _parse_backstress_from_path(path: Path) -> str:
    s = str(path)
    if re.search(r"(^|[/_])bs[_-]?on($|[/_])", s):
        return "on"
    if re.search(r"(^|[/_])bs[_-]?off($|[/_])", s):
        return "off"
    return ""


def _flag_join(flags: Iterable[str]) -> str:
    flags = [f for f in flags if f]
    return ";".join(flags) if flags else ""


def find_history_files(roots: List[Path]) -> List[Path]:
    files: List[Path] = []
    for root in roots:
        if root.is_file() and root.name == HISTORY_NAME:
            files.append(root)
        elif root.is_dir():
            files.extend(root.rglob(HISTORY_NAME))
    # Unique, stable order.
    seen = set()
    out = []
    for f in sorted(files, key=lambda p: str(p)):
        rp = f.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(f)
    return out


def summarize_run(hist_path: Path, tail_frac: float, tail_min_rows: int) -> Dict[str, Any]:
    run_dir = hist_path.parent
    params = _read_params(run_dir)
    finished = (run_dir / FINISHED_NAME).exists()

    # Low-memory enough for these runs; keep simple and robust.
    df = pd.read_csv(hist_path)
    tail = _tail(df, tail_frac=tail_frac, tail_min_rows=tail_min_rows)

    eps_total_final = _last(_nan_if_missing(df, "eps_total"))
    eps_plastic_final = _last(_nan_if_missing(df, "eps_plastic"))
    epsp_over_epstotal = eps_plastic_final / eps_total_final if eps_total_final and math.isfinite(eps_total_final) else float("nan")

    tau = _nan_if_missing(df, "tau_MPa") if "tau_MPa" in df.columns else _nan_if_missing(df, "sigma_MPa")
    tau_tail = _nan_if_missing(tail, "tau_MPa") if "tau_MPa" in tail.columns else _nan_if_missing(tail, "sigma_MPa")

    # v12 displacement/projection diagnostics.
    dx_raw_tail = _nan_if_missing(tail, "mean_dx_raw_free_reduced")
    dx_after_tail = _nan_if_missing(tail, "mean_dx_after_projection_free_reduced")
    dx_app_tail = _nan_if_missing(tail, "dx_app_reduced")
    n_projected_tail = _nan_if_missing(tail, "n_com_projected_lines")
    n_pinned_raw_tail = _nan_if_missing(tail, "n_pinned_raw_lines")

    mean_dx_raw_tail = _median(dx_raw_tail)
    mean_dx_after_tail = _median(dx_after_tail)
    dx_app_tail_median = _median(dx_app_tail)
    dx_projection_deficit = dx_app_tail_median - mean_dx_after_tail if math.isfinite(dx_app_tail_median) and math.isfinite(mean_dx_after_tail) else float("nan")
    dx_projection_deficit_abs = abs(dx_projection_deficit) if math.isfinite(dx_projection_deficit) else float("nan")
    if math.isfinite(dx_app_tail_median) and abs(dx_app_tail_median) > 1e-300 and math.isfinite(dx_projection_deficit):
        dx_projection_deficit_over_dx_app = dx_projection_deficit / dx_app_tail_median
    else:
        dx_projection_deficit_over_dx_app = float("nan")

    n_pinned_raw_med = _median(n_pinned_raw_tail)
    n_pinned_raw_max = _max(n_pinned_raw_tail)
    n_projected_med = _median(n_projected_tail)
    n_projected_max = _max(n_projected_tail)

    v12_cols = [
        "mean_dx_raw_free_reduced",
        "mean_dx_after_projection_free_reduced",
        "dx_app_reduced",
        "n_com_projected_lines",
        "n_pinned_raw_lines",
    ]
    v12_fields_present = all(c in df.columns for c in v12_cols)
    v12_branch_active_tail = bool(math.isfinite(n_pinned_raw_max) and n_pinned_raw_max > 0.0)

    # Metadata: prefer params, fall back to history/path.
    T_K = _safe_float(params.get("temperature_K"), _last(_nan_if_missing(df, "T_K"), _parse_T_from_path(run_dir)))
    rho_m2 = _safe_float(params.get("forest_rho_m2"), _last(_nan_if_missing(df, "forest_rho_actual_m2"), _parse_rho_from_path(run_dir)))
    backstress = str(params.get("backstress_mobility", "") or _parse_backstress_from_path(run_dir))
    project_pinned = _safe_bool(params.get("v12_project_backstress_on_pinned_lines"), default=False)

    flags: List[str] = []
    if not finished:
        flags.append("not_finished")
    if not v12_fields_present:
        flags.append("v12_fields_missing")
    if _median(_nan_if_missing(tail, "n_live_pins")) <= 0:
        flags.append("no_live_pins_tail")
    if v12_fields_present and not project_pinned and not v12_branch_active_tail:
        flags.append("v12_branch_not_active_tail")
    if v12_fields_present and project_pinned and v12_branch_active_tail:
        flags.append("unexpected_pinned_raw_in_compat_mode")
    if math.isfinite(epsp_over_epstotal) and epsp_over_epstotal > 0.995:
        flags.append("epsp_near_total")
    if v12_branch_active_tail and math.isfinite(dx_projection_deficit_over_dx_app):
        if dx_projection_deficit_over_dx_app > 0.05:
            flags.append("dx_suppressed_vs_dx_app")
        elif abs(dx_projection_deficit_over_dx_app) < 0.02:
            flags.append("little_dx_suppression_vs_dx_app")

    row: Dict[str, Any] = {
        "run_dir": str(run_dir),
        "backstress": backstress,
        "T_K": T_K,
        "rho_m2": rho_m2,
        "finished": bool(finished),
        "n_steps": int(_last(_nan_if_missing(df, "step"), len(df))),
        "n_rows": int(len(df)),
        "tail_rows": int(len(tail)),
        "eps_total_final": eps_total_final,
        "eps_plastic_final": eps_plastic_final,
        "epsp_over_epstotal_final": epsp_over_epstotal,
        "plastic_deficit_final": 1.0 - epsp_over_epstotal if math.isfinite(epsp_over_epstotal) else float("nan"),
        "tau_tail_median_MPa": _median(tau_tail),
        "tau_tail_p10_MPa": _q(tau_tail, 0.10),
        "tau_tail_p90_MPa": _q(tau_tail, 0.90),
        "tau_tail_abs_median_MPa": _median(tau_tail.abs()),
        "tau_final_MPa": _last(tau),
        "tau_max_abs_MPa": _max(tau.abs()),
        "n_obstacles": _safe_float(params.get("n_obstacles"), _last(_nan_if_missing(df, "n_obstacles_active"))),
        "n_crossed_total_final": _last(_nan_if_missing(df, "n_crossed_total")),
        "n_capture_total": _sum(_nan_if_missing(df, "n_capture")),
        "n_depin_total": _sum(_nan_if_missing(df, "n_depin")),
        "n_pinned_tail_median": _median(_nan_if_missing(tail, "n_pinned_nodes")),
        "n_live_pins_tail_median": _median(_nan_if_missing(tail, "n_live_pins")),
        "tau_local_median_tail_MPa": _median(_nan_if_missing(tail, "tau_local_median_MPa")),
        "tau_local_p90_tail_MPa": _median(_nan_if_missing(tail, "tau_local_p90_MPa")),
        "tau_local_max_tail_MPa": _max(_nan_if_missing(tail, "tau_local_max_MPa")),
        "tau_local_uncapped_max_tail_MPa": _max(_nan_if_missing(tail, "tau_local_uncapped_max_MPa")),
        "frac_tau_local_capped_tail": _median(_nan_if_missing(tail, "frac_tau_local_capped")),
        "phi_median_tail": _median(_nan_if_missing(tail, "phi_median")),
        "phi_p90_tail": _median(_nan_if_missing(tail, "phi_p90")),
        "F_line_tension_median_tail_N": _median(_nan_if_missing(tail, "F_line_tension_median_N")),
        "crossing_rate_max_tail_s": _max(_nan_if_missing(tail, "crossing_rate_max_s")),
        "crossing_expected_events_tail": _sum(_nan_if_missing(tail, "crossing_expected_events_step")),
        "free_node_fraction_tail_median": _median(_nan_if_missing(tail, "free_node_fraction")),
        "pinned_node_fraction_tail_median": _median(_nan_if_missing(tail, "pinned_node_fraction")),
        "tau_back_abs_mean_tail_MPa": _median(_nan_if_missing(tail, "tau_back_abs_mean_MPa")),
        "tau_back_abs_p90_tail_MPa": _median(_nan_if_missing(tail, "tau_back_abs_p90_MPa")),
        "peierls_only_baseline_tau_MPa": _median(_nan_if_missing(tail, "peierls_only_baseline_tau_MPa")),
        "required_free_glide_net_rate_s": _median(_nan_if_missing(tail, "required_free_glide_net_rate_s")),
        "walltime_total_s": _sum(_nan_if_missing(df, "step_walltime_s")),
        # v12 COM projection diagnostics.
        "v12_project_backstress_on_pinned_lines": project_pinned,
        "v12_backstress_projection_rule": params.get("v12_backstress_projection_rule", ""),
        "v12_fields_present": v12_fields_present,
        "v12_branch_active_tail": v12_branch_active_tail,
        "mean_dx_raw_free_tail_reduced": mean_dx_raw_tail,
        "mean_dx_after_projection_tail_reduced": mean_dx_after_tail,
        "dx_app_tail_reduced": dx_app_tail_median,
        "dx_projection_deficit_tail_reduced": dx_projection_deficit,
        "dx_projection_deficit_abs_tail_reduced": dx_projection_deficit_abs,
        "dx_projection_deficit_over_dx_app_tail": dx_projection_deficit_over_dx_app,
        "n_com_projected_lines_tail_median": n_projected_med,
        "n_com_projected_lines_tail_max": n_projected_max,
        "n_pinned_raw_lines_tail_median": n_pinned_raw_med,
        "n_pinned_raw_lines_tail_max": n_pinned_raw_max,
        "mean_dx_raw_minus_after_tail_reduced": mean_dx_raw_tail - mean_dx_after_tail if math.isfinite(mean_dx_raw_tail) and math.isfinite(mean_dx_after_tail) else float("nan"),
        "flags": _flag_join(flags),
    }
    return row


def make_plots(summary: pd.DataFrame, outdir: Path) -> None:
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    if len(summary) == 0:
        return

    # Plot each (T, backstress, project flag) group separately to avoid mixing runs.
    group_cols = ["T_K", "backstress", "v12_project_backstress_on_pinned_lines"]
    for keys, g in summary.groupby(group_cols, dropna=False):
        g = g.sort_values("rho_m2")
        label = f"T{keys[0]:g}_bs{keys[1]}_projectPinned{keys[2]}".replace("/", "_")

        fig, ax = plt.subplots()
        ax.plot(g["rho_m2"], g["tau_tail_abs_median_MPa"], marker="o")
        ax.set_xscale("log")
        ax.set_xlabel(r"forest density $\\rho$ (m$^{-2}$)")
        ax.set_ylabel(r"tail median $|\\tau|$ (MPa)")
        ax.set_title(label)
        fig.tight_layout()
        fig.savefig(outdir / f"tau_vs_rho_{label}.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.plot(g["rho_m2"], g["epsp_over_epstotal_final"], marker="o")
        ax.set_xscale("log")
        ax.set_xlabel(r"forest density $\\rho$ (m$^{-2}$)")
        ax.set_ylabel(r"$\\epsilon_p/\\epsilon_{total}$ final")
        ax.set_title(label)
        fig.tight_layout()
        fig.savefig(outdir / f"epsp_ratio_vs_rho_{label}.png", dpi=200)
        plt.close(fig)

        if g["v12_fields_present"].any():
            fig, ax = plt.subplots()
            ax.plot(g["rho_m2"], g["n_pinned_raw_lines_tail_median"], marker="o")
            ax.set_xscale("log")
            ax.set_xlabel(r"forest density $\\rho$ (m$^{-2}$)")
            ax.set_ylabel("tail median pinned-raw lines")
            ax.set_title(label)
            fig.tight_layout()
            fig.savefig(outdir / f"pinned_raw_lines_vs_rho_{label}.png", dpi=200)
            plt.close(fig)

            fig, ax = plt.subplots()
            ax.plot(g["rho_m2"], g["dx_projection_deficit_over_dx_app_tail"], marker="o")
            ax.set_xscale("log")
            ax.set_xlabel(r"forest density $\\rho$ (m$^{-2}$)")
            ax.set_ylabel(r"$(dx_{app}-dx_{after})/dx_{app}$ tail")
            ax.set_title(label)
            fig.tight_layout()
            fig.savefig(outdir / f"dx_suppression_vs_rho_{label}.png", dpi=200)
            plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Summarize Arrhenius Taylor DDD results with v12 projection diagnostics.")
    ap.add_argument("--root", nargs="+", required=True, help="One or more result roots, run dirs, or history CSV files.")
    ap.add_argument("--out", default="", help="Output CSV path. Default: <root>/summary_v12.csv for one root, else ./summary_v12.csv")
    ap.add_argument("--tail-frac", type=float, default=0.20, help="Fraction of the run used for tail statistics.")
    ap.add_argument("--tail-min-rows", type=int, default=50, help="Minimum rows used for tail statistics, capped by run length.")
    ap.add_argument("--show-table", action="store_true", help="Print compact table to stdout.")
    ap.add_argument("--make-plots", action="store_true", help="Write simple diagnostic PNG plots next to the summary CSV.")
    ap.add_argument("--only-finished", action="store_true", help="Skip runs without run.finished.")
    args = ap.parse_args(argv)

    roots = [Path(r) for r in args.root]
    histories = find_history_files(roots)
    if not histories:
        raise SystemExit("No single_glider_history.csv files found under: " + ", ".join(map(str, roots)))

    rows = []
    errors = []
    for hist in histories:
        try:
            row = summarize_run(hist, tail_frac=args.tail_frac, tail_min_rows=args.tail_min_rows)
            if args.only_finished and not row.get("finished", False):
                continue
            rows.append(row)
        except Exception as exc:
            errors.append((str(hist), repr(exc)))

    summary = pd.DataFrame(rows)
    if len(summary):
        summary = summary.sort_values([
            "T_K",
            "backstress",
            "v12_project_backstress_on_pinned_lines",
            "rho_m2",
            "run_dir",
        ], na_position="last")

    if args.out:
        out = Path(args.out)
    elif len(roots) == 1 and roots[0].is_dir():
        out = roots[0] / "summary_v12.csv"
    else:
        out = Path("summary_v12.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, index=False)

    if args.show_table:
        cols = [
            "T_K",
            "rho_m2",
            "backstress",
            "v12_project_backstress_on_pinned_lines",
            "finished",
            "epsp_over_epstotal_final",
            "plastic_deficit_final",
            "tau_tail_abs_median_MPa",
            "n_live_pins_tail_median",
            "n_pinned_raw_lines_tail_median",
            "n_com_projected_lines_tail_median",
            "dx_app_tail_reduced",
            "mean_dx_after_projection_tail_reduced",
            "dx_projection_deficit_over_dx_app_tail",
            "flags",
        ]
        cols = [c for c in cols if c in summary.columns]
        if len(summary):
            with pd.option_context("display.max_rows", 200, "display.max_columns", 200, "display.width", 240):
                print(summary[cols].to_string(index=False))
        else:
            print("No runs summarized.")

    if args.make_plots and len(summary):
        make_plots(summary, out.parent / "summary_v12_plots")

    print(f"Wrote {len(summary)} summarized runs to {out}")
    if errors:
        print("\nErrors:")
        for hist, exc in errors:
            print(f"  {hist}: {exc}")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
