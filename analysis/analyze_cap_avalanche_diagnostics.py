#!/usr/bin/env python3
"""
Analyze cap activity and avalanche-like burst structure in OpenDiS v14/v15/v16 runs.

Inputs:
  --run-dir <directory containing single_glider_history.csv>
  --root <results root; recursively scans for single_glider_history.csv>

Outputs per run:
  cap_avalanche_summary.txt/json/csv
  cap_avalanche_events.csv

Outputs for --root:
  cap_avalanche_summary_all.csv

This is a postprocessor only; it does not modify simulation results.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def _to_num(df: pd.DataFrame, cols=None) -> pd.DataFrame:
    if cols is None:
        cols = df.columns
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    x = np.abs(x)
    s = np.sum(x)
    if s <= 0:
        return 0.0
    x = np.sort(x)
    n = x.size
    return float((2.0 * np.sum((np.arange(1, n + 1) * x)) / (n * s)) - (n + 1.0) / n)


def _safe(v, default=float("nan")):
    try:
        if pd.isna(v):
            return default
        return float(v)
    except Exception:
        return default


def _parse_run_metadata(run_dir: Path) -> Dict[str, float]:
    meta = {"rho_m2": float("nan"), "T_K": float("nan")}
    name = run_dir.name
    # common: T1100_rho5e15
    try:
        if "T" in name and "_rho" in name:
            tpart = name.split("_rho")[0]
            rpart = name.split("_rho", 1)[1]
            if tpart.startswith("T"):
                meta["T_K"] = float(tpart[1:].replace("p", "."))
            meta["rho_m2"] = float(rpart.replace("p", "."))
    except Exception:
        pass
    # read run summary if present
    js = run_dir / "run_summary.json"
    if js.exists():
        try:
            d = json.loads(js.read_text())
            for k in ["rho_m2", "T_K", "strain_rate_s", "dt_s"]:
                if k in d:
                    meta[k] = _safe(d[k])
        except Exception:
            pass
    return meta


def _cluster_steps(active_steps: np.ndarray, max_gap: int) -> List[Tuple[int, int]]:
    active_steps = np.asarray(active_steps, dtype=int)
    if active_steps.size == 0:
        return []
    active_steps = np.unique(np.sort(active_steps))
    clusters = []
    start = int(active_steps[0])
    prev = start
    for s in active_steps[1:]:
        s = int(s)
        if s - prev <= max_gap:
            prev = s
        else:
            clusters.append((start, prev))
            start = prev = s
    clusters.append((start, prev))
    return clusters


def analyze_run(run_dir: Path, tau_cap_mpa: float, cap_tol: float, cluster_gap_steps: int,
                active_plastic_ratio: float, window_steps: int) -> Dict[str, float]:
    run_dir = Path(run_dir)
    hist_path = run_dir / "single_glider_history.csv"
    event_path = run_dir / "single_glider_crossing_events.csv"

    if not hist_path.exists():
        raise FileNotFoundError(f"No single_glider_history.csv in {run_dir}")

    hist = pd.read_csv(hist_path, low_memory=False)
    num_cols = [
        "step", "time_s", "eps_total", "eps_plastic", "d_eps_total", "d_eps_p_actual",
        "d_tau_step_MPa", "tau_MPa", "tau_after_step_MPa", "n_depin", "n_capture",
        "crossing_rate_max_s", "frac_tau_local_capped", "tau_local_median_MPa",
        "tau_local_p90_MPa", "tau_local_max_MPa", "n_live_pins",
    ]
    hist = _to_num(hist, num_cols)
    hist = hist.dropna(subset=["step"]).copy()
    hist["step"] = hist["step"].astype(int)

    if "d_eps_total" not in hist.columns:
        # Infer for older histories.
        if "eps_total" in hist.columns:
            hist["d_eps_total"] = hist["eps_total"].diff().fillna(hist["eps_total"])
        else:
            hist["d_eps_total"] = np.nan
    if "d_eps_p_actual" not in hist.columns:
        if "eps_plastic" in hist.columns:
            hist["d_eps_p_actual"] = hist["eps_plastic"].diff().fillna(hist["eps_plastic"])
        else:
            hist["d_eps_p_actual"] = np.nan
    if "d_tau_step_MPa" not in hist.columns:
        tau_col = "tau_after_step_MPa" if "tau_after_step_MPa" in hist.columns else "tau_MPa"
        hist["d_tau_step_MPa"] = hist[tau_col].diff().fillna(0.0)

    depin_col = "n_depin" if "n_depin" in hist.columns else None
    if depin_col is None:
        hist["n_depin"] = 0.0
    else:
        hist["n_depin"] = hist["n_depin"].fillna(0.0)

    tau_col = "tau_after_step_MPa" if "tau_after_step_MPa" in hist.columns else "tau_MPa"
    hist[tau_col] = pd.to_numeric(hist[tau_col], errors="coerce")

    # Add event-level depin info if present.
    ev = None
    if event_path.exists():
        ev = pd.read_csv(event_path, low_memory=False)
        ev = _to_num(ev, ["step", "time_s", "tau_MPa", "tau_local_MPa", "rate_s", "barrier_eV",
                          "work_eV", "F_pin_N", "F_line_tension_N", "cross_force_ratio",
                          "cross_drive_work_eV", "pileup_contributors", "L_feed_reduced"])
        if "event" in ev.columns:
            dep_ev = ev[ev["event"].astype(str).eq("depin_cross")].copy()
        else:
            dep_ev = ev.copy()
    else:
        dep_ev = pd.DataFrame()

    cap_threshold = tau_cap_mpa * (1.0 - cap_tol)

    # History-level cap metrics.
    tail_n = max(10, int(0.2 * len(hist)))
    tail = hist.tail(tail_n)
    hist_cap_tail_med = _safe(tail["frac_tau_local_capped"].median()) if "frac_tau_local_capped" in tail.columns else float("nan")
    hist_tau_loc_tail_med = _safe(tail["tau_local_median_MPa"].median()) if "tau_local_median_MPa" in tail.columns else float("nan")
    hist_live_tail_med = _safe(tail["n_live_pins"].median()) if "n_live_pins" in tail.columns else float("nan")

    # Event-level cap metrics: distinguishes persistent capped survivors from actual capped depin events.
    if len(dep_ev):
        n_dep_event = int(len(dep_ev))
        frac_depin_at_cap = float(np.mean(dep_ev["tau_local_MPa"].to_numpy(dtype=float) >= cap_threshold)) if "tau_local_MPa" in dep_ev.columns else float("nan")
        dep_tau_local_median = _safe(dep_ev["tau_local_MPa"].median()) if "tau_local_MPa" in dep_ev.columns else float("nan")
        dep_tau_local_p90 = _safe(dep_ev["tau_local_MPa"].quantile(0.90)) if "tau_local_MPa" in dep_ev.columns else float("nan")
        dep_tau_local_p99 = _safe(dep_ev["tau_local_MPa"].quantile(0.99)) if "tau_local_MPa" in dep_ev.columns else float("nan")
        dep_rate_dt_max = _safe((dep_ev["rate_s"] * _safe(hist["time_s"].diff().median(), 0.0)).max()) if "rate_s" in dep_ev.columns else float("nan")
        dep_by_step = dep_ev.groupby("step").size().rename("n_depin_event").reset_index()
        hist = hist.merge(dep_by_step, on="step", how="left")
        hist["n_depin_event"] = hist["n_depin_event"].fillna(0.0)
    else:
        n_dep_event = 0
        frac_depin_at_cap = float("nan")
        dep_tau_local_median = float("nan")
        dep_tau_local_p90 = float("nan")
        dep_tau_local_p99 = float("nan")
        dep_rate_dt_max = float("nan")
        hist["n_depin_event"] = hist["n_depin"]

    # Active step definition.
    d_eps_total = hist["d_eps_total"].to_numpy(dtype=float)
    d_eps_p = hist["d_eps_p_actual"].to_numpy(dtype=float)
    ndep = hist["n_depin_event"].to_numpy(dtype=float)
    plastic_ratio = np.divide(d_eps_p, d_eps_total, out=np.zeros_like(d_eps_p), where=np.abs(d_eps_total) > 1e-300)
    active = (ndep > 0) | (plastic_ratio > active_plastic_ratio)
    active_steps = hist.loc[active, "step"].to_numpy(dtype=int)
    clusters = _cluster_steps(active_steps, cluster_gap_steps)

    event_rows = []
    total_depin = float(np.nansum(ndep))
    for a, b in clusters:
        seg = hist[(hist["step"] >= a) & (hist["step"] <= b)].copy()
        if len(seg) == 0:
            continue
        event_size_depin = float(np.nansum(seg["n_depin_event"].to_numpy(dtype=float)))
        plastic_size = float(np.nansum(seg["d_eps_p_actual"].to_numpy(dtype=float)))
        imposed_size = float(np.nansum(seg["d_eps_total"].to_numpy(dtype=float)))
        dtau = seg["d_tau_step_MPa"].to_numpy(dtype=float)
        stress_drop_negsum = float(np.nansum(np.maximum(-dtau, 0.0)))
        tau_vals = seg[tau_col].to_numpy(dtype=float)
        peak_to_valley = float(np.nanmax(tau_vals) - np.nanmin(tau_vals)) if np.isfinite(tau_vals).any() else float("nan")
        event_rows.append({
            "start_step": a,
            "end_step": b,
            "duration_steps": int(b - a + 1),
            "event_size_depin": event_size_depin,
            "plastic_size": plastic_size,
            "imposed_size": imposed_size,
            "plastic_over_imposed": plastic_size / imposed_size if imposed_size > 0 else float("nan"),
            "stress_drop_negsum_MPa": stress_drop_negsum,
            "peak_to_valley_stress_drop_MPa": peak_to_valley,
            "max_n_depin_step": float(np.nanmax(seg["n_depin_event"].to_numpy(dtype=float))),
            "max_plastic_ratio_step": float(np.nanmax(np.divide(
                seg["d_eps_p_actual"].to_numpy(dtype=float),
                seg["d_eps_total"].to_numpy(dtype=float),
                out=np.zeros(len(seg)), where=np.abs(seg["d_eps_total"].to_numpy(dtype=float)) > 1e-300
            ))),
            "mean_frac_tau_local_capped": _safe(seg["frac_tau_local_capped"].mean()) if "frac_tau_local_capped" in seg.columns else float("nan"),
        })

    events_df = pd.DataFrame(event_rows)
    events_df.to_csv(run_dir / "cap_avalanche_events.csv", index=False)

    # Burstiness metrics.
    if total_depin > 0 and len(hist):
        dep_step_counts = ndep[np.isfinite(ndep)]
        largest_event = float(events_df["event_size_depin"].max()) if len(events_df) else 0.0
        largest_fraction = largest_event / total_depin if total_depin > 0 else float("nan")
        sorted_steps = np.sort(dep_step_counts)[::-1]
        top_n = max(1, int(math.ceil(0.01 * len(sorted_steps))))
        top_1pct_fraction = float(np.sum(sorted_steps[:top_n]) / np.sum(sorted_steps)) if np.sum(sorted_steps) > 0 else 0.0
        depin_step_gini = _gini(dep_step_counts)
        # Fano factor using windows of fixed step count.
        steps = hist["step"].to_numpy(dtype=int)
        if len(steps):
            smin, smax = int(np.nanmin(steps)), int(np.nanmax(steps))
            bins = np.arange(smin, smax + window_steps + 1, window_steps)
            # Map per-row counts by step; assumes rows nearly one per step.
            counts, _ = np.histogram(hist["step"].to_numpy(dtype=int), bins=bins, weights=ndep)
            mean_counts = np.mean(counts) if len(counts) else 0.0
            fano = float(np.var(counts) / mean_counts) if mean_counts > 0 else 0.0
        else:
            fano = float("nan")
    else:
        largest_event = 0.0
        largest_fraction = 0.0
        top_1pct_fraction = 0.0
        depin_step_gini = 0.0
        fano = 0.0

    largest_drop = float(events_df["stress_drop_negsum_MPa"].max()) if len(events_df) else 0.0

    # Conservative flag; tune as needed.
    avalanche_like = bool(
        (largest_fraction > 0.03 and largest_drop > 0.5) or
        (fano > 5.0 and top_1pct_fraction > 0.25) or
        (largest_fraction > 0.05)
    )

    meta = _parse_run_metadata(run_dir)
    summary = {
        "run_dir": str(run_dir),
        **meta,
        "n_history_rows": int(len(hist)),
        "total_depin_history_or_event": total_depin,
        "total_depin_event_file": n_dep_event,
        "n_avalanches": int(len(events_df)),
        "largest_event_avalanche_depin": largest_event,
        "largest_event_avalanche_fraction": largest_fraction,
        "largest_stress_drop_negsum_MPa": largest_drop,
        "event_count_window_fano": fano,
        "top_1pct_steps_event_fraction": top_1pct_fraction,
        "depin_step_gini": depin_step_gini,
        "history_tail_frac_tau_local_capped_median": hist_cap_tail_med,
        "history_tail_tau_local_median_MPa": hist_tau_loc_tail_med,
        "history_tail_n_live_pins_median": hist_live_tail_med,
        "frac_depin_events_at_tau_cap": frac_depin_at_cap,
        "depin_tau_local_median_MPa": dep_tau_local_median,
        "depin_tau_local_p90_MPa": dep_tau_local_p90,
        "depin_tau_local_p99_MPa": dep_tau_local_p99,
        "depin_event_rate_dt_max_est": dep_rate_dt_max,
        "avalanche_like": avalanche_like,
    }

    # Interpretation flag.
    if hist_cap_tail_med >= 0.5 and (not math.isfinite(frac_depin_at_cap) or frac_depin_at_cap < 0.1):
        summary["cap_interpretation"] = "persistent_live_pin_cap_not_depin_event_cap"
    elif hist_cap_tail_med >= 0.5 and math.isfinite(frac_depin_at_cap) and frac_depin_at_cap >= 0.5:
        summary["cap_interpretation"] = "depin_events_cap_dominated"
    elif hist_cap_tail_med >= 0.5:
        summary["cap_interpretation"] = "live_pin_cap_active"
    else:
        summary["cap_interpretation"] = "cap_not_dominant_in_tail"

    (run_dir / "cap_avalanche_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    pd.DataFrame([summary]).to_csv(run_dir / "cap_avalanche_summary.csv", index=False)

    lines = ["Cap + avalanche diagnostic", "==========================", ""]
    for k, v in summary.items():
        lines.append(f"{k}: {v}")
    (run_dir / "cap_avalanche_summary.txt").write_text("\n".join(lines) + "\n")
    return summary


def discover_runs(root: Path) -> List[Path]:
    root = Path(root)
    return sorted(p.parent for p in root.rglob("single_glider_history.csv"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--run-dir", type=Path, default=None)
    ap.add_argument("--tau-cap-MPa", type=float, default=8000.0)
    ap.add_argument("--cap-tol", type=float, default=1e-6)
    ap.add_argument("--cluster-gap-steps", type=int, default=1)
    ap.add_argument("--active-plastic-ratio", type=float, default=1.0)
    ap.add_argument("--window-steps", type=int, default=500)
    ap.add_argument("--show-table", action="store_true")
    args = ap.parse_args()

    if args.run_dir is None and args.root is None:
        ap.error("Provide --run-dir or --root")

    if args.run_dir is not None:
        runs = [args.run_dir]
        out_root = args.run_dir
    else:
        runs = discover_runs(args.root)
        out_root = args.root

    summaries = []
    for rd in runs:
        try:
            summaries.append(analyze_run(rd, args.tau_cap_MPa, args.cap_tol,
                                         args.cluster_gap_steps, args.active_plastic_ratio,
                                         args.window_steps))
        except Exception as e:
            summaries.append({"run_dir": str(rd), "error": str(e)})
            print(f"ERROR in {rd}: {e}")

    all_df = pd.DataFrame(summaries)
    if args.root is not None:
        out = Path(out_root) / "cap_avalanche_summary_all.csv"
        all_df.to_csv(out, index=False)
        print(f"Wrote: {out}")

    if args.show_table and len(all_df):
        cols = [
            "rho_m2",
            "total_depin_history_or_event",
            "largest_event_avalanche_depin",
            "largest_event_avalanche_fraction",
            "largest_stress_drop_negsum_MPa",
            "event_count_window_fano",
            "top_1pct_steps_event_fraction",
            "history_tail_frac_tau_local_capped_median",
            "frac_depin_events_at_tau_cap",
            "depin_tau_local_median_MPa",
            "depin_tau_local_p90_MPa",
            "cap_interpretation",
            "avalanche_like",
        ]
        cols = [c for c in cols if c in all_df.columns]
        print(all_df[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
