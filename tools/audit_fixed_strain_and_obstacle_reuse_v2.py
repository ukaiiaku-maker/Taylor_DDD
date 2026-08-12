#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import pandas as pd
import numpy as np
import sys

def parse_case(name):
    mT = re.search(r"T([0-9.]+)", name)
    mr = re.search(r"rho([0-9.eE+-]+)", name)
    return float(mT.group(1)) if mT else np.nan, float(mr.group(1)) if mr else np.nan

def find_first_existing(d, names):
    for n in names:
        p = d / n
        if p.exists():
            return p
    # fallback
    cands = sorted(d.glob("*history*.csv"))
    if cands:
        return cands[0]
    return None

def read_csv_safe(p):
    if p is None or not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(p)
    except Exception:
        return pd.DataFrame()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--target-strain", type=float, required=True)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--max-depin-per-line-obstacle", type=int, default=1)
    args = ap.parse_args()

    root = Path(args.root)
    failures = []
    rows = []

    case_dirs = [d for d in sorted(root.glob("T*")) if d.is_dir()]
    if not case_dirs:
        raise SystemExit(f"No case directories found under {root}")

    for case_dir in case_dirs:
        hist_path = find_first_existing(case_dir, [
            "single_glider_history.csv",
            "history.csv",
            "run_history.csv",
        ])
        ev_path = find_first_existing(case_dir, [
            "single_glider_crossing_events.csv",
            "crossing_events.csv",
            "events.csv",
        ])

        h = read_csv_safe(hist_path)
        ev = read_csv_safe(ev_path)

        T, rho = parse_case(case_dir.name)

        eps_final = np.nan
        if not h.empty:
            for c in ["eps_total", "total_strain", "strain_total", "eps"]:
                if c in h.columns:
                    eps_final = float(h[c].iloc[-1])
                    break

        strain_ok = np.isfinite(eps_final) and (
            abs(eps_final - args.target_strain) <= max(args.tol, 1e-6*args.target_strain)
        )

        direct_cross = 0
        repeated_pairs = 0
        n_depin = 0
        unique_pairs = 0
        unique_obs = 0
        n_capture = 0

        if not ev.empty and "event" in ev.columns:
            estr = ev["event"].astype(str)
            direct_cross = int((estr == "direct_cross").sum())
            n_capture = int(estr.str.contains("capture", case=False, na=False).sum())
            dep = ev[estr.str.contains("depin|cross", case=False, na=False)].copy()

            # Exclude direct_cross from depin count if present.
            if len(dep):
                dep = dep[dep["event"].astype(str) != "direct_cross"]

            n_depin = len(dep)

            if len(dep) and {"line_id", "obs_id"}.issubset(dep.columns):
                unique_pairs = dep[["line_id", "obs_id"]].drop_duplicates().shape[0]
                counts = dep.groupby(["line_id", "obs_id"]).size()
                repeated_pairs = int((counts > args.max_depin_per_line_obstacle).sum())
                unique_obs = dep["obs_id"].nunique()

        if not strain_ok:
            failures.append((case_dir.name, "final strain missing or not target", eps_final))
        if direct_cross:
            failures.append((case_dir.name, "direct_cross events present", direct_cross))
        if repeated_pairs:
            failures.append((case_dir.name, "repeated line/obstacle depin pairs", repeated_pairs))

        rows.append({
            "case": case_dir.name,
            "T_K": T,
            "rho_m2": rho,
            "history_file": str(hist_path) if hist_path else "MISSING",
            "event_file": str(ev_path) if ev_path else "MISSING",
            "history_rows": len(h),
            "event_rows": len(ev),
            "eps_total_final": eps_final,
            "target_strain": args.target_strain,
            "direct_cross_events": direct_cross,
            "capture_events": n_capture,
            "depin_cross_events": n_depin,
            "unique_line_obstacle_depin_pairs": unique_pairs,
            "unique_depin_obstacles": unique_obs,
            "repeated_line_obstacle_pairs": repeated_pairs,
        })

    out = pd.DataFrame(rows)
    if {"T_K", "rho_m2"}.issubset(out.columns):
        out = out.sort_values(["T_K", "rho_m2"])

    print(out.to_string(index=False))

    (root / "analysis").mkdir(exist_ok=True)
    out.to_csv(root / "analysis" / "fixed_strain_obstacle_reuse_audit_v2.csv", index=False)
    print("\nwrote:", root / "analysis" / "fixed_strain_obstacle_reuse_audit_v2.csv")

    if failures:
        print("\nAUDIT FAILED")
        for f in failures[:100]:
            print(f)
        sys.exit(2)

    print("\nAUDIT PASSED")

if __name__ == "__main__":
    main()
