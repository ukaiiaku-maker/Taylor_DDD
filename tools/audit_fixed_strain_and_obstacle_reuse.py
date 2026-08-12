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

    for case_dir in sorted(root.glob("T*")):
        hist = case_dir / "single_glider_history.csv"
        evp = case_dir / "single_glider_crossing_events.csv"

        if not hist.exists():
            failures.append((case_dir.name, "missing history", ""))
            continue

        h = pd.read_csv(hist)
        if len(h) == 0:
            failures.append((case_dir.name, "empty history", ""))
            continue

        eps_final = float(h["eps_total"].iloc[-1]) if "eps_total" in h.columns else np.nan
        strain_ok = abs(eps_final - args.target_strain) <= max(args.tol, 1e-6*args.target_strain)

        direct_cross = 0
        repeated_pairs = 0
        n_depin = 0
        unique_pairs = 0
        unique_obs = 0

        if evp.exists():
            ev = pd.read_csv(evp)
            if len(ev):
                if "event" in ev.columns:
                    direct_cross = int((ev["event"].astype(str) == "direct_cross").sum())
                    dep = ev[ev["event"].astype(str).str.contains("depin|cross", case=False, na=False)].copy()
                else:
                    dep = pd.DataFrame()

                if len(dep) and {"line_id", "obs_id"}.issubset(dep.columns):
                    n_depin = len(dep)
                    unique_pairs = dep[["line_id", "obs_id"]].drop_duplicates().shape[0]
                    counts = dep.groupby(["line_id", "obs_id"]).size()
                    repeated_pairs = int((counts > args.max_depin_per_line_obstacle).sum())
                    unique_obs = dep["obs_id"].nunique()

        if not strain_ok:
            failures.append((case_dir.name, "final eps_total != target", eps_final))
        if direct_cross:
            failures.append((case_dir.name, "direct_cross events present", direct_cross))
        if repeated_pairs:
            failures.append((case_dir.name, "repeated line/obstacle depin pairs", repeated_pairs))

        T, rho = parse_case(case_dir.name)
        rows.append({
            "case": case_dir.name,
            "T_K": T,
            "rho_m2": rho,
            "eps_total_final": eps_final,
            "target_strain": args.target_strain,
            "direct_cross_events": direct_cross,
            "depin_cross_events": n_depin,
            "unique_line_obstacle_depin_pairs": unique_pairs,
            "unique_depin_obstacles": unique_obs,
            "repeated_line_obstacle_pairs": repeated_pairs,
        })

    out = pd.DataFrame(rows).sort_values(["T_K", "rho_m2"])
    print(out.to_string(index=False))

    (root / "analysis").mkdir(exist_ok=True)
    out.to_csv(root / "analysis" / "fixed_strain_obstacle_reuse_audit.csv", index=False)
    print("\nwrote:", root / "analysis" / "fixed_strain_obstacle_reuse_audit.csv")

    if failures:
        print("\nAUDIT FAILED")
        for f in failures[:50]:
            print(f)
        sys.exit(2)

    print("\nAUDIT PASSED")

if __name__ == "__main__":
    main()
