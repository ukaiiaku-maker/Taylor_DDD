#!/usr/bin/env python3
"""Compare a stock native ExaDiS trajectory with a staged Arrhenius run."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np


def _load_case(path: Path) -> tuple[dict, np.ndarray, Path]:
    case = path / "audit_enabled" if (path / "audit_enabled").is_dir() else path
    summary = json.loads((case / "final_summary.json").read_text())
    curve = np.atleast_2d(np.loadtxt(case / "stress_strain_dens.dat", comments="#"))
    if not np.isclose(
        float(curve[-1, 1]), float(summary["strain"]), rtol=0.0, atol=1.0e-20
    ):
        final_row = np.zeros(curve.shape[1], dtype=float)
        final_row[0] = float(summary["istep"])
        final_row[1] = float(summary["strain"])
        final_row[2] = float(summary["stress_Pa"])
        final_row[3] = float(summary["density_m2"])
        if curve.shape[1] > 4:
            final_row[4] = float(summary["Nnodes"])
        if curve.shape[1] > 5:
            final_row[5] = float(summary["Nsegs"])
        if curve.shape[1] > 6:
            final_row[6] = float(summary["dt_s"])
        if curve.shape[1] > 7:
            final_row[7] = float(summary["time_s"])
        curve = np.vstack((curve, final_row))
    return summary, curve, case / "event_audit.jsonl"


def _relative(left: float, right: float) -> float:
    return abs(right - left) / max(abs(left), np.finfo(float).tiny)


def _curve_metrics(stock: np.ndarray, candidate: np.ndarray) -> dict:
    upper = min(float(stock[:, 1].max()), float(candidate[:, 1].max()))
    stock_common = stock[stock[:, 1] <= upper]
    if len(stock_common) < 2:
        raise ValueError("trajectories do not have two overlapping strain samples")
    strain = stock_common[:, 1]
    stress_candidate = np.interp(strain, candidate[:, 1], candidate[:, 2])
    density_candidate = np.interp(strain, candidate[:, 1], candidate[:, 3])
    stress_scale = max(float(np.max(np.abs(stock_common[:, 2]))), 1.0)
    density_scale = np.maximum(np.abs(stock_common[:, 3]), 1.0)
    return {
        "overlap_max_strain": upper,
        "overlap_rows": len(stock_common),
        "stress_normalized_rmse": float(
            np.sqrt(np.mean((stress_candidate - stock_common[:, 2]) ** 2)) / stress_scale
        ),
        "stress_max_abs_difference_over_stock_peak": float(
            np.max(np.abs(stress_candidate - stock_common[:, 2])) / stress_scale
        ),
        "density_max_relative_difference": float(
            np.max(np.abs(density_candidate - stock_common[:, 3]) / density_scale)
        ),
    }


def _hazard_metrics(path: Path) -> dict:
    values = {key: [] for key in ("Rdt", "P", "G_used_eV")}
    events: Counter[str] = Counter()
    rows = 0
    with path.open() as handle:
        for line in handle:
            if '"mechanism":"mobility_fcc0_arrhenius"' not in line:
                continue
            row = json.loads(line)
            rows += 1
            events[str(row.get("event_class", "unknown"))] += 1
            for key in values:
                values[key].append(float(row[key]))
    if not rows:
        raise ValueError(f"no native Arrhenius mobility rows in {path}")
    arrays = {key: np.asarray(value, float) for key, value in values.items()}
    finite = all(np.isfinite(value).all() for value in arrays.values())
    return {
        "rows": rows,
        "event_classes": dict(sorted(events.items())),
        "all_finite": finite,
        "Rdt_median": float(np.median(arrays["Rdt"])),
        "Rdt_p95": float(np.percentile(arrays["Rdt"], 95)),
        "Rdt_max": float(np.max(arrays["Rdt"])),
        "P_median": float(np.median(arrays["P"])),
        "P_p95": float(np.percentile(arrays["P"], 95)),
        "G_used_eV_min": float(np.min(arrays["G_used_eV"])),
        "G_used_eV_max": float(np.max(arrays["G_used_eV"])),
        "transparent_barrier_fraction": float(np.mean(arrays["G_used_eV"] <= 0.0)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stock", type=Path, required=True)
    parser.add_argument("--arrhenius", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stage", default="A1_arrhenius_peierls_only")
    parser.add_argument("--trajectory-tolerance", type=float, default=0.20)
    parser.add_argument(
        "--adaptive-event-integration", action="store_true",
        help="declare that the native nonlinear subcycling integrator handled multi-hit rates",
    )
    args = parser.parse_args()

    stock_summary, stock_curve, _ = _load_case(args.stock)
    arr_summary, arr_curve, arr_audit = _load_case(args.arrhenius)
    curve = _curve_metrics(stock_curve, arr_curve)
    final = {
        "stress_relative_difference": _relative(
            float(stock_summary["stress_Pa"]), float(arr_summary["stress_Pa"])
        ),
        "density_relative_difference": _relative(
            float(stock_summary["density_m2"]), float(arr_summary["density_m2"])
        ),
        "node_count_relative_difference": _relative(
            float(stock_summary["Nnodes"]), float(arr_summary["Nnodes"])
        ),
        "segment_count_relative_difference": _relative(
            float(stock_summary["Nsegs"]), float(arr_summary["Nsegs"])
        ),
        "plastic_strain_relative_difference": _relative(
            float(stock_summary["pstrain"]), float(arr_summary["pstrain"])
        ),
    }
    hazard = _hazard_metrics(arr_audit)
    tol = args.trajectory_tolerance
    gates = {
        "stress_curve_within_tolerance": (
            curve["stress_max_abs_difference_over_stock_peak"] <= tol
        ),
        "final_stress_within_tolerance": final["stress_relative_difference"] <= tol,
        "density_curve_within_tolerance": curve["density_max_relative_difference"] <= tol,
        "final_density_within_tolerance": final["density_relative_difference"] <= tol,
        "node_count_within_tolerance": final["node_count_relative_difference"] <= tol,
        "segment_count_within_tolerance": final["segment_count_relative_difference"] <= tol,
        "hazard_audit_finite": hazard["all_finite"],
        "no_unintended_transparent_barrier": hazard["transparent_barrier_fraction"] == 0.0,
        "Rdt_stable_or_adaptive": (
            hazard["Rdt_p95"] < 0.2 or args.adaptive_event_integration
        ),
    }
    output = {
        "stage": args.stage,
        "status": "passed" if all(gates.values()) else "failed",
        "trajectory_tolerance": tol,
        "stock": str(args.stock),
        "arrhenius": str(args.arrhenius),
        "final_state": final,
        "curve": curve,
        "hazard": hazard,
        "adaptive_event_integration": args.adaptive_event_integration,
        "stock_event_counts": stock_summary.get("audit_rows_by_mechanism", {}),
        "arrhenius_event_counts": arr_summary.get("audit_rows_by_mechanism", {}),
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0 if output["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
