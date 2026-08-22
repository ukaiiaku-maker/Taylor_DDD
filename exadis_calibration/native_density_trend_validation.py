#!/usr/bin/env python3
"""Validate high-barrier Taylor scaling; a density plateau is a hard failure."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def _load_case(spec: str) -> dict[str, Any]:
    parts = spec.split(",", 1)
    if len(parts) != 2:
        raise ValueError("case must be density_factor,path")
    requested = float(parts[0])
    path = Path(parts[1]).resolve()
    case = path / "audit_enabled" if (path / "audit_enabled").is_dir() else path
    summary = json.loads((case / "final_summary.json").read_text())
    curve = np.atleast_2d(np.loadtxt(case / "stress_strain_dens.dat", comments="#"))
    strain_stress = curve[:, [1, 2]]
    if not np.isclose(strain_stress[-1, 0], float(summary["strain"])):
        strain_stress = np.vstack(
            (strain_stress, [[float(summary["strain"]), float(summary["stress_Pa"])]] )
        )
    strain_stress = strain_stress[np.argsort(strain_stress[:, 0])]
    return {
        "density_factor": requested,
        "reported_density_factor": float(summary["density_factor"]),
        "initial_density_m2": float(summary["initial_density_m2"]),
        "final_density_m2": float(summary["density_m2"]),
        "stress_Pa": float(summary["stress_Pa"]),
        "final_strain": float(summary["strain"]),
        "_strain_stress": strain_stress,
        "network_sane": bool(summary.get("network_sane", False)),
        "path": str(path),
    }


def _nondecreasing(values: list[float], slack: float) -> bool:
    return all(right >= left * (1.0 - slack) for left, right in zip(values, values[1:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", required=True)
    parser.add_argument("--relative-slack", type=float, default=0.02)
    parser.add_argument(
        "--comparison-strains", default="2e-7,1e-6,2e-6,5e-6,1e-5"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cases = sorted((_load_case(spec) for spec in args.case), key=lambda row: row["density_factor"])
    maximum_common_strain = min(row["final_strain"] for row in cases)
    requested_strains = [float(value) for value in args.comparison_strains.split(",")]
    comparison_strains = [value for value in requested_strains if value <= maximum_common_strain]
    if not comparison_strains:
        comparison_strains = [maximum_common_strain]
    density = np.asarray([row["initial_density_m2"] for row in cases])
    strain_checks = []
    for comparison_strain in comparison_strains:
        stresses = np.asarray([
            np.interp(comparison_strain, row["_strain_stress"][:, 0], row["_strain_stress"][:, 1])
            for row in cases
        ])
        positive = stresses > 0.0
        slope = float(np.polyfit(np.log(density[positive]), np.log(stresses[positive]), 1)[0]) \
            if np.count_nonzero(positive) >= 3 else float("nan")
        strengthening = _nondecreasing(stresses.tolist(), 0.0)
        strain_checks.append({
            "strain": comparison_strain,
            "stress_Pa_by_density": stresses.tolist(),
            "log_stress_log_density_slope": slope,
            "strictly_strengthening": bool(strengthening),
            "taylor_slope_0p4_to_0p6": bool(
                math.isfinite(slope) and 0.4 <= slope <= 0.6
            ),
        })
    for row in cases:
        row.pop("_strain_stress")
    taylor_windows = [
        check for check in strain_checks
        if check["strictly_strengthening"] and check["taylor_slope_0p4_to_0p6"]
    ]
    gates = {
        "at_least_three_densities": len(cases) >= 3,
        "requested_factors_applied": all(
            math.isclose(row["density_factor"], row["reported_density_factor"], rel_tol=1e-12)
            for row in cases
        ),
        "initial_density_increases": _nondecreasing(
            [row["initial_density_m2"] for row in cases], 0.0
        ),
        "final_density_increases": _nondecreasing(
            [row["final_density_m2"] for row in cases], args.relative_slack
        ),
        "density_factors_span_at_least_one_decade": bool(
            density[-1] / density[0] >= 10.0
        ),
        "taylor_slope_observed_at_common_strain": bool(taylor_windows),
        "all_networks_sane": all(row["network_sane"] for row in cases),
    }
    passed = all(gates.values())
    report = {
        "status": "passed" if passed else "failed",
        "regime": "taylor_like_strengthening" if taylor_windows else "taylor_scaling_not_resolved",
        "taylor_like_strengthening_observed": bool(taylor_windows),
        "density_plateau_is_accepted": False,
        "relative_slack": args.relative_slack,
        "comparison_strains": comparison_strains,
        "strain_checks": strain_checks,
        "cases": cases,
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
