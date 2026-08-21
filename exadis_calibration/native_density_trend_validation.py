#!/usr/bin/env python3
"""Validate density scaling and a low-temperature Taylor-like strength trend."""

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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cases = sorted((_load_case(spec) for spec in args.case), key=lambda row: row["density_factor"])
    comparison_strain = min(row["final_strain"] for row in cases)
    for row in cases:
        curve = row.pop("_strain_stress")
        row["comparison_strain"] = comparison_strain
        row["comparison_stress_Pa"] = float(
            np.interp(comparison_strain, curve[:, 0], curve[:, 1])
        )
    stresses = [row["comparison_stress_Pa"] for row in cases]
    strictly_strengthening = _nondecreasing(stresses, 0.0)
    stress_scale = max(1.0, max(abs(value) for value in stresses))
    density_plateau = (max(stresses) - min(stresses)) / stress_scale <= args.relative_slack
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
        "density_response_is_taylor_like_or_plateau": (
            strictly_strengthening or density_plateau
        ),
        "all_networks_sane": all(row["network_sane"] for row in cases),
    }
    passed = all(gates.values())
    report = {
        "status": "passed" if passed else "failed",
        "regime": (
            "taylor_like_strengthening" if strictly_strengthening
            else "athermal_density_plateau_without_resolved_taylor_strengthening"
        ),
        "taylor_like_strengthening_observed": strictly_strengthening,
        "density_plateau_observed": density_plateau,
        "relative_slack": args.relative_slack,
        "comparison_strain": comparison_strain,
        "cases": cases,
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
