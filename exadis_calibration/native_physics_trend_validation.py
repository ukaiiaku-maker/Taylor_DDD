#!/usr/bin/env python3
"""Validate native Arrhenius ExaDiS temperature/rate trends and audit coverage."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


def _case_dir(path: Path) -> Path:
    enabled = path / "audit_enabled"
    return enabled if enabled.is_dir() else path


def _load_case(spec: str) -> dict[str, Any]:
    parts = spec.split(",", 2)
    if len(parts) != 3:
        raise ValueError("case must be T_K,strain_rate_s,path")
    temperature = float(parts[0])
    strain_rate = float(parts[1])
    path = Path(parts[2]).resolve()
    case = _case_dir(path)
    summary = json.loads((case / "final_summary.json").read_text())
    curve = np.atleast_2d(np.loadtxt(case / "stress_strain_dens.dat", comments="#"))
    final_point = np.array([[summary["istep"], summary["strain"], summary["stress_Pa"]]])
    strain_stress = curve[:, [1, 2]]
    if not np.isclose(strain_stress[-1, 0], final_point[0, 1]):
        strain_stress = np.vstack((strain_stress, final_point[:, 1:3]))
    strain_stress = strain_stress[np.argsort(strain_stress[:, 0])]

    events: Counter[str] = Counter()
    accepted: Counter[str] = Counter()
    finite = True
    audited_decisions = True
    large_hazard_without_scheme = 0
    audit_path = case / "event_audit.jsonl"
    with audit_path.open() as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("record_type") != "event":
                continue
            mechanism = str(row.get("mechanism", "unknown"))
            if mechanism not in {
                "mobility_fcc0_arrhenius", "topology_split", "cross_slip"
            }:
                continue
            event_class = str(row.get("event_class", "unknown"))
            events[f"{mechanism}:{event_class}"] += 1
            if row.get("accepted_arrhenius") == 1:
                accepted[f"{mechanism}:{event_class}"] += 1
            for key in ("R_s", "Rdt", "P", "G_used_eV"):
                value = float(row.get(key, 0.0))
                finite = finite and math.isfinite(value)
            if mechanism in {"topology_split", "cross_slip"} and row.get(
                "geometry_admissible"
            ) == 1:
                audited_decisions = audited_decisions and row.get(
                    "accepted_arrhenius"
                ) in (0, 1)
                if float(row.get("Rdt", 0.0)) > 1.0:
                    has_scheme = (
                        float(row.get("accumulated_hazard_after", 0.0)) >= 0.0
                        and float(row.get("selection_threshold", -1.0)) >= 0.0
                    )
                    if not has_scheme:
                        large_hazard_without_scheme += 1

    return {
        "temperature_K": temperature,
        "strain_rate_s": strain_rate,
        "path": str(path),
        "stress_Pa": float(summary["stress_Pa"]),
        "final_strain": float(summary["strain"]),
        "_strain_stress": strain_stress,
        "density_m2": float(summary["density_m2"]),
        "Nnodes": int(summary["Nnodes"]),
        "Nsegs": int(summary["Nsegs"]),
        "network_sane": bool(summary.get("network_sane", True)),
        "audit_finite": finite,
        "arrhenius_decisions_audited": audited_decisions,
        "large_hazard_without_scheme": large_hazard_without_scheme,
        "event_rows": dict(sorted(events.items())),
        "accepted_events": dict(sorted(accepted.items())),
    }


def _nondecreasing(values: list[float], relative_slack: float) -> bool:
    return all(
        right >= left * (1.0 - relative_slack)
        for left, right in zip(values, values[1:])
    )


def _nonincreasing(values: list[float], relative_slack: float) -> bool:
    return all(
        right <= left * (1.0 + relative_slack)
        for left, right in zip(values, values[1:])
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--relative-slack", type=float, default=0.05)
    parser.add_argument("--athermal-plateau-slack", type=float, default=0.03)
    args = parser.parse_args()

    cases = [_load_case(spec) for spec in args.case]
    comparison_strain = min(row["final_strain"] for row in cases)
    for row in cases:
        curve = row.pop("_strain_stress")
        row["comparison_strain"] = comparison_strain
        row["comparison_stress_Pa"] = float(
            np.interp(comparison_strain, curve[:, 0], curve[:, 1])
        )
    temperatures = sorted({row["temperature_K"] for row in cases})
    rates = sorted({row["strain_rate_s"] for row in cases})

    rate_checks: dict[str, bool] = {}
    for temperature in temperatures:
        subset = sorted(
            (row for row in cases if row["temperature_K"] == temperature),
            key=lambda row: row["strain_rate_s"],
        )
        if len(subset) > 1:
            rate_checks[str(temperature)] = _nondecreasing(
                [row["comparison_stress_Pa"] for row in subset], args.relative_slack
            )

    temperature_checks: dict[str, bool] = {}
    temperature_regimes: dict[str, str] = {}
    for rate in rates:
        subset = sorted(
            (row for row in cases if row["strain_rate_s"] == rate),
            key=lambda row: row["temperature_K"],
        )
        if len(subset) > 1:
            stresses = [row["comparison_stress_Pa"] for row in subset]
            softening = _nonincreasing(stresses, args.relative_slack)
            scale = max(1.0, max(abs(value) for value in stresses))
            spread = (max(stresses) - min(stresses)) / scale
            athermal_plateau = spread <= args.athermal_plateau_slack
            temperature_checks[str(rate)] = softening or athermal_plateau
            temperature_regimes[str(rate)] = (
                "thermal_softening" if softening else
                "athermal_plateau" if athermal_plateau else
                "thermal_strengthening_failure"
            )

    hard_gates = {
        "all_networks_sane": all(row["network_sane"] for row in cases),
        "all_hazards_finite": all(row["audit_finite"] for row in cases),
        "all_arrhenius_decisions_audited": all(
            row["arrhenius_decisions_audited"] for row in cases
        ),
        "large_Rdt_has_cumulative_or_high_hazard_scheme": all(
            row["large_hazard_without_scheme"] == 0 for row in cases
        ),
    }
    physics_gates = {
        "higher_rate_strengthens": bool(rate_checks) and all(rate_checks.values()),
        "higher_temperature_softens_or_is_athermal": bool(temperature_checks)
        and all(temperature_checks.values()),
    }
    passed = all(hard_gates.values()) and all(physics_gates.values())
    report = {
        "status": "passed" if passed else "failed",
        "cases": cases,
        "rate_checks_by_temperature_K": rate_checks,
        "temperature_checks_by_strain_rate_s": temperature_checks,
        "temperature_regime_by_strain_rate_s": temperature_regimes,
        "hard_gates": hard_gates,
        "physics_gates": physics_gates,
        "relative_slack": args.relative_slack,
        "athermal_plateau_slack": args.athermal_plateau_slack,
        "comparison_strain": comparison_strain,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
