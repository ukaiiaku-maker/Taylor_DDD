#!/usr/bin/env python3
"""Aggregate required Taylor contact diagnostics from native ExaDiS sweeps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def finite_median(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else None


def load_case(path: Path) -> dict:
    case = path / "audit_enabled" if (path / "audit_enabled").is_dir() else path
    summary = json.loads((case / "final_summary.json").read_text())
    rows = []
    accepted = []
    with (case / "event_audit.jsonl").open() as handle:
        for line in handle:
            row = json.loads(line)
            if not row.get("interaction_class") or row.get("kinetically_eligible") != 1:
                continue
            rows.append(row)
            if row.get("accepted_arrhenius") == 1:
                accepted.append(row)
    volume = float(summary.get("volume_m3", 0.0))
    unique_candidates = {int(r.get("residence_key", 0)) for r in rows}
    load_bearing = [r for r in rows if abs(float(r.get("F_event_used_N", 0.0))) > 0.0]
    transparent = [
        r for r in rows
        if abs(float(r.get("G_used_eV", 0.0)) - float(r.get("G_floor_eV", 0.0))) <= 1.0e-12
    ]
    return {
        "path": str(path),
        "density_factor": float(summary["density_factor"]),
        "total_line_density_m2": float(summary["density_m2"]),
        "mobile_line_density_m2": summary.get("mobile_line_density_m2"),
        "forest_intersecting_line_density_m2": summary.get("forest_intersecting_line_density_m2"),
        "junction_density_m3": summary.get("junction_density_m3"),
        "interaction_candidate_density_m3": len(unique_candidates) / volume if volume > 0.0 else None,
        "accepted_interaction_density_m3": len(accepted) / volume if volume > 0.0 else None,
        "mean_segment_length_m": summary.get("mean_segment_length_m"),
        "median_L_eff_m": finite_median([r.get("L_eff_m") for r in rows]),
        "median_L_eff_over_b": finite_median([r.get("phi_geom_L_over_b") for r in rows]),
        "median_phi_eff": finite_median([
            r.get("phi_eff")
            for r in rows
            if abs(float(r.get("tau_app_resolved_Pa", 0.0))) > 0.0
        ]),
        "median_tau_eff_Taylor_Pa": finite_median([r.get("tau_eff_Pa") for r in rows]),
        "median_tau_eff_over_tau_app": finite_median([
            float(r.get("tau_eff_Pa", 0.0)) / float(r.get("tau_app_resolved_Pa", 1.0))
            for r in rows if abs(float(r.get("tau_app_resolved_Pa", 0.0))) > 0.0
        ]),
        "load_bearing_candidate_density_m3": len({int(r.get("residence_key", 0)) for r in load_bearing}) / volume if volume > 0.0 else None,
        "transparent_candidate_density_m3": len({int(r.get("residence_key", 0)) for r in transparent}) / volume if volume > 0.0 else None,
        "interaction_rows": len(rows),
        "accepted_interaction_rows": len(accepted),
        "network_sane": bool(summary["network_sane"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = sorted((load_case(path) for path in args.cases), key=lambda x: x["density_factor"])
    result = {"status": "passed", "cases": cases}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
