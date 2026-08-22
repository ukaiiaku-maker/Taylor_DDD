#!/usr/bin/env python3
"""Fit native EXP-floor mechanisms across a stock ExaDiS campaign manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from anisotropic_hazard_fit import (
    _mobility_node_targets,
    fit_mobility_equivalence,
)


def read_mobility_rows(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            if '"mechanism":"mobility_fcc0"' in line:
                rows.append(json.loads(line))
    return pd.json_normalize(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("results/exadis_native_arrhenius_campaign_fit"))
    parser.add_argument("--burgers-m", type=float, default=2.55e-10)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    temperatures = [float(value) for value in manifest["temperatures_K"]]
    frames = []
    source_rows = 0
    for run in manifest["runs"]:
        path = Path(run["audit_jsonl"])
        frame = read_mobility_rows(path)
        source_rows += len(frame)
        frame["run_id"] = str(run["run_id"])
        frame["initial_state_id"] = str(run.get("initial_state_id", "state_1"))
        frame["strain_rate_s"] = float(run["strain_rate_s"])
        frame["state_id"] = frame["run_id"] + ":step=" + frame["step"].astype(str)
        arm_rows = len(frame)
        frame = _mobility_node_targets(frame, 900.0, args.burgers_m)
        frame["source_arm_rows"] = arm_rows / max(len(frame), 1)
        # Stock FCC_0 has no temperature input.  Its measured response is the
        # calibration target at every requested temperature; this expansion is
        # explicit rather than pretending that duplicate stock simulations are
        # independent temperature-sensitive observations.
        for temperature in temperatures:
            expanded = frame.copy()
            expanded["temperature_K"] = temperature
            expanded["condition_id"] = (
                "rate=" + expanded["strain_rate_s"].astype(str) +
                ":T=" + str(temperature) + ":state=" + expanded["state_id"]
            )
            frames.append(expanded)

    campaign = pd.concat(frames, ignore_index=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    mobility = fit_mobility_equivalence(
        campaign, temperatures[len(temperatures) // 2], args.burgers_m, args.outdir,
        adaptive_event_integration=bool(manifest.get("adaptive_event_integration", False)),
    )
    summary = {
        "status": "passed" if mobility.get("replacement_eligible", False) else "failed",
        "manifest": str(args.manifest),
        "source_rows": source_rows,
        "expanded_rows": len(campaign),
        "temperatures_K": temperatures,
        "stock_temperature_semantics": (
            "FCC_0 has no temperature input; each measured rate/state target is evaluated "
            "against the fitted law at every campaign temperature"
        ),
        "initial_state_ids": sorted(campaign["initial_state_id"].unique().tolist()),
        "eta0_default_s": 1.0e12,
        "mobility": mobility,
        "native_replacement_authorized": bool(mobility.get("replacement_eligible", False)),
    }
    output = args.outdir / "native_arrhenius_campaign_fit_summary.json"
    output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if summary["native_replacement_authorized"]:
        fitted = mobility["exp_floor_params"]
        config = {
            "temperature_K": 900.0,
            "eta0_default_s": 1.0e12,
            "mobility_peierls": {
                "replacement_eligible": True,
                "H_eV": fitted["H_eV"],
                "S_kB": fitted["S_kB"],
                "sigma_c_GPa": fitted["sigma_c_GPa"],
                "f": fitted["f"],
                "a": fitted["a"],
                "n": fitted["n"],
                "jump_b": mobility["jump_b"],
                "vstar_b3": mobility["vstar_characteristic_b3"],
                "anisotropic_coupling": mobility["anisotropic_coupling"],
            },
            "topology": {"replacement_eligible": False},
            "cross_slip": {"replacement_eligible": False},
            "collision": {"replacement_eligible": False},
        }
        (args.outdir / "exadis_arrhenius_gate_passed.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n"
        )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["native_replacement_authorized"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
