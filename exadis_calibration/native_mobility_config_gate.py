#!/usr/bin/env python3
"""Evaluate an exact native mobility JSON block on a frozen campaign table."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from anisotropic_hazard_fit import (
    ExpFloorParams,
    exp_floor_free_energy_eV_array,
    exp_floor_signed_velocity_m_s,
)

KB_EV_K = 8.617333262145e-5


def _metrics(observed: np.ndarray, predicted: np.ndarray, selected: np.ndarray) -> dict:
    log_error = (
        np.log10(np.abs(predicted[selected]) + 1e-30)
        - np.log10(np.abs(observed[selected]) + 1e-30)
    )
    ratio = np.abs(predicted[selected]) / np.maximum(np.abs(observed[selected]), 1e-30)
    return {
        "rows": int(selected.sum()),
        "rmse_log10_abs_velocity": float(np.sqrt(np.mean(log_error**2))),
        "sign_accuracy": float(
            np.mean(np.sign(predicted[selected]) == np.sign(observed[selected]))
        ),
        "median_velocity_ratio": float(np.median(ratio)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--burgers-m", type=float, default=2.55e-10)
    parser.add_argument("--adaptive-event-integration", action="store_true")
    args = parser.parse_args()

    table = pd.read_csv(args.table)
    config = json.loads(args.config.read_text())
    block = config["mobility_peierls"]
    p = ExpFloorParams(
        H_eV=float(block["H_eV"]),
        S_kB=float(block["S_kB"]),
        sigma_c_GPa=float(block["sigma_c_GPa"]),
        f=float(block["f"]),
        a=float(block["a"]),
        n=float(block["n"]),
        eta0_s=float(config.get("eta0_default_s", 1.0e12)),
    )
    tau = table["tau_event_Pa"].to_numpy(float)
    observed = table["velocity_event_m_s"].to_numpy(float)
    temperature = table["temperature_K"].to_numpy(float)
    dt = table["dt_s"].to_numpy(float)
    held = table["held_out"].astype(str).str.lower().isin(("true", "1")).to_numpy()
    train = ~held
    predicted = exp_floor_signed_velocity_m_s(
        tau, temperature, p, float(block["jump_b"]), args.burgers_m
    )
    gp = exp_floor_free_energy_eV_array(np.maximum(tau, 0.0), temperature, p)
    gm = exp_floor_free_energy_eV_array(np.maximum(-tau, 0.0), temperature, p)
    rate = np.maximum(
        p.eta0_s * np.exp(-gp / (KB_EV_K * temperature)),
        p.eta0_s * np.exp(-gm / (KB_EV_K * temperature)),
    )
    rdt = rate * dt
    training = _metrics(observed, predicted, train)
    held_out = _metrics(observed, predicted, held)
    ratio = held_out["median_velocity_ratio"]
    bounded_values = {
        "H_eV": (float(block["H_eV"]), 0.005, 1.5),
        "S_kB": (float(block["S_kB"]), -4.0, 10.0),
        "log10_sigma_c_GPa": (math.log10(float(block["sigma_c_GPa"])), -3.0, 1.3),
        "f": (float(block["f"]), 0.002, 0.95),
        "a": (float(block["a"]), 0.05, 60.0),
        "n": (float(block["n"]), 0.35, 4.0),
        "log10_jump_b": (math.log10(float(block["jump_b"])), -2.0, math.log10(50.0)),
    }
    on_bounds = []
    for name, (value, lower, upper) in bounded_values.items():
        if min(value - lower, upper - value) <= 0.005 * (upper - lower):
            on_bounds.append(name)
    initial_states = sorted(table["initial_state_id"].astype(str).unique().tolist())
    gates = {
        "training_rmse_lt_0p5": training["rmse_log10_abs_velocity"] < 0.5,
        "heldout_rmse_lt_0p75": held_out["rmse_log10_abs_velocity"] < 0.75,
        "sign_accuracy_gt_0p97": held_out["sign_accuracy"] > 0.97,
        "median_ratio_0p5_to_2": 0.5 <= ratio <= 2.0,
        "Rdt_stable_or_adaptive": (
            float(np.percentile(rdt, 95)) < 0.2 or args.adaptive_event_integration
        ),
        "no_parameter_on_bound": not on_bounds,
        "multiple_initial_network_states": len(initial_states) >= 2,
        "multiple_temperatures": table["temperature_K"].nunique() >= 3,
    }
    result = {
        "status": "passed" if all(gates.values()) else "failed",
        "config": str(args.config),
        "table": str(args.table),
        "initial_state_ids": initial_states,
        "training": training,
        "held_out": held_out,
        "Rdt_median": float(np.median(rdt)),
        "Rdt_p95": float(np.percentile(rdt, 95)),
        "parameters_on_bounds": on_bounds,
        "adaptive_event_integration": args.adaptive_event_integration,
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
