#!/usr/bin/env python3
"""Fit EXP-floor candidate hazards for native topology and cross-slip rows."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

from anisotropic_hazard_fit import EV_J, ExpFloorParams, exp_floor_rate_s, hazard_probability


def _note_number(note: str, key: str) -> float:
    match = re.search(rf"(?:^|;){re.escape(key)}=([-+0-9.eE]+)", note)
    return float(match.group(1)) if match else math.nan


def _topology_drive(row: dict, burgers_m: float) -> float:
    speeds = [
        np.linalg.norm(row.get(name, [0.0, 0.0, 0.0]))
        for name in ("before_velocity_m_s", "after_velocity0_m_s", "after_velocity1_m_s")
    ]
    speed = max(max(speeds), 1e-30)
    reference_length_m = 300.0 * burgers_m
    return max(float(row.get("delta_power_W", 0.0)), 0.0) / (
        speed * burgers_m * reference_length_m
    )


def _cross_slip_drive(row: dict) -> float:
    note = str(row.get("note", ""))
    primary = _note_number(note, "primary_force_internal")
    cross = _note_number(note, "cross_force_internal")
    threshold = _note_number(note, "stock_force_threshold_internal")
    denom = math.nan
    if primary and np.isfinite(primary):
        denom = primary / float(row.get("tau_local_Pa", 0.0))
    difference = cross - primary
    if (not np.isfinite(denom) or denom == 0.0) and difference and np.isfinite(difference):
        denom = difference / float(row.get("tau_eff_Pa", 0.0))
    if not np.isfinite(denom) or denom == 0.0:
        return 0.0
    return max((abs(cross) - abs(primary) - threshold) / abs(denom), 0.0)


def _auc(score: np.ndarray, label: np.ndarray) -> float:
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    ranks = np.empty(len(score), dtype=float)
    start = 0
    while start < len(score):
        end = start + 1
        while end < len(score) and sorted_score[end] == sorted_score[start]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + 1 + end)
        start = end
    n1 = int(label.sum())
    n0 = len(label) - n1
    return float((ranks[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def load_rows(manifest: dict, mechanism: str, burgers_m: float) -> dict[str, np.ndarray]:
    drives, labels, dts, conditions, initial_states = [], [], [], [], []
    token = f'"mechanism":"{mechanism}"'
    for run in manifest["runs"]:
        with Path(run["audit_jsonl"]).open() as handle:
            for line in handle:
                if token not in line:
                    continue
                row = json.loads(line)
                drive = (
                    _topology_drive(row, burgers_m)
                    if mechanism == "topology_split" else _cross_slip_drive(row)
                )
                drives.append(drive)
                labels.append(int(row["accepted_stock"]) == 1)
                dts.append(max(float(row.get("dt_s", 1e-9)), 1e-30))
                conditions.append(f"{run['run_id']}:step={row['step']}")
                initial_states.append(str(run.get("initial_state_id", "state_1")))
    return {
        "drive": np.asarray(drives, float),
        "label": np.asarray(labels, bool),
        "dt": np.asarray(dts, float),
        "condition": np.asarray(conditions, object),
        "initial_state": np.asarray(initial_states, object),
    }


def fit_mechanism(data: dict[str, np.ndarray], temperatures: list[float], name: str) -> dict:
    drive0, label0, dt0, condition0, initial_state0 = (
        data["drive"], data["label"], data["dt"], data["condition"],
        data["initial_state"],
    )
    drive = np.tile(drive0, len(temperatures))
    label = np.tile(label0, len(temperatures))
    dt = np.tile(dt0, len(temperatures))
    temperature = np.repeat(np.asarray(temperatures, float), len(drive0))
    condition = np.concatenate([
        np.asarray([f"{value}:T={T}" for value in condition0], object)
        for T in temperatures
    ])
    initial_state = np.tile(initial_state0, len(temperatures))
    unique_initial_states = sorted(np.unique(initial_state).tolist())
    if len(unique_initial_states) > 1:
        held = initial_state == unique_initial_states[-1]
    else:
        held = np.asarray([
            sum(str(value).encode("utf-8")) % 5 == 0 for value in condition0
        ])
        held = np.tile(held, len(temperatures))
    if held.sum() == 0 or (~held).sum() == 0:
        held = np.arange(len(label)) % 5 == 0
    train = ~held

    bounds = [
        (0.005, 1.5), (-4.0, 10.0), (-3.0, 1.3),
        (0.002, 0.95), (0.05, 60.0), (0.35, 4.0),
    ]

    def unpack(x: np.ndarray) -> ExpFloorParams:
        return ExpFloorParams(
            H_eV=float(x[0]), S_kB=float(x[1]), sigma_c_GPa=10.0 ** float(x[2]),
            f=float(x[3]), a=float(x[4]), n=float(x[5]), eta0_s=1.0e12,
        )

    positive_weight = max((~label[train]).sum() / max(label[train].sum(), 1), 1.0)

    def objective(x: np.ndarray) -> float:
        p = unpack(x)
        rate = exp_floor_rate_s(drive[train], temperature[train], p)
        probability = np.clip(hazard_probability(rate, dt[train]), 1e-12, 1.0 - 1e-12)
        weights = np.where(label[train], positive_weight, 1.0)
        return float(-np.mean(weights * (
            label[train] * np.log(probability) +
            (~label[train]) * np.log(1.0 - probability)
        )))

    global_fit = differential_evolution(
        objective, bounds, seed=2718, maxiter=24, popsize=7, polish=False
    )
    local_fit = minimize(
        objective, global_fit.x, method="L-BFGS-B", bounds=bounds,
        options={"maxiter": 2000, "ftol": 1e-12},
    )
    xbest = local_fit.x
    params = unpack(xbest)
    rate = exp_floor_rate_s(drive, temperature, params)
    rdt = rate * dt
    probability = hazard_probability(rate, dt)
    # Calibrate the deterministic quantile on training to preserve event count.
    threshold = float(np.quantile(probability[train], 1.0 - label[train].mean()))
    predicted = probability >= threshold

    parameter_names = ["H_eV", "S_kB", "log10_sigma_c_GPa", "f", "a", "n"]
    at_bounds = []
    for parameter, value, (lower, upper) in zip(parameter_names, xbest, bounds):
        if min(value - lower, upper - value) <= 0.005 * (upper - lower):
            at_bounds.append(parameter)

    held_accuracy = float(np.mean(predicted[held] == label[held]))
    observed_fraction = float(label[held].mean())
    predicted_fraction = float(predicted[held].mean())
    balance_error = abs(predicted_fraction - observed_fraction) / max(observed_fraction, 1e-12)
    gates = {
        "heldout_auc_gt_0p8": _auc(probability[held], label[held]) > 0.8,
        "heldout_accuracy_gt_0p8": held_accuracy > 0.8,
        "heldout_class_balance_within_20pct": balance_error <= 0.2,
        "no_parameter_on_bound": not at_bounds,
        "multiple_temperatures": len(temperatures) >= 3,
        "multiple_initial_network_states": len(unique_initial_states) >= 2,
        "Rdt_p95_lt_1": float(np.percentile(rdt, 95)) < 1.0,
    }
    if name == "topology":
        gates["event_classes_fully_resolved"] = False
    if name == "cross_slip":
        gates["rejected_candidates_have_executable_geometry"] = False
    return {
        "status": "fit",
        "mechanism": name,
        "source_rows": len(label0),
        "expanded_rows": len(label),
        "positive_rows": int(label0.sum()),
        "eta0_fixed_s": 1.0e12,
        "exp_floor_params": vars(params),
        "deterministic_quantile": threshold,
        "heldout_auc": _auc(probability[held], label[held]),
        "heldout_accuracy": held_accuracy,
        "heldout_observed_fraction": observed_fraction,
        "heldout_predicted_fraction": predicted_fraction,
        "heldout_balance_relative_error": balance_error,
        "Rdt_median": float(np.median(rdt)),
        "Rdt_p95": float(np.percentile(rdt, 95)),
        "parameters_on_bounds": at_bounds,
        "gates": gates,
        "replacement_eligible": all(gates.values()),
        "replacement_blockers": [key for key, passed in gates.items() if not passed],
        "drive_definition": (
            "positive trial delta power divided by candidate speed, b, and minseg reference length"
            if name == "topology" else
            "positive abs(tau_cross)-abs(tau_primary)-force_threshold stress margin"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("results/exadis_native_discrete_hazard_fit"))
    parser.add_argument("--burgers-m", type=float, default=2.55e-10)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    temperatures = [float(value) for value in manifest["temperatures_K"]]
    fits = []
    for mechanism, label in (("topology_split", "topology"), ("cross_slip", "cross_slip")):
        fits.append(fit_mechanism(
            load_rows(manifest, mechanism, args.burgers_m), temperatures, label
        ))
    summary = {
        "manifest": str(args.manifest),
        "temperatures_K": temperatures,
        "fits": fits,
        "native_replacement_authorized": all(fit["replacement_eligible"] for fit in fits),
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "native_discrete_hazard_fit_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["native_replacement_authorized"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
