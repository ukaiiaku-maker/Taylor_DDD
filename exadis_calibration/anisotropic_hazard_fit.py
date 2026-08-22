#!/usr/bin/env python3
"""Fit equivalent anisotropic Arrhenius hazards to instrumented ExaDiS audits.

This module does not replace ExaDiS physics.  It calibrates *equivalent* TST
hazard laws against the stock ExaDiS response recorded by
`exadis_audit.binding_event_audit` so that the fitted laws can later be promoted
into native mobility, cross-slip, and collision hooks.

Supported targets
-----------------
1. Mobility/FCC_0 equivalence:
   Fit a signed forward-minus-reverse Arrhenius glide law to the stock mobility
   velocity projected on each segment endpoint.  This is the rate-law target for
   item 2 in the instrumentation hierarchy.

2. Cross-slip acceptance equivalence:
   Fit an EXP-floor competing-event hazard to stock cross-slip candidate rows
   when the audit exposes accepted/rejected candidate labels.  This is the item
   5 target.

3. Collision acceptance equivalence:
   Fit an EXP-floor activated-collision hazard to rows explicitly labelled as
   activated or ambiguous collision candidates.  Deterministic core-overlap and
   numerical cleanup rows are excluded by default.  This is the item 6 target.

The fitter accepts partial audit output.  Missing mechanisms are reported as
`status=no_data` rather than silently producing parameters.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

try:
    from scipy.optimize import differential_evolution, minimize, brentq
except Exception:  # pragma: no cover
    differential_evolution = None
    minimize = None
    brentq = None

KB_EV_K = 8.617333262145e-5
EV_J = 1.602176634e-19


@dataclass
class ExpFloorParams:
    H_eV: float = 0.50
    S_kB: float = -9.0
    sigma_c_GPa: float = 14.5
    f: float = 0.20
    a: float = 6.65607
    n: float = 2.15276
    eta0_s: float = 1.0e12


@dataclass
class LinearPeierlsParams:
    H_eV: float = 0.05
    S_kB: float = 0.0
    vstar_b3: float = 10.0
    eta0_s: float = 1.0e12
    jump_b: float = 1.0


@dataclass
class AnisotropicCoupling:
    """Minimal anisotropic activation-work coupling.

    tau_eff = tau_s + a_nn sigma_nn + a_mm sigma_mm + a_np sigma_np

    The audit can provide these components directly.  If non-glide components
    are unavailable, the fitter reduces to Schmid-only coupling.
    """

    a_nn: float = 0.0
    a_mm: float = 0.0
    a_np: float = 0.0
    abs_effective_stress: bool = False


def _as_float_array(s: pd.Series, default: float = 0.0) -> np.ndarray:
    return pd.to_numeric(s, errors="coerce").fillna(default).to_numpy(dtype=float)


def _find_first_column(df: pd.DataFrame, names: Iterable[str]) -> Optional[str]:
    for n in names:
        if n in df.columns:
            return n
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def _find_informative_column(df: pd.DataFrame, names: Iterable[str]) -> Optional[str]:
    """Find a numeric column with at least one finite, nonzero observation.

    Native audit rows share a stable schema, so fields that do not apply to a
    mechanism are present as zeros.  Selecting solely by column existence would
    otherwise choose (for example) mobility ``tau_eff_Pa=0`` ahead of its
    informative ``tau_local_Pa`` diagnostic.
    """
    fallback = None
    for name in names:
        column = _find_first_column(df, [name])
        if column is None:
            continue
        if fallback is None:
            fallback = column
        values = pd.to_numeric(df[column], errors="coerce").to_numpy(dtype=float)
        if np.any(np.isfinite(values) & (np.abs(values) > 0.0)):
            return column
    return fallback


def _vector_column(df: pd.DataFrame, name: str, width: int) -> Optional[np.ndarray]:
    column = _find_first_column(df, [name])
    if column is None:
        return None
    out = np.full((len(df), width), np.nan, dtype=float)
    for index, value in enumerate(df[column]):
        if isinstance(value, (list, tuple, np.ndarray)) and len(value) == width:
            try:
                out[index] = np.asarray(value, dtype=float)
            except (TypeError, ValueError):
                pass
    return out


def tensor_components_from_audit(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Resolve Schmid/non-Schmid components from native ``sigma`` and geometry.

    ExaDiS records stress in Voigt order ``xx, yy, zz, yz, xz, xy``.  The slip
    direction is the normalized Burgers vector, ``n`` is the candidate/current
    plane normal, and ``p = n x m`` is the in-plane transverse direction.
    """
    sigma6 = _vector_column(df, "applied_stress_Pa", 6)
    burg = _vector_column(df, "burg", 3)
    plane = _vector_column(df, "plane", 3)
    if sigma6 is None or burg is None or plane is None:
        return {}

    mnorm = np.linalg.norm(burg, axis=1)
    nnorm = np.linalg.norm(plane, axis=1)
    valid = np.isfinite(sigma6).all(axis=1) & (mnorm > 0.0) & (nnorm > 0.0)
    m = np.zeros_like(burg)
    n = np.zeros_like(plane)
    m[valid] = burg[valid] / mnorm[valid, None]
    n[valid] = plane[valid] / nnorm[valid, None]
    p = np.cross(n, m)
    pnorm = np.linalg.norm(p, axis=1)
    valid &= pnorm > 0.0
    p[valid] /= pnorm[valid, None]

    sigma = np.zeros((len(df), 3, 3), dtype=float)
    sigma[:, 0, 0] = sigma6[:, 0]
    sigma[:, 1, 1] = sigma6[:, 1]
    sigma[:, 2, 2] = sigma6[:, 2]
    sigma[:, 1, 2] = sigma[:, 2, 1] = sigma6[:, 3]
    sigma[:, 0, 2] = sigma[:, 2, 0] = sigma6[:, 4]
    sigma[:, 0, 1] = sigma[:, 1, 0] = sigma6[:, 5]

    def contract(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        values = np.einsum("ni,nij,nj->n", left, sigma, right)
        values[~valid] = np.nan
        return values

    return {
        "tau_s_Pa": contract(m, n),
        "sigma_nn_Pa": contract(n, n),
        "sigma_mm_Pa": contract(m, m),
        "sigma_np_Pa": contract(n, p),
    }


def read_jsonl(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    if not rows:
        return pd.DataFrame()
    return pd.json_normalize(rows)


def exp_floor_free_energy_eV(tau_eff_Pa: np.ndarray, T_K: float, p: ExpFloorParams) -> np.ndarray:
    """EXP-floor free energy with entropy outside the enthalpy floor."""
    tau = np.maximum(np.asarray(tau_eff_Pa, dtype=float), 0.0)
    x = tau / (p.sigma_c_GPa * 1.0e9)
    Hshape = p.H_eV * (p.f + (1.0 - p.f) * np.exp(-p.a * np.power(x, p.n)))
    G = Hshape - KB_EV_K * T_K * p.S_kB
    Gfloor = p.f * p.H_eV - KB_EV_K * T_K * p.S_kB
    return np.maximum(G, np.maximum(Gfloor, 0.0))


def exp_floor_rate_s(tau_eff_Pa: np.ndarray, T_K: float, p: ExpFloorParams) -> np.ndarray:
    G = exp_floor_free_energy_eV(tau_eff_Pa, T_K, p)
    return p.eta0_s * np.exp(-G / (KB_EV_K * T_K))


def exp_floor_signed_velocity_m_s(
    tau_signed_Pa: np.ndarray,
    temperature_K: np.ndarray | float,
    p: ExpFloorParams,
    jump_b: float,
    b_m: float,
) -> np.ndarray:
    """Directional forward-minus-reverse EXP-floor Peierls velocity."""
    tau = np.asarray(tau_signed_Pa, dtype=float)
    temperature = np.asarray(temperature_K, dtype=float)
    gp = exp_floor_free_energy_eV_array(np.maximum(tau, 0.0), temperature, p)
    gm = exp_floor_free_energy_eV_array(np.maximum(-tau, 0.0), temperature, p)
    rp = p.eta0_s * np.exp(-gp / (KB_EV_K * temperature))
    rm = p.eta0_s * np.exp(-gm / (KB_EV_K * temperature))
    return jump_b * b_m * (rp - rm)


def exp_floor_free_energy_eV_array(
    tau_eff_Pa: np.ndarray, temperature_K: np.ndarray | float, p: ExpFloorParams
) -> np.ndarray:
    tau = np.maximum(np.asarray(tau_eff_Pa, dtype=float), 0.0)
    temperature = np.asarray(temperature_K, dtype=float)
    x = tau / (p.sigma_c_GPa * 1.0e9)
    enthalpy = p.H_eV * (
        p.f + (1.0 - p.f) * np.exp(-p.a * np.power(x, p.n))
    )
    g = enthalpy - KB_EV_K * temperature * p.S_kB
    floor = p.f * p.H_eV - KB_EV_K * temperature * p.S_kB
    return np.maximum(g, np.maximum(floor, 0.0))


def linear_work_barrier_eV(tau_eff_Pa: np.ndarray, T_K: float, p: LinearPeierlsParams, b_m: float) -> np.ndarray:
    tau = np.maximum(np.asarray(tau_eff_Pa, dtype=float), 0.0)
    vstar = p.vstar_b3 * b_m**3
    enthalpy = np.maximum(0.0, p.H_eV - tau * vstar / EV_J)
    return np.maximum(0.0, enthalpy - KB_EV_K * T_K * p.S_kB)


def linear_signed_velocity_m_s(tau_signed_Pa: np.ndarray, T_K: float, p: LinearPeierlsParams, b_m: float) -> np.ndarray:
    tau = np.asarray(tau_signed_Pa, dtype=float)
    gp = linear_work_barrier_eV(np.maximum(tau, 0.0), T_K, p, b_m)
    gm = linear_work_barrier_eV(np.maximum(-tau, 0.0), T_K, p, b_m)
    rp = p.eta0_s * np.exp(-gp / (KB_EV_K * T_K))
    rm = p.eta0_s * np.exp(-gm / (KB_EV_K * T_K))
    return p.jump_b * b_m * (rp - rm)


def hazard_probability(rate_s: np.ndarray, dt_s: np.ndarray) -> np.ndarray:
    return -np.expm1(-np.clip(rate_s * dt_s, 0.0, 50.0))


def effective_stress_from_audit(df: pd.DataFrame, coupling: AnisotropicCoupling) -> np.ndarray:
    # Preferred direct columns from binding/native audit.
    tau_col = _find_informative_column(df, [
        "tau_eff_Pa",
        "tau_resolved_from_applied_Pa",
        "tau_external_pk_Pa",
        "tau_resolved_external_Pa",
        "tau_local_Pa",
        "tau_from_total_nodal_force_Pa",
    ])
    derived = tensor_components_from_audit(df)
    if tau_col is not None:
        tau = _as_float_array(df[tau_col])
    elif "tau_s_Pa" in derived:
        tau = np.nan_to_num(derived["tau_s_Pa"], nan=0.0)
    else:
        return np.zeros(len(df), dtype=float)
    nn_col = _find_first_column(df, ["sigma_nn_Pa", "non_glide_sigma_nn_Pa"])
    mm_col = _find_first_column(df, ["sigma_mm_Pa", "non_glide_sigma_mm_Pa"])
    np_col = _find_first_column(df, ["sigma_np_Pa", "tau_non_planar_Pa", "secondary_shear_Pa"])

    if nn_col is not None:
        tau = tau + coupling.a_nn * _as_float_array(df[nn_col])
    elif "sigma_nn_Pa" in derived:
        tau = tau + coupling.a_nn * np.nan_to_num(derived["sigma_nn_Pa"], nan=0.0)
    if mm_col is not None:
        tau = tau + coupling.a_mm * _as_float_array(df[mm_col])
    elif "sigma_mm_Pa" in derived:
        tau = tau + coupling.a_mm * np.nan_to_num(derived["sigma_mm_Pa"], nan=0.0)
    if np_col is not None:
        tau = tau + coupling.a_np * _as_float_array(df[np_col])
    elif "sigma_np_Pa" in derived:
        tau = tau + coupling.a_np * np.nan_to_num(derived["sigma_np_Pa"], nan=0.0)

    if coupling.abs_effective_stress:
        tau = np.abs(tau)
    return tau


def _mechanism_mask(df: pd.DataFrame, tokens: Iterable[str]) -> np.ndarray:
    if "mechanism" not in df.columns:
        return np.zeros(len(df), dtype=bool)
    s = df["mechanism"].astype(str).str.lower()
    mask = np.zeros(len(df), dtype=bool)
    for tok in tokens:
        mask |= s.str.contains(tok.lower(), regex=False).to_numpy()
    return mask


def _truthy_mask(series: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    text = series.astype(str).str.lower().isin(["true", "yes", "accepted", "selected"]).to_numpy()
    return (numeric == 1.0) | text


def _mobility_node_targets(m: pd.DataFrame, fallback_T_K: float, b_m: float) -> pd.DataFrame:
    """Collapse arm rows to the actual event-conjugate nodal degree of freedom.

    FCC_0 solves one constrained nodal velocity from one projected nodal force.
    Treating each arm projection as an independent signed event produced the old
    91.5% sign ceiling.  The nodal generalized force is reconstructed from
    power/velocity and divided by half the attached line length, exactly matching
    the stress definition used by the native Arrhenius mobility.
    """
    velocity = _vector_column(m, "velocity_m_s", 3)
    if velocity is None or "node_id" not in m or "step" not in m:
        return pd.DataFrame()
    work = m.copy()
    work["_vx"] = velocity[:, 0]
    work["_vy"] = velocity[:, 1]
    work["_vz"] = velocity[:, 2]
    group_columns = ["step", "node_id"]
    for optional in (
        "condition_id", "state_id", "initial_state_id", "strain_rate_s",
        "temperature_K",
    ):
        if optional in work.columns:
            group_columns.append(optional)

    work["_L_m"] = pd.to_numeric(work["L_m"], errors="coerce")
    work["_dt_s"] = pd.to_numeric(
        work["dt_s"] if "dt_s" in work else 1e-9, errors="coerce"
    )
    if "after_power_W" in work:
        work["_power_W"] = pd.to_numeric(work["after_power_W"], errors="coerce")
    else:
        force = _vector_column(work, "force_N", 3)
        if force is None:
            return pd.DataFrame()
        work["_power_W"] = (
            force[:, 0] * work["_vx"].to_numpy(float)
            + force[:, 1] * work["_vy"].to_numpy(float)
            + force[:, 2] * work["_vz"].to_numpy(float)
        )

    # This aggregation is algebraically identical to the former per-group
    # loop, but keeps multi-state native campaigns out of Python object/swap
    # overhead.  Velocity, power, and dt are nodal values repeated on arm rows;
    # attached arm lengths are the only quantity summed.
    grouped = work.groupby(group_columns, sort=False, dropna=False, as_index=False).agg(
        _vx=("_vx", "first"),
        _vy=("_vy", "first"),
        _vz=("_vz", "first"),
        _power_W=("_power_W", "first"),
        _dt_s=("_dt_s", "first"),
        _total_length_m=("_L_m", "sum"),
    )
    speed = np.sqrt(grouped["_vx"] ** 2 + grouped["_vy"] ** 2 + grouped["_vz"] ** 2)
    half_length = 0.5 * grouped["_total_length_m"].to_numpy(float)
    power = grouped["_power_W"].to_numpy(float)
    dt = grouped["_dt_s"].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = np.abs(power) / speed.to_numpy(float) / (b_m * half_length)
    valid = (
        np.isfinite(speed.to_numpy(float)) & (speed.to_numpy(float) > 0.0)
        & np.isfinite(half_length) & (half_length > 0.0)
        & np.isfinite(power) & (power != 0.0)
        & np.isfinite(tau) & (tau > 0.0)
    )
    nodes = grouped.loc[valid, group_columns].copy()
    sign = np.sign(power[valid])
    nodes["temperature_K"] = (
        pd.to_numeric(nodes["temperature_K"], errors="coerce")
        if "temperature_K" in nodes else fallback_T_K
    )
    nodes["dt_s"] = np.where(dt[valid] > 0.0, dt[valid], 1e-9)
    nodes["tau_event_Pa"] = sign * tau[valid]
    nodes["velocity_event_m_s"] = sign * speed.to_numpy(float)[valid]
    return nodes


def fit_mobility_equivalence(
    df: pd.DataFrame, T_K: float, b_m: float, outdir: Path,
    adaptive_event_integration: bool = False,
) -> dict:
    if {"tau_event_Pa", "velocity_event_m_s"}.issubset(df.columns):
        m = df
        nodes = df.copy()
        arm_row_count = int(df.get("source_arm_rows", pd.Series([len(df)])).sum())
    else:
        mask = _mechanism_mask(df, ["mobility_fcc0"])
        m = df[mask].copy()
        if m.empty:
            return {"status": "no_data", "mechanism": "mobility"}
        nodes = _mobility_node_targets(m, T_K, b_m)
        arm_row_count = int(len(m))
    if len(nodes) < 20:
        return {"status": "too_few_nodal_rows", "mechanism": "mobility", "rows": int(len(nodes))}

    tau = nodes["tau_event_Pa"].to_numpy(float)
    velocity = nodes["velocity_event_m_s"].to_numpy(float)
    temperature = nodes["temperature_K"].to_numpy(float)
    dt = nodes["dt_s"].to_numpy(float)
    condition_count = int(nodes.get("condition_id", pd.Series(["single"])).nunique())
    temperature_count = int(np.unique(temperature).size)
    state_count = int(nodes.get("state_id", pd.Series(["single"])).nunique())
    initial_state_count = int(
        nodes.get("initial_state_id", pd.Series(["single"])).nunique()
    )

    if "initial_state_id" in nodes and initial_state_count > 1:
        # Keep an entire originating network out of the optimizer.  All of its
        # nodes, rates, evolved snapshots, and temperature evaluations remain
        # held out, avoiding row-level or duplicate-temperature leakage.
        held_state = sorted(nodes["initial_state_id"].astype(str).unique())[-1]
        held = (nodes["initial_state_id"].astype(str) == held_state).to_numpy()
    elif "state_id" in nodes and state_count > 1:
        held = nodes["state_id"].astype(str).map(
            lambda value: sum(value.encode("utf-8")) % 5 == 0
        ).to_numpy()
    else:
        held = ((nodes["node_id"].to_numpy(int) * 37 + nodes["step"].to_numpy(int)) % 5) == 0
    if held.sum() < 5 or (~held).sum() < 10:
        held = np.arange(len(nodes)) % 5 == 0
    train = ~held
    optimize = train.copy()
    train_indices = np.flatnonzero(train)
    if len(train_indices) > 50_000:
        rng = np.random.default_rng(1701)
        keep = rng.choice(train_indices, size=50_000, replace=False)
        optimize[:] = False
        optimize[keep] = True

    bounds = [
        (0.005, 1.5),      # H_eV
        (-4.0, 10.0),      # S_kB
        (-3.0, 1.3),       # log10 sigma_c_GPa
        (0.002, 0.95),     # floor fraction
        (0.05, 60.0),      # shape a
        (0.35, 4.0),       # shape n
        (-2.0, math.log10(50.0)),  # log10 jump_b
    ]

    def unpack(x: np.ndarray) -> tuple[ExpFloorParams, float]:
        return ExpFloorParams(
            H_eV=float(x[0]), S_kB=float(x[1]), sigma_c_GPa=10.0 ** float(x[2]),
            f=float(x[3]), a=float(x[4]), n=float(x[5]), eta0_s=1.0e12,
        ), 10.0 ** float(x[6])

    def objective(x: np.ndarray, selected: np.ndarray = optimize) -> float:
        p, jump = unpack(x)
        pred = exp_floor_signed_velocity_m_s(tau[selected], temperature[selected], p, jump, b_m)
        error = np.log10(np.abs(pred) + 1e-30) - np.log10(np.abs(velocity[selected]) + 1e-30)
        sign_penalty = 4.0 * (np.sign(pred) != np.sign(velocity[selected]))
        return float(np.mean(error * error + sign_penalty))

    x0 = np.array([0.08, 1.0, -0.3, 0.08, 4.0, 1.0, 0.5])
    if differential_evolution is not None:
        result_global = differential_evolution(
            objective, bounds=bounds, maxiter=24, popsize=7, polish=False,
            seed=1701, workers=1,
        )
        start = result_global.x
        method = "differential_evolution+L-BFGS-B"
    else:
        start = x0
        method = "initial_only_no_scipy"
    if minimize is not None:
        result = minimize(objective, start, bounds=bounds, method="L-BFGS-B",
                          options={"maxiter": 3000, "ftol": 1e-12})
        xbest = result.x
    else:
        xbest = start
    pbest, jump_best = unpack(xbest)
    implied_vstar_m3 = (
        pbest.H_eV * (1.0 - pbest.f) * pbest.a * pbest.n /
        (pbest.sigma_c_GPa * 1.0e9) * EV_J
    )
    implied_vstar_b3 = implied_vstar_m3 / (b_m ** 3)
    prediction = exp_floor_signed_velocity_m_s(tau, temperature, pbest, jump_best, b_m)

    def metrics(selected: np.ndarray) -> dict:
        observed = velocity[selected]
        predicted = prediction[selected]
        log_error = np.log10(np.abs(predicted) + 1e-30) - np.log10(np.abs(observed) + 1e-30)
        ratio = np.abs(predicted) / np.maximum(np.abs(observed), 1e-30)
        return {
            "rows": int(selected.sum()),
            "rmse_log10_abs_velocity": float(np.sqrt(np.mean(log_error**2))),
            "sign_accuracy": float(np.mean(np.sign(predicted) == np.sign(observed))),
            "median_velocity_ratio": float(np.median(ratio)),
        }

    train_metrics = metrics(train)
    held_metrics = metrics(held)
    gp = exp_floor_free_energy_eV_array(np.maximum(tau, 0.0), temperature, pbest)
    gm = exp_floor_free_energy_eV_array(np.maximum(-tau, 0.0), temperature, pbest)
    rate_max = np.maximum(
        pbest.eta0_s * np.exp(-gp / (KB_EV_K * temperature)),
        pbest.eta0_s * np.exp(-gm / (KB_EV_K * temperature)),
    )
    rdt = rate_max * dt
    at_bounds = []
    parameter_names = ["H_eV", "S_kB", "log10_sigma_c_GPa", "f", "a", "n", "log10_jump_b"]
    for name, value, (lower, upper) in zip(parameter_names, xbest, bounds):
        span = upper - lower
        if min(value - lower, upper - value) <= 0.005 * span:
            at_bounds.append(name)

    multi_condition = temperature_count >= 3 and state_count >= 2 and condition_count >= 8
    gates = {
        "training_rmse_lt_0p5": train_metrics["rmse_log10_abs_velocity"] < 0.5,
        "heldout_rmse_lt_0p75": held_metrics["rmse_log10_abs_velocity"] < 0.75,
        "sign_accuracy_gt_0p97": held_metrics["sign_accuracy"] > 0.97,
        "median_ratio_0p5_to_2": 0.5 <= held_metrics["median_velocity_ratio"] <= 2.0,
        "rdt_stable_or_adaptive": (
            float(np.percentile(rdt, 95.0)) < 0.2 or adaptive_event_integration
        ),
        "no_parameter_on_bound": not at_bounds,
        "multi_temperature_state_campaign": multi_condition,
        "multiple_initial_network_states": initial_state_count >= 2,
    }
    replacement_eligible = all(gates.values())
    blockers = [name for name, passed in gates.items() if not passed]

    table = nodes.copy()
    table["fit_velocity_m_s"] = prediction
    table["held_out"] = held
    table["Rdt_max_direction"] = rdt
    table.to_csv(outdir / "mobility_fit_observed_vs_predicted.csv", index=False)
    return {
        "status": "fit",
        "replacement_eligible": replacement_eligible,
        "replacement_blockers": blockers,
        "mechanism": "mobility",
        "arm_rows": arm_row_count,
        "nodal_rows": int(len(nodes)),
        "optimization_rows": int(optimize.sum()),
        "method": method,
        "temperature_count": temperature_count,
        "state_count": state_count,
        "initial_state_count": initial_state_count,
        "condition_count": condition_count,
        "burgers_m": b_m,
        "eta0_fixed_s": 1.0e12,
        "exp_floor_params": asdict(pbest),
        "jump_b": jump_best,
        "vstar_characteristic_b3": implied_vstar_b3,
        "vstar_definition": "H*(1-f)*a*n/sigma_c characteristic phi*V* scale; local -dG/dtau also carries x^(n-1)*exp(-a*x^n)",
        "anisotropic_coupling": asdict(AnisotropicCoupling()),
        "training": train_metrics,
        "held_out": held_metrics,
        "Rdt_median": float(np.median(rdt)),
        "Rdt_p95": float(np.percentile(rdt, 95.0)),
        "adaptive_event_integration": adaptive_event_integration,
        "parameters_on_bounds": at_bounds,
        "gates": gates,
        "stress_definition": "abs(projected_nodal_power)/speed/(b*half_attached_line_length)",
    }


def fit_binary_hazard(
    df: pd.DataFrame,
    T_K: float,
    outdir: Path,
    mechanism_name: str,
    tokens: Iterable[str],
    exclude_deterministic: bool = True,
) -> dict:
    mask = _mechanism_mask(df, tokens)
    m = df[mask].copy()
    rows_before_exclusion = int(len(m))
    if exclude_deterministic and "deterministic_geometry_only" in m.columns:
        det = _truthy_mask(m["deterministic_geometry_only"])
        m = m[~det].copy()
    if m.empty:
        status = "insufficient_candidate_labels" if rows_before_exclusion else "no_data"
        return {
            "status": status,
            "mechanism": mechanism_name,
            "rows_before_deterministic_exclusion": rows_before_exclusion,
            "reason": "all observed candidates were deterministic geometry/cleanup" if rows_before_exclusion else "no mechanism rows",
            "replacement_eligible": False,
        }

    acc_col = _find_first_column(m, ["accepted_stock", "accepted", "selected", "event_accepted"])
    if acc_col is None:
        return {"status": "no_acceptance_column", "mechanism": mechanism_name, "rows": int(len(m))}
    numeric_labels = pd.to_numeric(m[acc_col], errors="coerce").to_numpy(dtype=float)
    text_labels = m[acc_col].astype(str).str.lower()
    valid = np.isin(numeric_labels, [0.0, 1.0]) | text_labels.isin(
        ["true", "false", "yes", "no", "accepted", "rejected", "selected"]
    ).to_numpy()
    m = m.iloc[np.where(valid)[0]].copy()
    if m.empty:
        return {
            "status": "insufficient_candidate_labels",
            "mechanism": mechanism_name,
            "rows_before_label_filter": rows_before_exclusion,
            "replacement_eligible": False,
        }
    y = _truthy_mask(m[acc_col]).astype(float)
    if np.unique(y).size < 2:
        return {
            "status": "single_acceptance_class",
            "mechanism": mechanism_name,
            "rows": int(len(m)),
            "replacement_eligible": False,
        }

    dt_col = _find_first_column(m, ["dt_s", "dt", "realdt_s"])
    if dt_col is None:
        dt = np.full(len(m), 1e-9)
    else:
        dt = _as_float_array(m[dt_col], default=1e-9)
        dt = np.where(dt > 0.0, dt, 1e-9)

    derived = tensor_components_from_audit(m)
    has_nn = _find_first_column(m, ["sigma_nn_Pa", "non_glide_sigma_nn_Pa"]) is not None or "sigma_nn_Pa" in derived
    has_mm = _find_first_column(m, ["sigma_mm_Pa", "non_glide_sigma_mm_Pa"]) is not None or "sigma_mm_Pa" in derived
    has_np = _find_first_column(m, ["sigma_np_Pa", "tau_non_planar_Pa", "secondary_shear_Pa"]) is not None or "sigma_np_Pa" in derived

    def nll(x):
        # x = [H, log10_sigma_c_GPa, f_logit-ish, a, n, log10_eta0, a_nn, a_mm, a_np]
        f = 1.0 / (1.0 + np.exp(-x[2]))
        p = ExpFloorParams(
            H_eV=max(x[0], 1e-6),
            S_kB=-9.0,
            sigma_c_GPa=10.0 ** x[1],
            f=f,
            a=max(x[3], 1e-6),
            n=max(x[4], 1e-6),
            eta0_s=10.0 ** x[5],
        )
        c = AnisotropicCoupling(
            a_nn=x[6] if has_nn else 0.0,
            a_mm=x[7] if has_mm else 0.0,
            a_np=x[8] if has_np else 0.0,
            abs_effective_stress=True,
        )
        tau = effective_stress_from_audit(m, c)
        rate = exp_floor_rate_s(tau, T_K, p)
        prob = np.clip(hazard_probability(rate, dt), 1e-12, 1.0 - 1e-12)
        return float(-np.nanmean(y * np.log(prob) + (1.0 - y) * np.log(1.0 - prob)))

    x0 = np.array([0.50, math.log10(14.5), math.log(0.20/0.80), 6.65607, 2.15276, 12.0, 0.0, 0.0, 0.0])
    bounds = [(0.01, 3.0), (-2.0, 2.0), (-6.0, 3.0), (0.01, 20.0), (0.25, 5.0), (6.0, 14.5), (-3.0, 3.0), (-3.0, 3.0), (-3.0, 3.0)]
    if differential_evolution is not None and len(m) >= 20:
        res0 = differential_evolution(nll, bounds=bounds, maxiter=50, polish=False, seed=17)
        start = res0.x
        method = "differential_evolution+minimize"
    else:
        start = x0
        method = "minimize_only_or_initial"
    if minimize is not None:
        res = minimize(nll, start, bounds=bounds, method="Nelder-Mead", options={"maxiter": 3000})
        xbest = res.x
    else:
        xbest = start
        method = "initial_only_no_scipy"

    fbest = 1.0 / (1.0 + np.exp(-xbest[2]))
    pbest = ExpFloorParams(
        H_eV=max(xbest[0], 1e-6),
        S_kB=-9.0,
        sigma_c_GPa=10.0 ** xbest[1],
        f=fbest,
        a=max(xbest[3], 1e-6),
        n=max(xbest[4], 1e-6),
        eta0_s=10.0 ** xbest[5],
    )
    cbest = AnisotropicCoupling(
        a_nn=xbest[6] if has_nn else 0.0,
        a_mm=xbest[7] if has_mm else 0.0,
        a_np=xbest[8] if has_np else 0.0,
        abs_effective_stress=True,
    )
    tau = effective_stress_from_audit(m, cbest)
    rate = exp_floor_rate_s(tau, T_K, pbest)
    rdt = rate * dt
    prob = hazard_probability(rate, dt)
    pred = prob >= 0.5
    out = pd.DataFrame({
        "accepted_stock": y,
        "hazard_probability": prob,
        "hazard_rate_s": rate,
        "Rdt": rdt,
        "tau_eff_fit_Pa": tau,
        "dt_s": dt,
    })
    out.to_csv(outdir / f"{mechanism_name}_fit_observed_vs_predicted.csv", index=False)

    return {
        "status": "fit",
        "replacement_eligible": False,
        "replacement_blockers": [
            "single-temperature deterministic stock decisions do not identify H, S, eta0, and activation volume uniquely",
            "site multiplicity and event strain increment are not identified by this audit",
            "held-out multi-temperature trajectory validation has not been run",
        ],
        "mechanism": mechanism_name,
        "rows": int(len(m)),
        "positive_fraction": float(np.mean(y)),
        "method": method,
        "negative_log_likelihood": nll(xbest),
        "temperature_K": T_K,
        "exp_floor_params": asdict(pbest),
        "anisotropic_coupling": asdict(cbest),
        "classification_accuracy_p05": float(np.mean(pred == y.astype(bool))),
        "probability_median": float(np.nanmedian(prob)),
        "Rdt_median": float(np.nanmedian(rdt)),
        "Rdt_p95": float(np.nanpercentile(rdt, 95.0)),
        "tau_eff_median_GPa": float(np.nanmedian(tau) / 1e9),
    }


def calibrate_native_stress_anchor(stress_strain_path: Optional[Path], strain_rate_s: float) -> dict:
    if stress_strain_path is None or not stress_strain_path.exists():
        return {"status": "no_stress_strain_file"}
    try:
        dat = pd.read_csv(stress_strain_path, comment="#", sep=r"\s+", header=None)
    except Exception as exc:
        return {"status": "read_failed", "error": str(exc)}
    if dat.shape[1] < 4 or dat.empty:
        return {"status": "empty_or_bad_format"}
    native_columns = ["step", "strain", "stress", "density", "Nnodes", "Nsegs", "dt", "time"]
    if dat.shape[1] > len(native_columns):
        return {"status": "unsupported_column_count", "columns": int(dat.shape[1])}
    dat.columns = native_columns[:dat.shape[1]]
    tail = dat.tail(max(3, len(dat)//5))
    return {
        "status": "read",
        "strain_rate_s": strain_rate_s,
        "tail_stress_native": float(tail["stress"].median()),
        "tail_density_native": float(tail["density"].median()),
        "final_stress_native": float(dat["stress"].iloc[-1]),
        "final_density_native": float(dat["density"].iloc[-1]),
        "rows": int(len(dat)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("audit_jsonl", type=Path, help="event_audit.jsonl from binding/native ExaDiS audit")
    ap.add_argument("--outdir", type=Path, default=Path("results/exadis_hazard_fit"))
    ap.add_argument("--temperature-K", type=float, default=900.0)
    ap.add_argument("--burgers-m", type=float, default=2.55e-10)
    ap.add_argument("--strain-rate-s", type=float, default=1.0e3)
    ap.add_argument("--stress-strain-dens", type=Path, default=None)
    ap.add_argument("--no-cross-slip", action="store_true")
    ap.add_argument("--no-collision", action="store_true")
    ap.add_argument("--adaptive-event-integration", action="store_true")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = read_jsonl(args.audit_jsonl)
    summary = {
        "input_audit_jsonl": str(args.audit_jsonl),
        "rows_total": int(len(df)),
        "temperature_K": args.temperature_K,
        "burgers_m": args.burgers_m,
        "strain_rate_s": args.strain_rate_s,
        "arrhenius_replacements_connected": False,
        "native_stress_anchor": calibrate_native_stress_anchor(args.stress_strain_dens, args.strain_rate_s),
        "fits": [],
        "notes": [
            "These are calibrated mechanism surrogates, not universal barrier constants.",
            "No fit is eligible for native replacement until H, S, effective activation volume phi*V*, site multiplicity, and event strain increment are identified and validated out of sample.",
            "Mobility fit uses directional forward-minus-reverse EXP-floor Peierls kinetics at fixed eta0=1e12 s^-1.",
            "Cross-slip and collision fits use binary stock acceptance labels when available.",
            "Independent pathways combine by summing hazards; sequential obstacles require renewal/residence-time treatment and must not be collapsed into a hazard sum.",
        ],
    }

    if df.empty:
        summary["fits"].append({"status": "no_data", "mechanism": "all"})
    else:
        summary["fits"].append(fit_mobility_equivalence(
            df, args.temperature_K, args.burgers_m, args.outdir,
            adaptive_event_integration=args.adaptive_event_integration,
        ))
        if not args.no_cross_slip:
            summary["fits"].append(fit_binary_hazard(df, args.temperature_K, args.outdir, "cross_slip", ["cross_slip", "cross-slip", "xslip"]))
        if not args.no_collision:
            summary["fits"].append(fit_binary_hazard(df, args.temperature_K, args.outdir, "collision", ["collision", "annihilation"]))

    summary["native_replacement_authorized"] = bool(summary["fits"]) and all(
        fit.get("replacement_eligible", False) for fit in summary["fits"]
    )

    with (args.outdir / "anisotropic_hazard_fit_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
