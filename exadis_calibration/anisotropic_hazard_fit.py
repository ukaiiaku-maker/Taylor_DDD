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


def fit_mobility_equivalence(df: pd.DataFrame, T_K: float, b_m: float, outdir: Path) -> dict:
    mask = _mechanism_mask(df, ["mobility", "mobility_force", "fcc0"])
    m = df[mask].copy()
    if m.empty:
        return {"status": "no_data", "mechanism": "mobility"}

    vel_col = _find_first_column(m, [
        "velocity_glide_m_s",
        "glide_velocity_m_s",
        "v_glide_m_s",
        "velocity_m_s",
        "node_velocity_projected_m_s",
        "speed_m_s",
    ])
    if vel_col is None:
        # Fallback: infer speed from vector components if present.
        vx = _find_first_column(m, ["node_velocity_x_m_s", "velocity_x_m_s", "v_x"])
        vy = _find_first_column(m, ["node_velocity_y_m_s", "velocity_y_m_s", "v_y"])
        vz = _find_first_column(m, ["node_velocity_z_m_s", "velocity_z_m_s", "v_z"])
        if vx and vy and vz:
            y = np.linalg.norm(np.column_stack([_as_float_array(m[vx]), _as_float_array(m[vy]), _as_float_array(m[vz])]), axis=1)
        else:
            return {"status": "no_velocity_column", "mechanism": "mobility", "rows": int(len(m))}
    else:
        y = _as_float_array(m[vel_col])

    base_coupling = AnisotropicCoupling()
    tau0 = effective_stress_from_audit(m, base_coupling)
    finite = np.isfinite(y) & np.isfinite(tau0) & (np.abs(y) > 0.0) & (np.abs(tau0) > 0.0)
    y = y[finite]
    tau0 = tau0[finite]
    mfit = m.iloc[np.where(finite)[0]].copy()
    if len(y) < 8:
        return {"status": "too_few_rows", "mechanism": "mobility", "rows": int(len(y))}

    derived = tensor_components_from_audit(mfit)
    has_nn = _find_first_column(mfit, ["sigma_nn_Pa", "non_glide_sigma_nn_Pa"]) is not None or "sigma_nn_Pa" in derived
    has_mm = _find_first_column(mfit, ["sigma_mm_Pa", "non_glide_sigma_mm_Pa"]) is not None or "sigma_mm_Pa" in derived
    has_np = _find_first_column(mfit, ["sigma_np_Pa", "tau_non_planar_Pa", "secondary_shear_Pa"]) is not None or "sigma_np_Pa" in derived

    def objective(x):
        # x = [log10_eta0, log10_vstar_b3, jump_b, a_nn, a_mm, a_np]
        p = LinearPeierlsParams(
            H_eV=0.05,
            S_kB=0.0,
            vstar_b3=10.0 ** x[1],
            eta0_s=10.0 ** x[0],
            jump_b=max(x[2], 1e-6),
        )
        coupling = AnisotropicCoupling(
            a_nn=x[3] if has_nn else 0.0,
            a_mm=x[4] if has_mm else 0.0,
            a_np=x[5] if has_np else 0.0,
            abs_effective_stress=False,
        )
        tau = effective_stress_from_audit(mfit, coupling)
        pred = linear_signed_velocity_m_s(tau, T_K, p, b_m)
        # Fit logarithmic magnitude and sign consistency.
        eps = 1e-30
        err_mag = np.log10(np.abs(pred) + eps) - np.log10(np.abs(y) + eps)
        sign_penalty = 2.0 * (np.sign(pred) != np.sign(y)).astype(float)
        return float(np.nanmean(err_mag**2 + sign_penalty))

    x0 = np.array([12.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    bounds = [(8.0, 14.5), (-1.0, 3.0), (1e-3, 20.0), (-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)]
    if minimize is None:
        xbest, loss = x0, objective(x0)
        method = "initial_only_no_scipy"
    else:
        res = minimize(objective, x0, bounds=bounds, method="Nelder-Mead", options={"maxiter": 2000})
        xbest = res.x
        loss = objective(xbest)
        method = "scipy_minimize"

    pbest = LinearPeierlsParams(
        H_eV=0.05,
        S_kB=0.0,
        vstar_b3=10.0 ** xbest[1],
        eta0_s=10.0 ** xbest[0],
        jump_b=max(xbest[2], 1e-6),
    )
    cbest = AnisotropicCoupling(
        a_nn=xbest[3] if has_nn else 0.0,
        a_mm=xbest[4] if has_mm else 0.0,
        a_np=xbest[5] if has_np else 0.0,
    )
    tau_best = effective_stress_from_audit(mfit, cbest)
    pred_best = linear_signed_velocity_m_s(tau_best, T_K, pbest, b_m)

    table = pd.DataFrame({
        "native_velocity_m_s": y,
        "fit_velocity_m_s": pred_best,
        "tau_eff_fit_Pa": tau_best,
        "tau_eff_base_Pa": tau0,
    })
    table.to_csv(outdir / "mobility_fit_observed_vs_predicted.csv", index=False)

    return {
        "status": "fit",
        "replacement_eligible": False,
        "replacement_blockers": [
            "single-temperature audit cannot identify activation enthalpy and entropy separately",
            "FCC_0 arm projection is a calibrated surrogate of shared nodal mobility, not an event-conjugate barrier",
            "site multiplicity and event strain increment are not identified by this audit",
        ],
        "mechanism": "mobility",
        "rows": int(len(y)),
        "method": method,
        "loss": loss,
        "temperature_K": T_K,
        "burgers_m": b_m,
        "peierls_linear_work_params": asdict(pbest),
        "anisotropic_coupling": asdict(cbest),
        "rmse_log10_abs_velocity": float(np.sqrt(np.nanmean((np.log10(np.abs(pred_best)+1e-30)-np.log10(np.abs(y)+1e-30))**2))),
        "sign_accuracy": float(np.mean(np.sign(pred_best) == np.sign(y))),
        "native_velocity_median_abs_m_s": float(np.nanmedian(np.abs(y))),
        "fit_velocity_median_abs_m_s": float(np.nanmedian(np.abs(pred_best))),
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
            "Mobility fit uses signed linear-work Peierls kinetics to match FCC_0 velocities.",
            "Cross-slip and collision fits use binary stock acceptance labels when available.",
            "Independent pathways combine by summing hazards; sequential obstacles require renewal/residence-time treatment and must not be collapsed into a hazard sum.",
        ],
    }

    if df.empty:
        summary["fits"].append({"status": "no_data", "mechanism": "all"})
    else:
        summary["fits"].append(fit_mobility_equivalence(df, args.temperature_K, args.burgers_m, args.outdir))
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
