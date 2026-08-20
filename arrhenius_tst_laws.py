#!/usr/bin/env python3
"""Transition-state/Arrhenius kinetic laws for DDD mechanism audits.

This module is intentionally independent of ParaDiS, ExaDiS, OpenDiS, and the
reduced Taylor drivers.  It provides the common rate and barrier functions that
should be used by every mechanism-specific adapter so that Peierls glide,
forest depinning, junction reactions, cross slip, and source activation are all
expressed through the same transition-state theory convention.

Stress convention
-----------------
The barrier input is an effective local stress conjugate to the activated event.
For force-work depinning this should be

    tau_eff = F_PK * x_dagger / v_star

not the geometric diagnostic F_PK/(b L_eff).  This matches the convention used
in the corrected continuous-contact Taylor front work.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Tuple

try:
    import numpy as np
except Exception:  # pragma: no cover - permits documentation imports without numpy
    np = None

KB_EV_K = 8.617333262145e-5
EV_J = 1.602176634e-19


@dataclass(frozen=True)
class ExpFloorBarrier:
    """Exponential-floor activation barrier.

    G(tau, T) = G0(T) [ f + (1-f) exp{-a (tau/sigma_c)^n} ]
    G0(T) = H - k_B T S, where S is supplied in units of k_B.

    sigma_c is in Pa.  The input tau is also in Pa.
    """

    H_eV: float
    S_kB: float
    sigma_c_Pa: float
    f: float
    a: float
    n: float

    def validate(self) -> None:
        if self.sigma_c_Pa <= 0:
            raise ValueError("sigma_c_Pa must be positive")
        if not (0.0 <= self.f <= 1.0):
            raise ValueError("floor fraction f must be in [0, 1]")
        if self.a <= 0.0 or self.n <= 0.0:
            raise ValueError("EXP-floor shape parameters a and n must be positive")

    def G0_eV(self, T_K: float) -> float:
        if T_K <= 0:
            raise ValueError("T_K must be positive")
        return self.H_eV - KB_EV_K * T_K * self.S_kB

    def floor_eV(self, T_K: float) -> float:
        return self.f * self.G0_eV(T_K)

    def barrier_eV(self, tau_eff_Pa, T_K: float):
        self.validate()
        if np is None:
            tau = max(float(tau_eff_Pa), 0.0)
            x = tau / self.sigma_c_Pa
            shape = self.f + (1.0 - self.f) * math.exp(-self.a * (x ** self.n))
            return max(self.G0_eV(T_K) * shape, self.floor_eV(T_K), 0.0)

        tau = np.asarray(tau_eff_Pa, dtype=float)
        x = np.maximum(tau, 0.0) / self.sigma_c_Pa
        shape = self.f + (1.0 - self.f) * np.exp(-self.a * np.power(x, self.n))
        return np.maximum(self.G0_eV(T_K) * shape, max(self.floor_eV(T_K), 0.0))

    def inverse_tau_eff_Pa(self, G_req_eV: float, T_K: float) -> Tuple[float, str]:
        """Invert the EXP-floor branch.

        Returns (tau_eff_Pa, regime).  Regime is one of:
        - zero-stress-transparent: required barrier is at or above G0(T)
        - floor-limited: required barrier is below the residual floor
        - finite: finite stress solution exists
        """
        G0 = self.G0_eV(T_K)
        floor = self.floor_eV(T_K)
        if G_req_eV >= G0:
            return 0.0, "zero-stress-transparent"
        if G_req_eV <= floor:
            return math.inf, "floor-limited"
        y = (G_req_eV / G0 - self.f) / (1.0 - self.f)
        y = min(max(y, 1e-300), 1.0 - 1e-15)
        x = (-math.log(y) / self.a) ** (1.0 / self.n)
        return self.sigma_c_Pa * x, "finite"


@dataclass(frozen=True)
class ArrheniusHazard:
    """Attempt frequency and temperature wrapper for an activation barrier."""

    barrier: ExpFloorBarrier
    eta0_s: float = 1.0e12

    def rate_s(self, tau_eff_Pa, T_K: float):
        G = self.barrier.barrier_eV(tau_eff_Pa, T_K)
        if np is None:
            return self.eta0_s * math.exp(-float(G) / (KB_EV_K * T_K))
        return self.eta0_s * np.exp(-G / (KB_EV_K * T_K))

    def probability(self, tau_eff_Pa, T_K: float, dt_s: float):
        rdt = self.rate_s(tau_eff_Pa, T_K) * dt_s
        if np is None:
            return -math.expm1(-min(max(float(rdt), 0.0), 50.0))
        return -np.expm1(-np.minimum(np.maximum(rdt, 0.0), 50.0))


def force_work_tau_eff_Pa(F_PK_N, x_dagger_m: float, v_star_m3: float):
    """Convert force-work bias to an effective stress."""
    if v_star_m3 <= 0.0:
        raise ValueError("v_star_m3 must be positive")
    return F_PK_N * x_dagger_m / v_star_m3


def signed_forward_minus_reverse_rate_s(
    tau_eff_Pa: float,
    law: ArrheniusHazard,
    T_K: float,
) -> float:
    """Signed Peierls-style forward-minus-reverse rate."""
    return float(law.rate_s(max(tau_eff_Pa, 0.0), T_K) - law.rate_s(max(-tau_eff_Pa, 0.0), T_K))


def required_barrier_for_taylor_flow_eV(
    rho_m2: float,
    T_K: float,
    strain_rate_s: float,
    eta0_s: float,
    b_m: float,
    prefactor: float = 16.0,
) -> float:
    """Analytical residual barrier required by the ideal Arrhenius-Taylor closure.

    This is the same barrier level used to compare an ideal independent-site
    Taylor branch with explicit-contact DDD runs.
    """
    arg = eta0_s * prefactor * rho_m2**2 * b_m**4 / strain_rate_s
    if arg <= 0.0:
        raise ValueError("Taylor flow argument must be positive")
    return KB_EV_K * T_K * math.log(arg)


def ideal_taylor_amplification(rho_m2: float, b_m: float) -> float:
    """Ideal tau_eff/tau_app for L_eff = 1/(2 sqrt(rho))."""
    return 1.0 / (2.0 * b_m * math.sqrt(rho_m2))


def analytical_taylor_exp_floor_stress_MPa(
    rho_m2: float,
    T_K: float,
    strain_rate_s: float,
    b_m: float,
    law: ArrheniusHazard,
) -> tuple[float, str, float, float]:
    """Predict ideal applied Taylor stress for an EXP-floor barrier.

    Returns (tau_app_MPa, regime, G_req_eV, tau_eff_GPa).
    """
    G_req = required_barrier_for_taylor_flow_eV(
        rho_m2=rho_m2,
        T_K=T_K,
        strain_rate_s=strain_rate_s,
        eta0_s=law.eta0_s,
        b_m=b_m,
    )
    tau_eff_Pa, regime = law.barrier.inverse_tau_eff_Pa(G_req, T_K)
    if math.isinf(tau_eff_Pa):
        return math.inf, regime, G_req, math.inf
    M = ideal_taylor_amplification(rho_m2, b_m)
    return tau_eff_Pa / M / 1.0e6, regime, G_req, tau_eff_Pa / 1.0e9


def scan_peak_by_temperature(
    temperatures_K: Iterable[float],
    densities_m2: Iterable[float],
    strain_rate_s: float,
    b_m: float,
    law: ArrheniusHazard,
) -> list[dict]:
    """Return analytical peak summaries for a density grid."""
    out = []
    densities = list(densities_m2)
    for T in temperatures_K:
        vals = []
        regimes = []
        for rho in densities:
            tau, regime, G_req, tau_eff = analytical_taylor_exp_floor_stress_MPa(
                rho_m2=rho,
                T_K=T,
                strain_rate_s=strain_rate_s,
                b_m=b_m,
                law=law,
            )
            vals.append(tau)
            regimes.append(regime)
        finite = [v if math.isfinite(v) else math.nan for v in vals]
        imax = int(np.nanargmax(finite)) if np is not None else max(range(len(finite)), key=lambda i: finite[i])
        out.append({
            "T_K": T,
            "peak_rho_m2": densities[imax],
            "peak_tau_MPa": finite[imax],
            "interior_peak": 0 < imax < len(densities) - 1,
            "monotone_up": all(finite[i] <= finite[i + 1] for i in range(len(finite) - 1)),
            "regimes": ",".join(regimes),
        })
    return out
