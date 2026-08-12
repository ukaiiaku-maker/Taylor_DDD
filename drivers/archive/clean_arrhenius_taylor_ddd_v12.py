#!/usr/bin/env python3
"""
clean_arrhenius_taylor_ddd_v8.py

Self-contained explicit-obstacle Arrhenius Taylor / Peierls DDD-style test driver.

v10: swept-capture plus center-of-mass-projected backstress mobility. This version is intentionally not a patch of previous scripts.  It keeps the
useful framework/output style, but rebuilds the mechanics consistently.

Core physics:
  - fixed imposed total strain rate and fixed target strain;
  - plastic strain from swept area of connected mobile dislocation lines;
  - explicit forest obstacles; density enters only through obstacle count/positions;
  - no analytic Taylor prefactor, no tau_app*X/b, no density multiplier in hazard;
  - first obstacle contact creates a pinned junction;
  - only pinned lead junctions can depin/cross;
  - queued same-line/same-obstacle nodes contribute pile-up count only;
  - one mobile line crosses a fixed obstacle at most once per imposed-strain pass;
  - Peierls/free glide and forest crossing both use:
        DeltaG(tau,T) = G_fit(tau,T) - kB*T*S_kB
    with no extra -tau*v*(tau) subtraction;
  - forest crossing local stress is obtained from local pin force:
        tau_local * vstar(tau_local,T) = F_pin * b.

Important mechanical change relative to the old prototype:
  free nodes do not move as completely independent particles.  After Peierls
  glide, a line-tension smoothing/relaxation step is applied while pinned nodes
  are held fixed.  This creates local bowing and an evolving pin reaction force.
  The pin_amplification_history.csv output is intended to verify whether each
  pin loads with time after capture.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np

KB_EV = 8.617333262145e-5
EV_J = 1.602176634e-19


@dataclass
class ExpFitBarrier:
    scale: float
    entropy_kB: float
    floor_frac: float
    T0_K: float
    G00_eV: float
    gT: float
    sigc0_MPa: float
    sT: float
    a: float
    n: float

    def G0_eV(self, T: float) -> float:
        return self.scale * self.G00_eV * math.exp(-self.gT * (T - self.T0_K) / self.T0_K)

    def sigc_MPa(self, T: float) -> float:
        return self.sigc0_MPa * math.exp(-self.sT * (T - self.T0_K) / self.T0_K)

    def H_eV(self, tau_MPa: float, T: float) -> float:
        tau = max(float(tau_MPa), 0.0)
        G0 = self.G0_eV(T)
        Gf = self.floor_frac * G0
        sigc = max(self.sigc_MPa(T), 1.0e-300)
        r = tau / sigc
        return Gf + (G0 - Gf) * math.exp(-self.a * (r ** self.n))

    def DeltaG_eV(self, tau_MPa: float, T: float) -> float:
        return max(0.0, self.H_eV(tau_MPa, T) - KB_EV * float(T) * self.entropy_kB)

    def DeltaG_eV_vec(self, tau_MPa: np.ndarray, T: float) -> np.ndarray:
        """Vectorized DeltaG over a numpy array of stresses (for per-node glide)."""
        tau = np.maximum(np.asarray(tau_MPa, dtype=float), 0.0)
        G0 = self.G0_eV(T)
        Gf = self.floor_frac * G0
        sigc = max(self.sigc_MPa(T), 1.0e-300)
        r = tau / sigc
        H = Gf + (G0 - Gf) * np.exp(-self.a * (r ** self.n))
        return np.maximum(0.0, H - KB_EV * float(T) * self.entropy_kB)

    def H_force_eV(self, F_N: float, T: float, Fc_N: float) -> float:
        """EXP-floor barrier as a function of generalized junction force.

        v11 force/work variant.  The DDD variable is the pin reaction force F,
        with conjugate one-Burgers-vector work W = F*b.  The ratio F/Fc is
        equivalent to W/Wc, where Fc is tied outside this class to a physical
        scale such as line tension or to the stress-fit work relation.
        """
        F = max(float(F_N), 0.0)
        Fc = max(float(Fc_N), 1.0e-300)
        G0 = self.G0_eV(T)
        Gf = self.floor_frac * G0
        r = F / Fc
        return Gf + (G0 - Gf) * math.exp(-self.a * (r ** self.n))

    def DeltaG_force_eV(self, F_N: float, T: float, Fc_N: float) -> float:
        """Activation barrier using force/work as the drive variable."""
        return max(0.0, self.H_force_eV(F_N, T, Fc_N) - KB_EV * float(T) * self.entropy_kB)

    def force_scale_from_stress_work_N(self, tau_ref_MPa: float, T: float, b_m: float) -> float:
        """Force scale from the work-conjugacy relation F*b = tau*v*(tau)."""
        tau = max(float(tau_ref_MPa), 0.0)
        return tau * 1.0e6 * self.vstar_m3(tau, T) / max(float(b_m), 1.0e-300)

    def vstar_m3(self, tau_MPa: float, T: float) -> float:
        """v* = -dH/dtau in m^3.

        H is in eV and tau is supplied in MPa.
        dH/dtau_Pa = dH/dtau_MPa / 1e6, then eV/Pa -> J/Pa = m^3.
        """
        tau = max(float(tau_MPa), 1.0e-12)
        G0 = self.G0_eV(T)
        Gf = self.floor_frac * G0
        sigc = max(self.sigc_MPa(T), 1.0e-300)
        r = tau / sigc
        if G0 <= Gf:
            return 0.0
        dH_dtau_MPa = - (G0 - Gf) * math.exp(-self.a * r**self.n) * self.a * self.n * r**(self.n - 1.0) / sigc
        return max(0.0, -dH_dtau_MPa * EV_J / 1.0e6)

    def activation_length_m(self, tau_ref_MPa: float, b_m: float, T: float) -> float:
        """Characteristic activation length ell_act = v*(tau_ref)/b^2  [m].

        Evaluated at a FIXED reference stress (not self-consistently at tau_local),
        so the force->stress map below is explicit and monotone in F.
        """
        v = self.vstar_m3(tau_ref_MPa, T)
        return v / (b_m * b_m)

    def tau_local_direct_MPa(self, F_N: float, b_m: float, tau_ref_MPa: float, T: float,
                             ell_fixed_red: float = 0.0,
                             tau_cap_MPa: float = 8.0e3) -> Tuple[float, float, float]:
        """Direct force-per-area concentration:  tau_local = F_pin / (b * ell_act).

        ell_act comes from the activation volume at a fixed reference stress
        (or a user-fixed length in units of b).  This replaces the implicit,
        non-monotone solve of F*b = tau*v*(tau) and CANNOT saturate at the
        tau*v* peak.  Returns (tau_local_MPa, ell_act_m, vstar_at_ref_m3).
        """
        v_ref = self.vstar_m3(tau_ref_MPa, T)
        if ell_fixed_red and ell_fixed_red > 0.0:
            ell = ell_fixed_red * b_m
        else:
            ell = max(v_ref / (b_m * b_m), b_m)  # floor at one Burgers vector
        tau_Pa = max(float(F_N), 0.0) / (b_m * ell)
        tau_MPa = min(tau_Pa / 1.0e6, float(tau_cap_MPa))
        return tau_MPa, ell, v_ref

    def tau_from_force_MPa(self, F_N: float, b_m: float, T: float, tau_max_MPa: float = 2.0e4) -> Tuple[float, float, float]:
        """Solve tau*vstar(tau,T)=F*b.

        Returns:
          tau_MPa, vstar_m3, tau*vstar in eV.
        The EXP v*(tau) can be non-monotonic, so this uses a robust grid search.
        """
        target_eV = max(float(F_N), 0.0) * b_m / EV_J
        if target_eV <= 0.0:
            v = self.vstar_m3(1.0e-9, T)
            return 0.0, v, 0.0

        def work_eV(tau_MPa: float) -> float:
            return tau_MPa * 1.0e6 * self.vstar_m3(tau_MPa, T) / EV_J

        grid = np.r_[np.linspace(0.0, 10.0, 80), np.logspace(1.0, math.log10(tau_max_MPa), 360)]
        vals = np.array([work_eV(float(t)) for t in grid])
        if not np.any(np.isfinite(vals)):
            tau = tau_max_MPa
            v = self.vstar_m3(tau, T)
            return tau, v, tau * 1.0e6 * v / EV_J

        idx = int(np.nanargmin(np.abs(vals - target_eV)))
        lo = float(grid[max(idx - 2, 0)])
        hi = float(grid[min(idx + 2, len(grid) - 1)])
        fine = np.linspace(lo, hi, 120)
        vals2 = np.array([work_eV(float(t)) for t in fine])
        tau = float(fine[int(np.nanargmin(np.abs(vals2 - target_eV)))])
        v = self.vstar_m3(tau, T)
        return tau, v, tau * 1.0e6 * v / EV_J


def prob_from_rate(rate_s: float, dt: float) -> float:
    x = max(0.0, float(rate_s)) * float(dt)
    return 1.0 if x > 50.0 else 1.0 - math.exp(-x)


def minimum_image_delta(dx: np.ndarray | float, L: float):
    return (dx + 0.5 * L) % L - 0.5 * L


def periodic_distance_x(x_obs: np.ndarray, x_node: float, Lx: float) -> np.ndarray:
    return np.abs(minimum_image_delta(x_obs - x_node, Lx))


def generate_obstacles(rng: np.random.Generator, rho_m2: float, Lx_red: float, Lz_red: float, b_m: float, min_count: int = 0):
    area_m2 = (Lx_red * b_m) * (Lz_red * b_m)
    nobs = max(int(round(rho_m2 * area_m2)), int(min_count))
    if nobs <= 0:
        return np.zeros(0), np.zeros(0)
    return rng.random(nobs) * Lx_red, rng.random(nobs) * Lz_red


def local_feed_length_reduced(z_line: np.ndarray, pinned_mask: np.ndarray, j: int, Lz: float) -> float:
    """Distance to nearest other pinned node along the same line; Lz if isolated."""
    pins = np.where(pinned_mask)[0]
    pins = pins[pins != j]
    dz_node = float(Lz) / max(len(z_line), 1)
    if len(pins) == 0:
        return float(Lz)
    dz = np.abs(z_line[pins] - z_line[j])
    dz = np.minimum(dz, Lz - dz)
    return max(float(np.min(dz)), dz_node)


def line_tension_back_stress_MPa(x_lines: np.ndarray, line_tension_N: float, b_m: float,
                                 Lx: float, dz_node_m: float, cap_MPa: float) -> np.ndarray:
    """Local line-tension back stress resolved on the glide plane, per node.

    The restoring force per unit length on a bowed line is T*kappa, so the
    resolved back stress is tau_back = (T/b) * kappa, with the local curvature
    kappa = d2x/dz2 (x in meters, z along the line).  Discretely,

        d2x_red = (x_{j-1} - x_j) + (x_{j+1} - x_j)            [reduced units]
        kappa   = d2x_red * b / dz_node_m^2                     [1/m]
        tau_back = -(T/b) * kappa                               [Pa]

    A node bowed FORWARD of its neighbors (x_j larger) has kappa < 0, giving a
    positive tau_back that opposes glide; a node lagging its neighbors gets a
    negative tau_back that line tension uses to pull it forward.  v5: the value
    is SIGNED, clipped to [-cap, +cap] in magnitude only, so line tension drives
    relaxation in both directions (the cap is a stability bound on stiff steps).
    Minimum-image handles the x periodicity and the periodic z-closure.
    """
    xm = np.roll(x_lines, 1, axis=1)
    xp = np.roll(x_lines, -1, axis=1)
    d2x_red = minimum_image_delta(xm - x_lines, Lx) + minimum_image_delta(xp - x_lines, Lx)
    kappa = (d2x_red * b_m) / (dz_node_m * dz_node_m)
    tau_back_MPa = (-(line_tension_N / b_m) * kappa) / 1.0e6
    return np.clip(tau_back_MPa, -float(cap_MPa), float(cap_MPa))


def line_tension_reaction_N(x_line: np.ndarray, z_line: np.ndarray, j: int, line_tension_N: float, Lx: float) -> float:
    """Resolved glide-direction line-tension reaction at pinned node j."""
    n = len(x_line)
    jm = (j - 1) % n
    jp = (j + 1) % n

    def vec(a: int, b: int) -> np.ndarray:
        dx = minimum_image_delta(x_line[b] - x_line[a], Lx)
        dz = z_line[b] - z_line[a]
        # z is fixed and ordered but periodic line closure has one long segment;
        # minimum image gives more stable tangent at closure.
        dz = minimum_image_delta(dz, z_line[-1] + (z_line[1] - z_line[0] if n > 1 else 1.0))
        return np.array([dx, dz], dtype=float)

    vm = vec(j, jm)
    vp = vec(j, jp)
    nm = np.linalg.norm(vm)
    npv = np.linalg.norm(vp)
    if nm <= 0.0 or npv <= 0.0:
        return 0.0
    # Sum of unit tensions on the pinned node; x component drives local bypass.
    return float(line_tension_N * abs((vm / nm + vp / npv)[0]))


def line_tension_relaxation(x_lines: np.ndarray, pinned: np.ndarray, Lx: float, relax: float, substeps: int) -> np.ndarray:
    """Connected-line smoothing with pinned nodes held fixed.

    Returns dx due to relaxation in reduced units with minimum-image differences.
    """
    if relax <= 0.0 or substeps <= 0:
        return np.zeros_like(x_lines)

    total_dx = np.zeros_like(x_lines)
    for _ in range(substeps):
        old = x_lines.copy()
        for li in range(x_lines.shape[0]):
            x = x_lines[li]
            left = x + minimum_image_delta(np.roll(x, 1) - x, Lx)
            right = x + minimum_image_delta(np.roll(x, -1) - x, Lx)
            smooth_target = 0.5 * (left + right)
            dx = relax * minimum_image_delta(smooth_target - x, Lx)
            dx[pinned[li]] = 0.0
            x_lines[li] = (x + dx) % Lx
            total_dx[li] += dx
    return total_dx


def swept_strain_increment(dx_lines: np.ndarray, Lx: float, Lz: float, b_m: float, out_spacing_m: float) -> float:
    if dx_lines.size == 0:
        return 0.0
    # Integrate dx*dz over all lines.  v5: SIGNED -- with a detailed-balance (signed)
    # mobility, net reverse glide must be able to recover plastic strain, so we no
    # longer clamp the increment to be non-negative.
    dz = Lz / max(dx_lines.shape[1], 1)
    swept_red2 = float(np.sum(dx_lines)) * dz
    return swept_red2 / max(Lx * Lz, 1e-300) * b_m / max(out_spacing_m, 1e-300)




def free_glide_required_net_rate_s(strain_rate_s: float, Lx_red: float, b_m: float,
                                   out_spacing_m: float, nline: int,
                                   glide_jump_length_red: float) -> float:
    """Net Peierls attempt rate needed to carry the imposed strain rate if all nodes are free.

    For uniform glide, swept_strain_rate = (nline/Lx) * (b/s_out) * dxdt_red,
    and dxdt_red = glide_jump_length_red * net_rate.
    This is a diagnostic only; it does not feed back into the simulation.
    """
    denom = max(float(nline) * b_m * max(float(glide_jump_length_red), 1e-300), 1e-300)
    return max(float(strain_rate_s), 0.0) * float(Lx_red) * float(out_spacing_m) / denom


def peierls_net_rate_s(barrier: ExpFitBarrier, tau_MPa: float, T: float,
                       attempt_frequency_s: float, prefactor: float) -> float:
    """Signed detailed-balance Peierls net rate for a positive applied stress."""
    kT = max(KB_EV * float(T), 1e-300)
    tau = max(float(tau_MPa), 0.0)
    G_fwd = barrier.DeltaG_eV(tau, T)
    G_rev = barrier.DeltaG_eV(0.0, T)
    return float(attempt_frequency_s) * float(prefactor) * (math.exp(-G_fwd / kT) - math.exp(-G_rev / kT))


def solve_peierls_baseline_MPa(barrier: ExpFitBarrier, T: float, required_net_rate_s: float,
                               attempt_frequency_s: float, prefactor: float,
                               tau_max_MPa: float = 1.0e4) -> float:
    """Stress at which a completely free line would carry the imposed strain rate.

    This quantifies the Peierls/friction baseline that can mask forest-density effects.
    Returns NaN if the requested rate cannot be bracketed.
    """
    target = max(float(required_net_rate_s), 0.0)
    if target <= 0.0:
        return 0.0
    lo, hi = 0.0, float(tau_max_MPa)
    if peierls_net_rate_s(barrier, hi, T, attempt_frequency_s, prefactor) < target:
        return float("nan")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if peierls_net_rate_s(barrier, mid, T, attempt_frequency_s, prefactor) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def make_barrier(args, scale: float, entropy: float, sigc0_MPa: float) -> ExpFitBarrier:
    return ExpFitBarrier(
        scale=scale,
        entropy_kB=entropy,
        floor_frac=args.expfit_floor_frac,
        T0_K=args.expfit_T0_K,
        G00_eV=args.expfit_G00_eV,
        gT=args.expfit_gT,
        sigc0_MPa=sigc0_MPa,
        sT=args.expfit_sT,
        a=args.expfit_a,
        n=args.expfit_n,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Clean v8 explicit-obstacle Arrhenius Taylor DDD test.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--temperature-K", type=float, required=True)
    ap.add_argument("--strain-rate", type=float, required=True)
    ap.add_argument("--target-strain", type=float, default=2e-4)
    ap.add_argument("--dt", type=float, default=5e-10)
    ap.add_argument("--steps", type=int, default=0)
    ap.add_argument("--seed", type=int, default=11)

    ap.add_argument("--b-m", type=float, default=2.48e-10)
    ap.add_argument("--mu-Pa", type=float, default=80e9)
    ap.add_argument("--elastic-modulus-MPa", type=float, default=1e5)
    ap.add_argument("--line-tension-factor", type=float, default=0.5)
    ap.add_argument("--line-tension-relax", type=float, default=-1.0,
                    help="<=0 derives relax from line tension + drag; >0 overrides manually.")
    ap.add_argument("--line-tension-relax-substeps", type=int, default=3)

    ap.add_argument("--cell-lx-reduced", type=float, default=3060)
    ap.add_argument("--cell-lz-reduced", type=float, default=1530)
    ap.add_argument("--mobile-line-count", type=int, default=4)
    ap.add_argument("--mobile-line-nodes", type=int, default=512)
    ap.add_argument("--initial-mobile-x-fraction", type=float, default=0.25)

    ap.add_argument("--forest-rho-m2", type=float, required=True)
    ap.add_argument("--forest-min-count", type=int, default=0)
    ap.add_argument("--capture-radius-reduced", type=float, default=4.0)
    ap.add_argument("--capture-mode", choices=["radius", "swept_crossing"], default="radius",
                    help="radius: legacy nearest-node proximity capture; swept_crossing: capture only "
                         "when a node sweeps across an obstacle x-coordinate during driven glide. "
                         "Use swept_crossing to avoid artificial t=0 overlap pinning.")
    ap.add_argument("--snap-swept-capture-to-obstacle", action="store_true",
                    help="with --capture-mode swept_crossing, place the captured lead node at "
                         "the obstacle x-position. This is usually the physical choice for an "
                         "explicit point obstacle.")
    ap.add_argument("--rearm-radius-reduced", type=float, default=-1.0,
                    help="v7: a crossed (line,obstacle) pair can re-pin once the lead node "
                         "has glided this far (reduced units) from the obstacle in x. Must be "
                         "> capture radius. <=0 defaults to 2x capture radius. Makes the "
                         "periodic forest PERSISTENT instead of burned-through-once.")
    ap.add_argument("--min-pin-age-steps", type=int, default=5,
                    help="v7: a freshly captured pin cannot depin until it has aged this many "
                         "steps, so it can load by line bowing before crossing is attempted.")
    ap.add_argument("--allow-same-step-depin", action="store_true",
                    help="diagnostic override. By default v8 uses operator splitting: a pin "
                         "created by capture is not allowed to cross in the same timestep, "
                         "but it can cross on the next step if --min-pin-age-steps allows it.")

    ap.add_argument("--attempt-frequency-s", type=float, default=1e12)
    ap.add_argument("--mobile-glide-prefactor", type=float, default=1.0)
    ap.add_argument("--max-free-dx-reduced", type=float, default=0.5,
                    help="v5/v9: pure CFL clip on |dx| per step. Keep this well below the "
                         "node spacing (typically 0.25--1 b) when backstress mobility is on; "
                         "do NOT tie it to the capture radius.")
    ap.add_argument("--glide-jump-length-reduced", type=float, default=1.0,
                    help="v5: physical activated glide distance per net attempt "
                         "(dx = jump_length * net_rate * dt). Fixed; must NOT depend on "
                         "forest spacing or mesh, so forest density stays out of the "
                         "local Peierls prefactor.")
    ap.add_argument("--phi-denominator-floor-MPa", type=float, default=1.0,
                    help="floor on |tau_app| when forming phi=tau_local/tau_app, so phi "
                         "does not blow up near zero applied stress.")

    ap.add_argument("--out-of-plane-spacing-mode", choices=["fixed", "forest_spacing", "mobile_count"], default="fixed")
    ap.add_argument("--out-of-plane-spacing-m", type=float, default=1e-6)

    ap.add_argument("--expfit-T0-K", type=float, default=1100.0)
    ap.add_argument("--expfit-G00-eV", type=float, default=1.908192)
    ap.add_argument("--expfit-gT", type=float, default=1.241743865563325)
    ap.add_argument("--expfit-sigc0-MPa", type=float, default=1497.042242375928,
                    help="shared characteristic stress; used by a branch only if its "
                         "dedicated --expfit-{cross,peierls}-sigc0-MPa is not set.")
    ap.add_argument("--expfit-cross-sigc0-MPa", type=float, default=-1.0,
                    help="v6: crossing (forest) characteristic stress. Set HIGH so pins "
                         "hold and load up before depinning (Taylor strength scale).")
    ap.add_argument("--expfit-peierls-sigc0-MPa", type=float, default=-1.0,
                    help="v6: Peierls (lattice friction) characteristic stress. Set LOW "
                         "(~1/10 of crossing) so the bulk flows at low stress.")
    ap.add_argument("--expfit-sT", type=float, default=0.10850578873777168)
    ap.add_argument("--expfit-a", type=float, default=2.2056211004282904)
    ap.add_argument("--expfit-n", type=float, default=2.5207319790155385)
    ap.add_argument("--expfit-floor-frac", type=float, default=0.0)
    ap.add_argument("--expfit-cross-scale", type=float, default=0.40)
    ap.add_argument("--expfit-cross-entropy-kB", type=float, default=-9.0)
    ap.add_argument("--expfit-peierls-scale", type=float, default=0.02)
    ap.add_argument("--expfit-peierls-entropy-kB", type=float, default=-9.0)

    ap.add_argument("--pileup-force-mode", choices=["none", "same_obstacle_global"], default="same_obstacle_global")
    ap.add_argument("--max-pileup-contributors", type=int, default=16)

    # --- v3: force -> local-stress concentration -------------------------------
    ap.add_argument("--tau-local-mode", choices=["direct", "workrelation"], default="direct",
                    help="direct: tau_local=F_pin/(b*ell_act) (monotone, no saturation). "
                         "workrelation: legacy implicit solve of F*b=tau*v*(tau) (can saturate).")
    ap.add_argument("--crossing-drive-mode", choices=["local_stress", "force_work"], default="local_stress",
                    help="v11: local_stress uses the v10 tau_local stress barrier. "
                         "force_work uses pin reaction work F_pin*b as the Arrhenius drive variable.")
    ap.add_argument("--cross-force-scale-mode", choices=["line_tension", "stress_work", "manual"], default="line_tension",
                    help="v11 force_work: physical force scale Fc. line_tension uses alpha*T_line; "
                         "stress_work uses tau_ref*v*(tau_ref)/b from the EXP fit; manual uses --cross-force-scale-N.")
    ap.add_argument("--cross-force-scale-factor", type=float, default=1.0,
                    help="multiplicative factor for line_tension or stress_work force scale; independent of density.")
    ap.add_argument("--cross-force-scale-N", type=float, default=-1.0,
                    help="manual v11 force scale for force_work mode; used only when --cross-force-scale-mode manual.")
    ap.add_argument("--activation-ref-stress-MPa", type=float, default=-1.0,
                    help="reference stress for ell_act=v*(tau_ref)/b^2; <=0 uses sigma_c(T).")
    ap.add_argument("--activation-length-reduced", type=float, default=0.0,
                    help="fixed activation length in units of b; >0 overrides v*-derived ell_act.")
    ap.add_argument("--tau-local-cap-MPa", type=float, default=8.0e3,
                    help="hard cap on tau_local (~theoretical strength) to prevent runaway.")
    ap.add_argument("--max-feed-length-reduced", type=float, default=-1.0,
                    help="cap on F_PK feed length; <=0 uses forest spacing 1/sqrt(rho)/b.")

    # --- v3: line-tension relaxation tied to physical line tension -------------
    ap.add_argument("--line-tension-drag-Pa-s", type=float, default=1.0e-5,
                    help="dislocation drag B used to derive the relaxation coefficient.")
    ap.add_argument("--line-tension-relax-max", type=float, default=0.5,
                    help="stability clamp on the derived relaxation coefficient.")

    # --- v4: line-tension BACK STRESS in the Peierls glide law -----------------
    ap.add_argument("--backstress-mobility", choices=["on", "off"], default="on",
                    help="on: free-node Peierls glide is driven by the LOCAL effective "
                         "stress tau_eff=max(0,tau_app-tau_back), tau_back=(T/b)*kappa, so "
                         "bowing against pins slows/arrests the segment (force balance). "
                         "off: recover v3 behaviour (uniform glide from tau_app).")
    ap.add_argument("--max-back-stress-MPa", type=float, default=5.0e3,
                    help="cap on the per-node line-tension back stress (stability).")
    ap.add_argument("--backstress-com-projection", choices=["external_drive", "none"], default="external_drive",
                    help="v10: with backstress mobility on, project the internal line-tension \
                         self-force so it cannot create a spurious net center-of-mass friction \
                         on a free periodic line. external_drive preserves the line-mean \
                         glide set by the applied Peierls stress while retaining curvature-driven \
                         relative node motion. none reproduces the raw v9 local tau_app-tau_back law.")
    ap.add_argument("--project-backstress-on-pinned-lines", action="store_true",
                    help="v12 diagnostic/physics switch. By default the external-drive COM projection "
                         "is applied only to lines with no live pins. A pinned line keeps the raw "
                         "local tau_app - tau_back mobility, so curvature/line-length work can "
                         "reduce net swept area. Set this flag to recover old v10/v11 behavior.")
    ap.add_argument("--line-tension-smooth-frac", type=float, default=0.05,
                    help="with backstress on, the Laplacian step is demoted to a LIGHT "
                         "numerical regularizer with this coefficient (avoids double-"
                         "counting T, which now lives in the velocity law). If "
                         "--line-tension-relax>0 is given, that value overrides this.")

    ap.add_argument("--plastic-strain-source", choices=["free_glide", "total"], default="free_glide",
                    help="free_glide: only driven forward glide sweeps plastic area "
                         "(line-tension relaxation does not contaminate eps_p).")

    ap.add_argument("--pin-diagnostic-every", type=int, default=0)
    ap.add_argument("--pin-diagnostic-max-rows", type=int, default=200000)

    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    T = float(args.temperature_K)
    b = float(args.b_m)
    Lx = float(args.cell_lx_reduced)
    Lz = float(args.cell_lz_reduced)
    nline = int(args.mobile_line_count)
    nn = int(args.mobile_line_nodes)
    nsteps = args.steps if args.steps and args.steps > 0 else int(math.ceil(args.target_strain / (args.strain_rate * args.dt)))

    # v6: Peierls (lattice friction) and crossing (forest obstacle strength) are
    # PHYSICALLY DIFFERENT stress scales and must be decoupled.  Each falls back to
    # the shared --expfit-sigc0-MPa if its dedicated value is not given (<=0).
    cross_sigc0 = args.expfit_cross_sigc0_MPa if args.expfit_cross_sigc0_MPa > 0 else args.expfit_sigc0_MPa
    peierls_sigc0 = args.expfit_peierls_sigc0_MPa if args.expfit_peierls_sigc0_MPa > 0 else args.expfit_sigc0_MPa
    cross = make_barrier(args, args.expfit_cross_scale, args.expfit_cross_entropy_kB, cross_sigc0)
    peierls = make_barrier(args, args.expfit_peierls_scale, args.expfit_peierls_entropy_kB, peierls_sigc0)

    obs_x, obs_z = generate_obstacles(rng, args.forest_rho_m2, Lx, Lz, b, args.forest_min_count)
    nobs = len(obs_x)

    x_lines = np.full((nline, nn), args.initial_mobile_x_fraction * Lx, dtype=float)
    z_base = np.linspace(0.0, Lz, nn, endpoint=False)
    z_lines = np.vstack([(z_base + 0.37 * li * Lz / max(nline, 1) / max(nn, 1)) % Lz for li in range(nline)])

    pinned = np.zeros((nline, nn), dtype=bool)
    pinned_obs = np.full((nline, nn), -1, dtype=int)
    # v7: instead of a PERMANENT "passed" set (which burned the forest through once and
    # left friction-only flow), a (line,obstacle) pair is only BLOCKED from re-capture
    # until the lead node has glided > rearm_radius away in x.  Maps pair -> node id to
    # monitor.  This makes the periodic forest persistent across wraps.
    blocked_until_far: Dict[Tuple[int, int], int] = {}
    rearm_radius_red = (args.rearm_radius_reduced if args.rearm_radius_reduced > 0
                        else 2.0 * args.capture_radius_reduced)
    if rearm_radius_red <= args.capture_radius_reduced:
        rearm_radius_red = args.capture_radius_reduced * 1.5  # safety: must exceed capture
    pinned_pair_node: Dict[Tuple[int, int], int] = {}
    queued_pair_count: Dict[Tuple[int, int], int] = {}
    pin_birth_step: Dict[Tuple[int, int], int] = {}
    pin_birth_tau_MPa: Dict[Tuple[int, int], float] = {}

    if args.out_of_plane_spacing_mode == "forest_spacing":
        s_out = 1.0 / math.sqrt(max(args.forest_rho_m2, 1e-300))
    elif args.out_of_plane_spacing_mode == "mobile_count":
        s_out = Lz * b / max(nline, 1)
    else:
        s_out = float(args.out_of_plane_spacing_m)

    line_tension_N = args.line_tension_factor * args.mu_Pa * b * b

    # v3/v4: choose the Laplacian smoother coefficient.
    #   - manual override (--line-tension-relax>0) always wins;
    #   - with v4 back-stress mobility ON, line tension lives in the velocity law,
    #     so the Laplacian is demoted to a LIGHT regularizer (smooth_frac);
    #   - with back-stress OFF (v3 mode), the smoother is the shape driver and uses
    #     the physically-derived, stability-clamped coefficient.
    dz_node_m = (Lz / max(nn, 1)) * b
    B_drag = max(args.line_tension_drag_Pa_s, 1e-300)
    relax_raw = 2.0 * line_tension_N * args.dt / (B_drag * dz_node_m * dz_node_m)
    relax_derived = min(relax_raw, float(args.line_tension_relax_max))
    if args.line_tension_relax and args.line_tension_relax > 0.0:
        relax_eff = float(args.line_tension_relax)
        relax_source = "manual"
    elif args.backstress_mobility == "on":
        relax_eff = float(args.line_tension_smooth_frac)
        relax_source = f"light_regularizer(smooth_frac={relax_eff:.3f})"
    else:
        relax_eff = relax_derived
        relax_source = f"derived(raw={relax_raw:.3e},clamped={relax_eff:.3f})"

    # v3: feeding length for the PK term is bounded by the forest spacing, not the cell.
    forest_spacing_red = (1.0 / math.sqrt(max(args.forest_rho_m2, 1e-300))) / b
    if args.max_feed_length_reduced and args.max_feed_length_reduced > 0.0:
        feed_cap_red = float(args.max_feed_length_reduced)
    else:
        feed_cap_red = min(Lz, forest_spacing_red)

    # reference stress for the activation length
    tau_ref_MPa = args.activation_ref_stress_MPa
    if tau_ref_MPa is None or tau_ref_MPa <= 0.0:
        tau_ref_MPa = cross.sigc_MPa(T)

    # v11: work-conjugate crossing drive.  In force_work mode the barrier is
    # evaluated as DeltaG(F_pin,T), equivalently DeltaG(W=F_pin*b,T).  The force
    # scale is tied to a physical quantity and is independent of rho.
    if args.cross_force_scale_mode == "manual":
        if args.cross_force_scale_N <= 0.0:
            raise ValueError("--cross-force-scale-N must be >0 when --cross-force-scale-mode manual")
        cross_force_scale_N = float(args.cross_force_scale_N)
        cross_force_scale_source = "manual"
    elif args.cross_force_scale_mode == "stress_work":
        base_Fc = cross.force_scale_from_stress_work_N(tau_ref_MPa, T, b)
        cross_force_scale_N = max(args.cross_force_scale_factor * base_Fc, 1.0e-300)
        cross_force_scale_source = f"stress_work(tau_ref={tau_ref_MPa:.6g} MPa,base={base_Fc:.6e} N)"
    else:
        base_Fc = line_tension_N
        cross_force_scale_N = max(args.cross_force_scale_factor * base_Fc, 1.0e-300)
        cross_force_scale_source = f"line_tension(base={base_Fc:.6e} N)"
    cross_work_scale_eV = cross_force_scale_N * b / EV_J

    # v8 diagnostic: quantify the pure Peierls/free-glide stress that would be obtained
    # if every mobile node were free.  If the measured tail stress is close to this
    # value, the density sweep is still friction/Peierls dominated rather than
    # forest-depinning controlled.
    required_free_net_rate_s = free_glide_required_net_rate_s(
        args.strain_rate, Lx, b, s_out, nline, args.glide_jump_length_reduced
    )
    peierls_baseline_tau_MPa = solve_peierls_baseline_MPa(
        peierls, T, required_free_net_rate_s,
        args.attempt_frequency_s, args.mobile_glide_prefactor,
        tau_max_MPa=max(args.max_back_stress_MPa, 1.0e4),
    )

    params = vars(args).copy()
    params.update({
        "nsteps_effective": nsteps,
        "n_obstacles": nobs,
        "fixed_cell_area_m2": (Lx * b) * (Lz * b),
        "out_of_plane_spacing_m_used": s_out,
        "model": "clean_arrhenius_taylor_explicit_obstacles_v12_pinned_line_raw_backstress",
        "v4_backstress_mobility": args.backstress_mobility,
        "v4_max_back_stress_MPa": args.max_back_stress_MPa,
        "v10_backstress_com_projection": args.backstress_com_projection,
        "v12_project_backstress_on_pinned_lines": args.project_backstress_on_pinned_lines,
        "v12_backstress_projection_rule": "project_only_unpinned_lines_unless_flag_set",
        "v4_line_tension_smooth_frac": args.line_tension_smooth_frac,
        "v3_relax_eff": relax_eff,
        "v3_relax_source": relax_source,
        "v3_dz_node_m": dz_node_m,
        "v3_feed_cap_reduced": feed_cap_red,
        "v3_forest_spacing_reduced": forest_spacing_red,
        "v3_tau_ref_MPa": tau_ref_MPa,
        "v3_activation_length_m": cross.activation_length_m(tau_ref_MPa, b, T),
        "v3_tau_local_mode": args.tau_local_mode,
        "v11_crossing_drive_mode": args.crossing_drive_mode,
        "v11_cross_force_scale_mode": args.cross_force_scale_mode,
        "v11_cross_force_scale_factor": args.cross_force_scale_factor,
        "v11_cross_force_scale_N": cross_force_scale_N,
        "v11_cross_force_scale_source": cross_force_scale_source,
        "v11_cross_work_scale_eV": cross_work_scale_eV,
        "v3_plastic_strain_source": args.plastic_strain_source,
        "v5_signed_mobility": True,
        "v5_glide_jump_length_reduced": args.glide_jump_length_reduced,
        "v5_max_free_dx_reduced_cfl": args.max_free_dx_reduced,
        "v5_phi_denominator_floor_MPa": args.phi_denominator_floor_MPa,
        "v6_cross_sigc0_MPa": cross_sigc0,
        "v6_peierls_sigc0_MPa": peierls_sigc0,
        "v6_sigc0_ratio_cross_over_peierls": cross_sigc0 / max(peierls_sigc0, 1e-300),
        "v7_rearm_radius_reduced": rearm_radius_red,
        "v7_min_pin_age_steps": args.min_pin_age_steps,
        "v7_persistent_forest": True,
        "v8_operator_split_no_same_step_depin": (not args.allow_same_step_depin),
        "v8_required_free_glide_net_rate_s": required_free_net_rate_s,
        "v8_peierls_only_baseline_tau_MPa": peierls_baseline_tau_MPa,
        "v8_peierls_effective_barrier_zero_stress_eV": peierls.DeltaG_eV(0.0, T),
        "v8_cross_effective_barrier_zero_stress_eV": cross.DeltaG_eV(0.0, T),
    })
    (outdir / "clean_arrhenius_params.json").write_text(json.dumps(params, indent=2))

    with open(outdir / "fixed_forest_obstacles.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["obs_id", "x_reduced", "z_reduced", "active_initial"])
        for oid in range(nobs):
            w.writerow([oid, obs_x[oid], obs_z[oid], 1])

    hist_cols = [
        "step", "time_s", "eps_total", "eps_plastic", "sigma_MPa", "tau_MPa",
        "forest_rho_actual_m2", "n_obstacles_active", "n_pinned_nodes", "n_crossed_total",
        "d_eps_p", "d_eps_p_swept_area", "line_length_reduced",
        "mean_dx_free_reduced", "mean_dx_relax_reduced", "max_dx_total_reduced",
        "glide_v_reduced_s", "crossing_rate_s", "crossing_barrier_model",
        "cross_expfloor_G0_eV", "cross_expfloor_floor_eV",
        "forest_stress_concentration_phi", "crossing_tau_local_MPa",
        "pin_force_N", "pin_work_eV", "vstar_b3",
        "crossing_drive_mode", "cross_force_scale_N", "cross_work_scale_eV",
        "cross_drive_work_eV", "cross_force_ratio",
        "n_live_pins", "tau_local_median_MPa", "tau_local_p90_MPa", "tau_local_max_MPa",
        "tau_local_uncapped_max_MPa", "frac_tau_local_capped",
        "phi_median", "phi_p90", "phi_max", "F_line_tension_median_N",
        "crossing_rate_max_s", "crossing_expected_events_step",
        "n_capture", "n_depin", "n_candidate_tests",
        "n_free_nodes", "free_node_fraction", "pinned_node_fraction",
        "mean_dx_raw_free_reduced", "mean_dx_after_projection_free_reduced", "dx_app_reduced",
        "n_com_projected_lines", "n_pinned_raw_lines",
        "blocked_pair_count", "tau_back_abs_mean_MPa", "tau_back_abs_p90_MPa",
        "peierls_only_baseline_tau_MPa", "required_free_glide_net_rate_s",
        "step_walltime_s",
    ]
    event_cols = [
        "step", "time_s", "line_id", "node_id", "obs_id", "event",
        "x_node_reduced", "z_node_reduced", "tau_MPa", "tau_local_MPa",
        "rate_s", "barrier_eV", "work_eV", "prefactor",
        "F_pin_N", "F_pk_N", "F_line_tension_N", "vstar_b3",
        "crossing_drive_mode", "cross_force_scale_N", "cross_drive_work_eV", "cross_force_ratio",
        "pileup_contributors", "L_feed_reduced",
    ]
    pin_cols = [
        "step", "time_s", "age_steps", "age_s", "line_id", "node_id", "obs_id",
        "tau_app_MPa", "tau_local_MPa", "phi_local",
        "F_pin_N", "F_pk_N", "F_line_tension_N",
        "crossing_drive_mode", "cross_force_scale_N", "cross_drive_work_eV", "cross_force_ratio",
        "pileup_contributors", "L_feed_reduced",
        "barrier_eV", "rate_s", "p_event", "vstar_b3", "work_eV",
        "x_node_reduced", "z_node_reduced", "birth_tau_MPa",
    ]

    eps_total = 0.0
    eps_p = 0.0
    time_s = 0.0
    n_crossed_total = 0
    pin_rows_written = 0

    with open(outdir / "single_glider_history.csv", "w", newline="") as hf, \
         open(outdir / "single_glider_crossing_events.csv", "w", newline="") as ef, \
         open(outdir / "pin_amplification_history.csv", "w", newline="") as pf:

        hw = csv.DictWriter(hf, fieldnames=hist_cols)
        ew = csv.DictWriter(ef, fieldnames=event_cols)
        pw = csv.DictWriter(pf, fieldnames=pin_cols)
        hw.writeheader()
        ew.writeheader()
        pw.writeheader()

        for step in range(1, nsteps + 1):
            t0 = time.time()
            time_s += args.dt
            eps_total = min(args.target_strain, eps_total + args.strain_rate * args.dt)
            # v5: SIGNED applied stress (no max(0,...) clamp).  If plastic strain
            # overshoots, a reverse elastic stress develops and pulls the system back,
            # which (with signed glide below) prevents the zero-stress ratchet.
            tau_MPa = args.elastic_modulus_MPa * (eps_total - eps_p)
            kT = max(KB_EV * T, 1e-300)
            pref_glide = args.attempt_frequency_s * args.mobile_glide_prefactor

            # Peierls/free glide for unpinned nodes, as a SIGNED detailed-balance process.
            #   tau_eff,j = tau_app - tau_back,j     (tau_back signed, = (T/b)*kappa_j)
            #   dx_j = jump_length * (R_fwd - R_rev) * dt,   R = nu0*pref*exp(-DeltaG(|.|)/kT)
            # At zero effective stress R_fwd = R_rev -> zero net glide (no ratchet).  As a
            # segment bows against its pins, tau_back grows, tau_eff -> 0, the segment
            # arrests; after a depin the curvature relaxes and it springs forward.  The
            # jump length is FIXED (does not carry forest spacing/mesh into the prefactor);
            # --max-free-dx-reduced is only a CFL clip.  --backstress-mobility off -> v3.
            free = ~pinned
            if args.backstress_mobility == "on":
                tau_back_MPa = line_tension_back_stress_MPa(
                    x_lines, line_tension_N, b, Lx, dz_node_m, args.max_back_stress_MPa
                )
            else:
                tau_back_MPa = np.zeros_like(x_lines)
            tau_eff = tau_MPa - tau_back_MPa
            G_fwd = peierls.DeltaG_eV_vec(np.maximum(tau_eff, 0.0), T)
            G_rev = peierls.DeltaG_eV_vec(np.maximum(-tau_eff, 0.0), T)
            net_rate = pref_glide * (np.exp(-G_fwd / kT) - np.exp(-G_rev / kT))
            dx_trial = args.glide_jump_length_reduced * net_rate * args.dt
            dx_all_raw = np.clip(dx_trial, -args.max_free_dx_reduced, args.max_free_dx_reduced)

            # v10: Internal line-tension forces on a periodic line should not create a net
            # center-of-mass friction when no pins are present.  The raw local law
            # tau_eff=tau_app-tau_back can artificially saturate forward and reverse node
            # hops and cancel the mean glide, leaving a high-stress, no-live-pin tail.
            # Project the curvature-driven part to zero mean over the free nodes on each
            # mobile line, then restore the line-mean glide corresponding to the applied
            # Peierls stress alone.  This keeps the local shape/backstress diagnostic and
            # relative node motion, but prevents self-stress from acting as an artificial
            # density-independent Peierls floor.  --backstress-com-projection none
            # reproduces the v9 raw local mobility.
            dx_all = dx_all_raw.copy()
            dx_app = 0.0
            n_com_projected_lines = 0
            n_pinned_raw_lines = 0
            mean_dx_raw_free = float(np.mean(dx_all_raw[free])) if np.any(free) else 0.0
            if args.backstress_mobility == "on" and args.backstress_com_projection == "external_drive":
                G_app_fwd = peierls.DeltaG_eV(max(tau_MPa, 0.0), T)
                G_app_rev = peierls.DeltaG_eV(max(-tau_MPa, 0.0), T)
                net_rate_app = pref_glide * (math.exp(-G_app_fwd / kT) - math.exp(-G_app_rev / kT))
                dx_app = float(np.clip(args.glide_jump_length_reduced * net_rate_app * args.dt,
                                       -args.max_free_dx_reduced, args.max_free_dx_reduced))
                for li_proj in range(nline):
                    fmask = free[li_proj]
                    if not np.any(fmask):
                        continue

                    line_has_live_pin = bool(np.any(pinned[li_proj]))
                    # v12: COM projection is only appropriate for an entirely free
                    # periodic line.  When a line has a live pin, the pin is an
                    # external constraint; curvature/backstress must be allowed to
                    # reduce mean swept-area rate.  Otherwise the old projection
                    # restores Peierls-only COM glide and projects away Taylor hardening.
                    if line_has_live_pin and not args.project_backstress_on_pinned_lines:
                        n_pinned_raw_lines += 1
                        continue

                    mean_raw = float(np.mean(dx_all_raw[li_proj, fmask]))
                    dx_all[li_proj, fmask] = dx_all_raw[li_proj, fmask] - mean_raw + dx_app
                    n_com_projected_lines += 1

                dx_all = np.clip(dx_all, -args.max_free_dx_reduced, args.max_free_dx_reduced)
            mean_dx_after_projection_free = float(np.mean(dx_all[free])) if np.any(free) else 0.0

            dx_free = np.where(free, dx_all, 0.0)
            v_free = float(np.mean(np.abs(dx_free[free])) / args.dt) if np.any(free) else 0.0
            x_before = x_lines.copy()
            x_lines[free] = (x_lines[free] + dx_free[free]) % Lx
            # v9: retain the driven-glide endpoint before relaxation.  Swept obstacle
            # capture should be based on the driven Peierls glide path, not on subsequent
            # line-tension smoothing.
            x_after_free = x_lines.copy()

            # Connected-line relaxation with pins held fixed.  With v4 back-stress
            # mobility ON this is a LIGHT numerical regularizer only (line tension lives
            # in the velocity law above); with it OFF it is the physical shape driver.
            dx_relax = line_tension_relaxation(
                x_lines, pinned, Lx,
                relax=relax_eff,
                substeps=int(args.line_tension_relax_substeps),
            )
            dx_total = dx_free + dx_relax

            # v7: RE-ARM the forest.  A blocked (line,obstacle) pair becomes capturable
            # again once its monitored lead node has glided > rearm_radius away in x, so
            # the periodic forest persists across wraps instead of burning through once.
            if blocked_until_far and nobs > 0:
                for (bli, boid) in list(blocked_until_far.keys()):
                    if boid < 0 or boid >= nobs:
                        blocked_until_far.pop((bli, boid), None)
                        continue
                    jw = blocked_until_far[(bli, boid)]
                    if periodic_distance_x(obs_x[boid], x_lines[bli, jw], Lx) > rearm_radius_red:
                        blocked_until_far.pop((bli, boid), None)

            # Capture after driven line motion.  v9 adds swept_crossing capture to avoid
            # artificial t=0 pinning from mere proximity.  In swept mode, a node can pin
            # only if its driven-glide path crosses an obstacle x-coordinate; z proximity
            # supplies the finite capture width.  This keeps density out of the hazard and
            # prevents the initial configuration from creating a large, nonphysical pinned
            # fraction before any plastic sweep has occurred.
            n_capture = 0
            n_candidate_tests = 0
            queued_pair_count = {k: 0 for k in pinned_pair_node}
            new_pins_this_step = set()
            if nobs > 0 and args.capture_radius_reduced > 0:
                rcap = float(args.capture_radius_reduced)
                for li in range(nline):
                    # node -> nearest eligible obstacle; list entries are (distance/order, node)
                    targets: Dict[int, List[Tuple[float, int]]] = {}
                    for j in range(nn):
                        if pinned[li, j]:
                            continue

                        dz_obs = minimum_image_delta(obs_z - z_lines[li, j], Lz)

                        if args.capture_mode == "swept_crossing":
                            dx_move = float(minimum_image_delta(x_after_free[li, j] - x_before[li, j], Lx))
                            if abs(dx_move) <= 1.0e-14:
                                continue
                            qx = minimum_image_delta(obs_x - x_before[li, j], Lx)
                            if dx_move > 0.0:
                                cands = np.where((qx > 0.0) & (qx <= dx_move) & (np.abs(dz_obs) <= rcap))[0]
                            else:
                                cands = np.where((qx < 0.0) & (qx >= dx_move) & (np.abs(dz_obs) <= rcap))[0]
                            n_candidate_tests += len(cands)
                            if len(cands) == 0:
                                continue
                            # closest first along the swept path, not closest by radial proximity
                            order_vals = np.abs(qx[cands])
                        else:
                            # Legacy proximity mode: nearest unblocked obstacle inside a radius.
                            d2 = periodic_distance_x(obs_x, x_lines[li, j], Lx)**2 + dz_obs**2
                            cands = np.where(d2 <= rcap * rcap)[0]
                            n_candidate_tests += len(cands)
                            if len(cands) == 0:
                                continue
                            order_vals = d2[cands]

                        chosen = None
                        bestd = None
                        for oid0 in cands[np.argsort(order_vals)]:
                            oid = int(oid0)
                            if (li, oid) not in blocked_until_far:
                                chosen = oid
                                if args.capture_mode == "swept_crossing":
                                    bestd = float(abs(minimum_image_delta(obs_x[oid] - x_before[li, j], Lx)))
                                else:
                                    # d2 exists only in radius mode
                                    bestd = float(periodic_distance_x(obs_x[oid:oid+1], x_lines[li, j], Lx)[0]**2 + dz_obs[oid]**2)
                                break
                        if chosen is None:
                            continue
                        targets.setdefault(chosen, []).append((bestd, j))

                    for oid, lst in targets.items():
                        key = (li, oid)
                        lst.sort()
                        if key in pinned_pair_node:
                            # lead already exists; everything here is pile-up this step
                            queued_pair_count[key] = queued_pair_count.get(key, 0) + len(lst)
                        else:
                            # nearest/crossing node becomes the lead junction.  In swept mode,
                            # optional snapping pins the node at the obstacle intersection rather
                            # than after it has overshot the obstacle during the CFL-limited step.
                            _, jlead = lst[0]
                            if args.capture_mode == "swept_crossing" and args.snap_swept_capture_to_obstacle:
                                x_lines[li, jlead] = obs_x[oid]
                            pinned[li, jlead] = True
                            pinned_obs[li, jlead] = oid
                            pinned_pair_node[key] = jlead
                            queued_pair_count[key] = len(lst) - 1
                            pin_birth_step[key] = step
                            pin_birth_tau_MPa[key] = tau_MPa
                            new_pins_this_step.add(key)
                            n_capture += 1
                            ew.writerow({
                                "step": step, "time_s": time_s, "line_id": li, "node_id": jlead,
                                "obs_id": oid, "event": "capture_pin",
                                "x_node_reduced": x_lines[li, jlead], "z_node_reduced": z_lines[li, jlead],
                                "tau_MPa": tau_MPa, "tau_local_MPa": 0.0,
                                "rate_s": 0.0, "barrier_eV": 0.0, "work_eV": 0.0,
                                "prefactor": 1.0, "F_pin_N": 0.0, "F_pk_N": 0.0,
                                "F_line_tension_N": 0.0, "vstar_b3": 0.0,
                                "crossing_drive_mode": args.crossing_drive_mode,
                                "cross_force_scale_N": cross_force_scale_N,
                                "cross_drive_work_eV": 0.0,
                                "cross_force_ratio": 0.0,
                                "pileup_contributors": queued_pair_count[key],
                                "L_feed_reduced": 0.0,
                            })

            n_depin = 0
            max_rate_cross = 0.0
            expected_cross = 0.0
            # v5: collect live-pin populations for proper statistics (median/p90/max),
            # not just the last dictionary entry.
            live_tau_local = []
            live_tau_local_uncapped = []
            live_phi = []
            live_F = []
            live_Flt = []
            live_vb3 = []
            live_drive_work = []
            live_force_ratio = []
            n_tau_capped = 0
            phi_floor = max(args.phi_denominator_floor_MPa, 1.0e-300)

            # Evaluate lead pinned junctions.
            for key, j in list(pinned_pair_node.items()):
                li, oid = key
                if (li, oid) in blocked_until_far or oid < 0 or oid >= nobs:
                    pinned_pair_node.pop(key, None)
                    queued_pair_count.pop(key, None)
                    pin_birth_step.pop(key, None)
                    pin_birth_tau_MPa.pop(key, None)
                    continue
                if not pinned[li, j] or pinned_obs[li, j] != oid:
                    pinned_pair_node.pop(key, None)
                    queued_pair_count.pop(key, None)
                    pin_birth_step.pop(key, None)
                    pin_birth_tau_MPa.pop(key, None)
                    continue

                Lfeed = local_feed_length_reduced(z_lines[li], pinned[li], j, Lz)
                Lfeed = min(Lfeed, feed_cap_red)  # v3: bounded by forest spacing, not the cell
                F_pk = tau_MPa * 1.0e6 * b * (Lfeed * b)
                F_lt = line_tension_reaction_N(x_lines[li], z_lines[li], j, line_tension_N, Lx)

                pile = 1
                if args.pileup_force_mode == "same_obstacle_global":
                    pile = 0
                    for (li2, oid2), _jj2 in pinned_pair_node.items():
                        if oid2 == oid:
                            pile += 1 + queued_pair_count.get((li2, oid2), 0)
                    pile = max(1, pile)
                    if args.max_pileup_contributors > 0:
                        pile = min(pile, args.max_pileup_contributors)

                F_pin = abs(F_lt) + abs(F_pk) * pile
                # Diagnostic local-stress map retained for reporting/compatibility.
                # In v11 force_work mode this stress does NOT drive the crossing barrier;
                # it only shows what the old v10 mapping would have implied.
                tau_local_uncapped = float("nan")
                if args.tau_local_mode == "direct":
                    tau_local, ell_act, vstar = cross.tau_local_direct_MPa(
                        F_pin, b, tau_ref_MPa, T,
                        ell_fixed_red=args.activation_length_reduced,
                        tau_cap_MPa=args.tau_local_cap_MPa,
                    )
                    tau_local_uncapped = (F_pin / (b * max(ell_act, 1e-300))) / 1.0e6
                    stress_work_eV = tau_local * 1.0e6 * vstar / EV_J
                else:
                    tau_local, vstar, stress_work_eV = cross.tau_from_force_MPa(F_pin, b, T)
                    tau_local_uncapped = tau_local

                drive_work_eV = max(F_pin, 0.0) * b / EV_J
                force_ratio = max(F_pin, 0.0) / max(cross_force_scale_N, 1.0e-300)
                if args.crossing_drive_mode == "force_work":
                    Gc = cross.DeltaG_force_eV(F_pin, T, cross_force_scale_N)
                    work_eV = drive_work_eV
                else:
                    Gc = cross.DeltaG_eV(tau_local, T)
                    work_eV = stress_work_eV
                rate_c = args.attempt_frequency_s * math.exp(-Gc / max(KB_EV * T, 1e-300))
                p_event = prob_from_rate(rate_c, args.dt)

                age_steps = step - int(pin_birth_step.get(key, step))
                # v5: phi against a floored |tau_app| so it does not blow up near zero stress
                phi = tau_local / max(abs(tau_MPa), phi_floor)
                vb3 = vstar / (b**3)
                if tau_local_uncapped >= args.tau_local_cap_MPa - 1e-9:
                    n_tau_capped += 1

                if args.pin_diagnostic_every and args.pin_diagnostic_every > 0:
                    if (step % args.pin_diagnostic_every == 0) and pin_rows_written < args.pin_diagnostic_max_rows:
                        pw.writerow({
                            "step": step, "time_s": time_s, "age_steps": age_steps, "age_s": age_steps * args.dt,
                            "line_id": li, "node_id": j, "obs_id": oid,
                            "tau_app_MPa": tau_MPa, "tau_local_MPa": tau_local, "phi_local": phi,
                            "F_pin_N": F_pin, "F_pk_N": F_pk, "F_line_tension_N": F_lt,
                            "crossing_drive_mode": args.crossing_drive_mode,
                            "cross_force_scale_N": cross_force_scale_N,
                            "cross_drive_work_eV": drive_work_eV,
                            "cross_force_ratio": force_ratio,
                            "pileup_contributors": pile, "L_feed_reduced": Lfeed,
                            "barrier_eV": Gc, "rate_s": rate_c, "p_event": p_event,
                            "vstar_b3": vb3, "work_eV": work_eV,
                            "x_node_reduced": x_lines[li, j], "z_node_reduced": z_lines[li, j],
                            "birth_tau_MPa": pin_birth_tau_MPa.get(key, float("nan")),
                        })
                        pin_rows_written += 1

                max_rate_cross = max(max_rate_cross, rate_c)
                expected_cross += p_event
                live_tau_local.append(tau_local)
                live_tau_local_uncapped.append(tau_local_uncapped)
                live_phi.append(phi)
                live_F.append(F_pin)
                live_Flt.append(abs(F_lt))
                live_vb3.append(vb3)
                live_drive_work.append(drive_work_eV)
                live_force_ratio.append(force_ratio)

                can_attempt_depin = (age_steps >= args.min_pin_age_steps)
                if (key in new_pins_this_step) and (not args.allow_same_step_depin):
                    can_attempt_depin = False

                if can_attempt_depin and rng.random() < p_event:
                    blocked_until_far[(li, oid)] = j  # v7: monitor this node for re-arm
                    n_crossed_total += 1
                    n_depin += 1

                    same = (pinned_obs[li] == oid)
                    pinned[li, same] = False
                    pinned_obs[li, same] = -1
                    pinned_pair_node.pop(key, None)
                    queued_pair_count.pop(key, None)
                    pin_birth_step.pop(key, None)
                    pin_birth_tau_MPa.pop(key, None)

                    ew.writerow({
                        "step": step, "time_s": time_s, "line_id": li, "node_id": j,
                        "obs_id": oid, "event": "depin_cross",
                        "x_node_reduced": x_lines[li, j], "z_node_reduced": z_lines[li, j],
                        "tau_MPa": tau_MPa, "tau_local_MPa": tau_local,
                        "rate_s": rate_c, "barrier_eV": Gc, "work_eV": work_eV,
                        "prefactor": 1.0, "F_pin_N": F_pin, "F_pk_N": F_pk,
                        "F_line_tension_N": F_lt, "vstar_b3": vb3,
                        "crossing_drive_mode": args.crossing_drive_mode,
                        "cross_force_scale_N": cross_force_scale_N,
                        "cross_drive_work_eV": drive_work_eV,
                        "cross_force_ratio": force_ratio,
                        "pileup_contributors": pile, "L_feed_reduced": Lfeed,
                    })

            # v3: plastic strain from the DRIVEN forward glide only by default, so the
            # line-tension relaxation (shape response) does not leak into eps_p.
            dx_plastic = dx_free if args.plastic_strain_source == "free_glide" else dx_total
            d_eps_p = swept_strain_increment(dx_plastic, Lx, Lz, b, s_out)
            eps_p += d_eps_p

            # Approximate line length diagnostic.
            line_len = 0.0
            for li in range(nline):
                dxseg = minimum_image_delta(np.roll(x_lines[li], -1) - x_lines[li], Lx)
                dzseg = np.roll(z_lines[li], -1) - z_lines[li]
                dzseg = minimum_image_delta(dzseg, Lz)
                line_len += float(np.sum(np.sqrt(dxseg**2 + dzseg**2)))

            # v5: population statistics over live lead pins (not the last dict entry).
            if live_tau_local:
                ltl = np.asarray(live_tau_local, dtype=float)
                lphi = np.asarray(live_phi, dtype=float)
                ltlu = np.asarray(live_tau_local_uncapped, dtype=float)
                n_live = int(ltl.size)
                tl_med = float(np.median(ltl)); tl_p90 = float(np.quantile(ltl, 0.9)); tl_max = float(np.max(ltl))
                tlu_max = float(np.nanmax(ltlu))
                phi_med = float(np.median(lphi)); phi_p90 = float(np.quantile(lphi, 0.9)); phi_max = float(np.max(lphi))
                Flt_med = float(np.median(np.asarray(live_Flt, dtype=float)))
                F_med = float(np.median(np.asarray(live_F, dtype=float)))
                vb3_med = float(np.median(np.asarray(live_vb3, dtype=float)))
                wd_med = float(np.median(np.asarray(live_drive_work, dtype=float)))
                fr_med = float(np.median(np.asarray(live_force_ratio, dtype=float)))
                frac_capped = float(n_tau_capped) / float(n_live)
            else:
                n_live = 0
                tl_med = tl_p90 = tl_max = tlu_max = 0.0
                phi_med = phi_p90 = phi_max = 0.0
                Flt_med = F_med = vb3_med = 0.0
                wd_med = fr_med = 0.0
                frac_capped = 0.0

            tau_back_abs = np.abs(tau_back_MPa) if 'tau_back_MPa' in locals() else np.zeros_like(x_lines)
            tau_back_abs_mean = float(np.mean(tau_back_abs)) if tau_back_abs.size else 0.0
            tau_back_abs_p90 = float(np.quantile(tau_back_abs, 0.9)) if tau_back_abs.size else 0.0
            n_free_nodes = int(np.sum(~pinned))
            n_total_nodes = int(pinned.size)

            hw.writerow({
                "step": step, "time_s": time_s,
                "eps_total": eps_total, "eps_plastic": eps_p,
                "sigma_MPa": tau_MPa, "tau_MPa": tau_MPa,
                "forest_rho_actual_m2": args.forest_rho_m2,
                "n_obstacles_active": nobs,
                "n_pinned_nodes": int(np.sum(pinned)),
                "n_crossed_total": n_crossed_total,
                "d_eps_p": d_eps_p,
                "d_eps_p_swept_area": d_eps_p,
                "line_length_reduced": line_len,
                "mean_dx_free_reduced": float(np.mean(dx_free)),
                "mean_dx_relax_reduced": float(np.mean(dx_relax)),
                "max_dx_total_reduced": float(np.max(dx_total)),
                "glide_v_reduced_s": v_free,
                "crossing_rate_s": max_rate_cross,
                "crossing_barrier_model": ("expfit_floor_force_work_v11" if args.crossing_drive_mode == "force_work" else "expfit_floor_direct_v5_signed"),
                "cross_expfloor_G0_eV": cross.G0_eV(T),
                "cross_expfloor_floor_eV": cross.floor_frac * cross.G0_eV(T),
                "forest_stress_concentration_phi": phi_med,
                "crossing_tau_local_MPa": tl_med,
                "pin_force_N": F_med,
                "pin_work_eV": wd_med,
                "vstar_b3": vb3_med,
                "crossing_drive_mode": args.crossing_drive_mode,
                "cross_force_scale_N": cross_force_scale_N,
                "cross_work_scale_eV": cross_work_scale_eV,
                "cross_drive_work_eV": wd_med,
                "cross_force_ratio": fr_med,
                "n_live_pins": n_live,
                "tau_local_median_MPa": tl_med,
                "tau_local_p90_MPa": tl_p90,
                "tau_local_max_MPa": tl_max,
                "tau_local_uncapped_max_MPa": tlu_max,
                "frac_tau_local_capped": frac_capped,
                "phi_median": phi_med,
                "phi_p90": phi_p90,
                "phi_max": phi_max,
                "F_line_tension_median_N": Flt_med,
                "crossing_rate_max_s": max_rate_cross,
                "crossing_expected_events_step": expected_cross,
                "n_capture": n_capture,
                "n_depin": n_depin,
                "n_candidate_tests": n_candidate_tests,
                "n_free_nodes": n_free_nodes,
                "free_node_fraction": n_free_nodes / max(n_total_nodes, 1),
                "pinned_node_fraction": 1.0 - n_free_nodes / max(n_total_nodes, 1),
                "mean_dx_raw_free_reduced": mean_dx_raw_free,
                "mean_dx_after_projection_free_reduced": mean_dx_after_projection_free,
                "dx_app_reduced": dx_app,
                "n_com_projected_lines": n_com_projected_lines,
                "n_pinned_raw_lines": n_pinned_raw_lines,
                "blocked_pair_count": len(blocked_until_far),
                "tau_back_abs_mean_MPa": tau_back_abs_mean,
                "tau_back_abs_p90_MPa": tau_back_abs_p90,
                "peierls_only_baseline_tau_MPa": peierls_baseline_tau_MPa,
                "required_free_glide_net_rate_s": required_free_net_rate_s,
                "step_walltime_s": time.time() - t0,
            })

            if eps_total >= args.target_strain - 1e-18:
                break

    with open(outdir / "single_glider_final_nodes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["line_id", "node_id", "x_reduced", "z_reduced", "pinned", "pinned_obs"])
        for li in range(nline):
            for j in range(nn):
                w.writerow([li, j, x_lines[li, j], z_lines[li, j], int(pinned[li, j]), int(pinned_obs[li, j])])

    (outdir / "run.finished").write_text("")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
