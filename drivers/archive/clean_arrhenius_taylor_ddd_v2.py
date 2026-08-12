#!/usr/bin/env python3
"""
clean_arrhenius_taylor_ddd_v2.py

Self-contained explicit-obstacle Arrhenius Taylor / Peierls DDD-style test driver.

This version is intentionally not a patch of previous scripts.  It keeps the
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
    # integrate dx dz over all lines; dx can be negative for relaxation.  Use net swept
    # glide displacement but do not allow negative plastic increment.
    dz = Lz / max(dx_lines.shape[1], 1)
    swept_red2 = float(np.sum(dx_lines)) * dz
    return max(0.0, swept_red2 / max(Lx * Lz, 1e-300) * b_m / max(out_spacing_m, 1e-300))


def make_barrier(args, scale: float, entropy: float) -> ExpFitBarrier:
    return ExpFitBarrier(
        scale=scale,
        entropy_kB=entropy,
        floor_frac=args.expfit_floor_frac,
        T0_K=args.expfit_T0_K,
        G00_eV=args.expfit_G00_eV,
        gT=args.expfit_gT,
        sigc0_MPa=args.expfit_sigc0_MPa,
        sT=args.expfit_sT,
        a=args.expfit_a,
        n=args.expfit_n,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Clean v2 explicit-obstacle Arrhenius Taylor DDD test.")
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
    ap.add_argument("--line-tension-relax", type=float, default=0.25)
    ap.add_argument("--line-tension-relax-substeps", type=int, default=3)

    ap.add_argument("--cell-lx-reduced", type=float, default=3060)
    ap.add_argument("--cell-lz-reduced", type=float, default=1530)
    ap.add_argument("--mobile-line-count", type=int, default=4)
    ap.add_argument("--mobile-line-nodes", type=int, default=512)
    ap.add_argument("--initial-mobile-x-fraction", type=float, default=0.25)

    ap.add_argument("--forest-rho-m2", type=float, required=True)
    ap.add_argument("--forest-min-count", type=int, default=0)
    ap.add_argument("--capture-radius-reduced", type=float, default=4.0)

    ap.add_argument("--attempt-frequency-s", type=float, default=1e12)
    ap.add_argument("--mobile-glide-prefactor", type=float, default=1.0)
    ap.add_argument("--max-free-dx-reduced", type=float, default=2.0)

    ap.add_argument("--out-of-plane-spacing-mode", choices=["fixed", "forest_spacing", "mobile_count"], default="fixed")
    ap.add_argument("--out-of-plane-spacing-m", type=float, default=1e-6)

    ap.add_argument("--expfit-T0-K", type=float, default=1100.0)
    ap.add_argument("--expfit-G00-eV", type=float, default=1.908192)
    ap.add_argument("--expfit-gT", type=float, default=1.241743865563325)
    ap.add_argument("--expfit-sigc0-MPa", type=float, default=1497.042242375928)
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

    cross = make_barrier(args, args.expfit_cross_scale, args.expfit_cross_entropy_kB)
    peierls = make_barrier(args, args.expfit_peierls_scale, args.expfit_peierls_entropy_kB)

    obs_x, obs_z = generate_obstacles(rng, args.forest_rho_m2, Lx, Lz, b, args.forest_min_count)
    nobs = len(obs_x)

    x_lines = np.full((nline, nn), args.initial_mobile_x_fraction * Lx, dtype=float)
    z_base = np.linspace(0.0, Lz, nn, endpoint=False)
    z_lines = np.vstack([(z_base + 0.37 * li * Lz / max(nline, 1) / max(nn, 1)) % Lz for li in range(nline)])

    pinned = np.zeros((nline, nn), dtype=bool)
    pinned_obs = np.full((nline, nn), -1, dtype=int)
    line_obstacle_passed = set()
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

    params = vars(args).copy()
    params.update({
        "nsteps_effective": nsteps,
        "n_obstacles": nobs,
        "fixed_cell_area_m2": (Lx * b) * (Lz * b),
        "out_of_plane_spacing_m_used": s_out,
        "model": "clean_arrhenius_taylor_explicit_obstacles_v2",
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
        "crossing_rate_max_s", "crossing_expected_events_step",
        "n_capture", "n_depin", "n_candidate_tests", "step_walltime_s",
    ]
    event_cols = [
        "step", "time_s", "line_id", "node_id", "obs_id", "event",
        "x_node_reduced", "z_node_reduced", "tau_MPa", "tau_local_MPa",
        "rate_s", "barrier_eV", "work_eV", "prefactor",
        "F_pin_N", "F_pk_N", "F_line_tension_N", "vstar_b3",
        "pileup_contributors", "L_feed_reduced",
    ]
    pin_cols = [
        "step", "time_s", "age_steps", "age_s", "line_id", "node_id", "obs_id",
        "tau_app_MPa", "tau_local_MPa", "phi_local",
        "F_pin_N", "F_pk_N", "F_line_tension_N",
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
            tau_MPa = max(0.0, args.elastic_modulus_MPa * (eps_total - eps_p))

            # Peierls/free glide for unpinned nodes.
            Gp = peierls.DeltaG_eV(tau_MPa, T)
            rate_p = args.attempt_frequency_s * args.mobile_glide_prefactor * math.exp(-Gp / max(KB_EV * T, 1e-300))
            v_free = args.max_free_dx_reduced * rate_p
            dx_free_scalar = min(args.max_free_dx_reduced, v_free * args.dt)

            dx_free = np.zeros_like(x_lines)
            free = ~pinned
            dx_free[free] = dx_free_scalar
            x_before = x_lines.copy()
            x_lines[free] = (x_lines[free] + dx_free[free]) % Lx

            # Connected-line relaxation/bowing with pins held fixed.
            dx_relax = line_tension_relaxation(
                x_lines, pinned, Lx,
                relax=float(args.line_tension_relax),
                substeps=int(args.line_tension_relax_substeps),
            )
            dx_total = dx_free + dx_relax

            # Capture after line motion/relaxation.
            n_capture = 0
            n_candidate_tests = 0
            if nobs > 0 and args.capture_radius_reduced > 0:
                rcap = float(args.capture_radius_reduced)
                for li in range(nline):
                    for j in range(nn):
                        if pinned[li, j]:
                            continue
                        d2 = periodic_distance_x(obs_x, x_lines[li, j], Lx)**2 + (obs_z - z_lines[li, j])**2
                        cands = np.where(d2 <= rcap * rcap)[0]
                        n_candidate_tests += len(cands)
                        if len(cands) == 0:
                            continue
                        chosen = None
                        for oid0 in sorted(cands, key=lambda k: d2[k]):
                            oid = int(oid0)
                            if (li, oid) not in line_obstacle_passed:
                                chosen = oid
                                break
                        if chosen is None:
                            continue

                        key = (li, chosen)
                        pinned[li, j] = True
                        pinned_obs[li, j] = chosen
                        # Snap to obstacle x for a well-defined local junction.
                        x_lines[li, j] = obs_x[chosen]
                        n_capture += 1

                        if key in pinned_pair_node:
                            queued_pair_count[key] = queued_pair_count.get(key, 0) + 1
                            evname = "capture_queued_same_line"
                        else:
                            pinned_pair_node[key] = j
                            queued_pair_count.setdefault(key, 0)
                            pin_birth_step[key] = step
                            pin_birth_tau_MPa[key] = tau_MPa
                            evname = "capture_pin"

                        ew.writerow({
                            "step": step, "time_s": time_s, "line_id": li, "node_id": j,
                            "obs_id": chosen, "event": evname,
                            "x_node_reduced": x_lines[li, j], "z_node_reduced": z_lines[li, j],
                            "tau_MPa": tau_MPa, "tau_local_MPa": 0.0,
                            "rate_s": 0.0, "barrier_eV": 0.0, "work_eV": 0.0,
                            "prefactor": 1.0, "F_pin_N": 0.0, "F_pk_N": 0.0,
                            "F_line_tension_N": 0.0, "vstar_b3": 0.0,
                            "pileup_contributors": queued_pair_count[key],
                            "L_feed_reduced": 0.0,
                        })

            n_depin = 0
            max_rate_cross = 0.0
            expected_cross = 0.0
            last_tau_local = 0.0
            last_phi = 0.0
            last_F = 0.0
            last_work = 0.0
            last_vb3 = 0.0

            # Evaluate lead pinned junctions.
            for key, j in list(pinned_pair_node.items()):
                li, oid = key
                if (li, oid) in line_obstacle_passed or oid < 0 or oid >= nobs:
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
                tau_local, vstar, work_eV = cross.tau_from_force_MPa(F_pin, b, T)
                Gc = cross.DeltaG_eV(tau_local, T)
                rate_c = args.attempt_frequency_s * math.exp(-Gc / max(KB_EV * T, 1e-300))
                p_event = prob_from_rate(rate_c, args.dt)

                age_steps = step - int(pin_birth_step.get(key, step))
                phi = tau_local / max(tau_MPa, 1.0e-300)
                vb3 = vstar / (b**3)

                if args.pin_diagnostic_every and args.pin_diagnostic_every > 0:
                    if (step % args.pin_diagnostic_every == 0) and pin_rows_written < args.pin_diagnostic_max_rows:
                        pw.writerow({
                            "step": step, "time_s": time_s, "age_steps": age_steps, "age_s": age_steps * args.dt,
                            "line_id": li, "node_id": j, "obs_id": oid,
                            "tau_app_MPa": tau_MPa, "tau_local_MPa": tau_local, "phi_local": phi,
                            "F_pin_N": F_pin, "F_pk_N": F_pk, "F_line_tension_N": F_lt,
                            "pileup_contributors": pile, "L_feed_reduced": Lfeed,
                            "barrier_eV": Gc, "rate_s": rate_c, "p_event": p_event,
                            "vstar_b3": vb3, "work_eV": work_eV,
                            "x_node_reduced": x_lines[li, j], "z_node_reduced": z_lines[li, j],
                            "birth_tau_MPa": pin_birth_tau_MPa.get(key, float("nan")),
                        })
                        pin_rows_written += 1

                max_rate_cross = max(max_rate_cross, rate_c)
                expected_cross += p_event
                last_tau_local = tau_local
                last_phi = phi
                last_F = F_pin
                last_work = work_eV
                last_vb3 = vb3

                if rng.random() < p_event:
                    line_obstacle_passed.add((li, oid))
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
                        "pileup_contributors": pile, "L_feed_reduced": Lfeed,
                    })

            d_eps_p = swept_strain_increment(dx_total, Lx, Lz, b, s_out)
            eps_p += d_eps_p

            # Approximate line length diagnostic.
            line_len = 0.0
            for li in range(nline):
                dxseg = minimum_image_delta(np.roll(x_lines[li], -1) - x_lines[li], Lx)
                dzseg = np.roll(z_lines[li], -1) - z_lines[li]
                dzseg = minimum_image_delta(dzseg, Lz)
                line_len += float(np.sum(np.sqrt(dxseg**2 + dzseg**2)))

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
                "crossing_barrier_model": "expfit_floor_direct_v2",
                "cross_expfloor_G0_eV": cross.G0_eV(T),
                "cross_expfloor_floor_eV": cross.floor_frac * cross.G0_eV(T),
                "forest_stress_concentration_phi": last_phi,
                "crossing_tau_local_MPa": last_tau_local,
                "pin_force_N": last_F,
                "pin_work_eV": last_work,
                "vstar_b3": last_vb3,
                "crossing_rate_max_s": max_rate_cross,
                "crossing_expected_events_step": expected_cross,
                "n_capture": n_capture,
                "n_depin": n_depin,
                "n_candidate_tests": n_candidate_tests,
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
