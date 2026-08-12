#!/usr/bin/env python3
"""
Patch clean_arrhenius_taylor_ddd_v13.py -> clean_arrhenius_taylor_ddd_v14.py.

v14 changes:
  1. Preflight kinetic diagnostics for the crossing barrier:
       - zero-force depinning rate/lifetime
       - high-drive/floor crossing rate and R*dt
       - forest spacing / free-sweep time / pin-survival ratio
       - warnings written to preflight_diagnostics.{json,txt}
  2. History stress bookkeeping fix:
       - tau_MPa/sigma_MPa are now written AFTER the plastic strain update,
         consistent with eps_total and eps_plastic in the same row.
       - tau_before_step_MPa and tau_after_step_MPa are both written.
       - d_tau_step_MPa, d_eps_total, d_eps_p_over_d_eps_total, and
         d_eps_p_per_depin_step are added.
  3. End-of-run summaries:
       - run_summary.json
       - run_summary.txt
       - run_summary.csv
     including stress peak/drop, monotonicity, depin/strain coupling, cap usage,
     and event-resolution diagnostics.

This is a diagnostics/accounting patch, not a physics retuning patch.
"""
from __future__ import annotations
from pathlib import Path
import shutil
import subprocess
import sys
import time

src = Path("clean_arrhenius_taylor_ddd_v13.py")
if not src.exists():
    alt = Path("clean_arrhenius_taylor_ddd_v13_actual_swept.py")
    if alt.exists():
        src = alt
    else:
        raise SystemExit("Could not find clean_arrhenius_taylor_ddd_v13.py")

dst = Path("clean_arrhenius_taylor_ddd_v14.py")
backup = src.with_suffix(src.suffix + f".bak_before_v14_{int(time.time())}")
shutil.copy2(src, backup)
print(f"backup: {backup}")
print(f"source: {src}")
print(f"dest:   {dst}")

s = src.read_text()

def fail(msg: str):
    raise SystemExit("ERROR: " + msg)

# -----------------------------------------------------------------------------
# 1. Add helper functions before main.
# -----------------------------------------------------------------------------
helper_marker = "\ndef main(argv=None) -> int:\n"
helper_code = r'''

def _safe_float(x, default=float("nan")) -> float:
    try:
        if x is None:
            return default
        y = float(x)
        return y if math.isfinite(y) else default
    except Exception:
        return default


def _safe_int(x, default=0) -> int:
    try:
        if x is None:
            return default
        return int(float(x))
    except Exception:
        return default


def _rate_from_barrier(deltaG_eV: float, T: float, attempt_frequency_s: float) -> float:
    kT = max(KB_EV * float(T), 1.0e-300)
    if not math.isfinite(deltaG_eV):
        return float("nan")
    x = -float(deltaG_eV) / kT
    if x < -745:
        return 0.0
    if x > 700:
        return float("inf")
    return float(attempt_frequency_s) * math.exp(x)


def write_run_summary(outdir: Path, params: dict, target_strain: float) -> None:
    """Write robust one-run diagnostics from single_glider_history.csv.

    The history CSV may be inspected while still being written or may contain a
    final malformed fragment after an interrupted run.  This reader therefore
    skips malformed rows and keeps only physically possible eps_total values.
    """
    hist_path = outdir / "single_glider_history.csv"
    if not hist_path.exists():
        return

    import csv as _csv

    rows = []
    with open(hist_path, newline="") as f:
        reader = _csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for row in reader:
            # DictReader stores extra fields under key None when a row has too many columns.
            if None in row:
                continue
            eps_total = _safe_float(row.get("eps_total"))
            if not math.isfinite(eps_total):
                continue
            if eps_total < -1e-15 or eps_total > 1.05 * max(float(target_strain), 1e-300):
                continue
            rows.append(row)

    if not rows:
        return

    def arr(name: str):
        return [_safe_float(r.get(name)) for r in rows]

    def finite(vals):
        return [v for v in vals if math.isfinite(v)]

    def median(vals):
        vals = sorted(finite(vals))
        if not vals:
            return float("nan")
        n = len(vals)
        mid = n // 2
        return vals[mid] if n % 2 else 0.5 * (vals[mid - 1] + vals[mid])

    def q(vals, quantile: float):
        vals = sorted(finite(vals))
        if not vals:
            return float("nan")
        if len(vals) == 1:
            return vals[0]
        x = quantile * (len(vals) - 1)
        lo = int(math.floor(x)); hi = int(math.ceil(x))
        if lo == hi:
            return vals[lo]
        return vals[lo] * (hi - x) + vals[hi] * (x - lo)

    final = rows[-1]
    n = len(rows)
    tail_n = max(20, n // 5)
    tail = rows[-tail_n:]

    eps_total = _safe_float(final.get("eps_total"))
    eps_p = _safe_float(final.get("eps_plastic"))
    tau = arr("tau_MPa")
    tau_tail = [_safe_float(r.get("tau_MPa")) for r in tail]
    d_tau = arr("d_tau_step_MPa")
    d_eps_p = arr("d_eps_p")
    d_eps_total = arr("d_eps_total")
    n_depin = [_safe_float(r.get("n_depin"), 0.0) for r in rows]
    n_capture = [_safe_float(r.get("n_capture"), 0.0) for r in rows]
    n_live_tail = [_safe_float(r.get("n_live_pins"), 0.0) for r in tail]
    cap_tail = [_safe_float(r.get("frac_tau_local_capped"), 0.0) for r in tail]
    cap_all = [_safe_float(r.get("frac_tau_local_capped"), 0.0) for r in rows]
    tau_local_tail = [_safe_float(r.get("tau_local_median_MPa"), 0.0) for r in tail]

    dt = float(params.get("dt", float("nan")))
    rdt_vals = []
    for r in rows:
        rate = _safe_float(r.get("crossing_rate_max_s"), float("nan"))
        if math.isfinite(rate) and math.isfinite(dt):
            rdt_vals.append(rate * dt)

    depin_steps = [i for i, v in enumerate(n_depin) if v > 0]
    no_depin_steps = [i for i, v in enumerate(n_depin) if v <= 0]
    d_eps_p_depin = [d_eps_p[i] for i in depin_steps if i < len(d_eps_p)]
    d_eps_p_no_depin = [d_eps_p[i] for i in no_depin_steps if i < len(d_eps_p)]
    d_tau_depin = [d_tau[i] for i in depin_steps if i < len(d_tau)]
    d_tau_no_depin = [d_tau[i] for i in no_depin_steps if i < len(d_tau)]

    tau_finite = finite(tau)
    peak_tau = max(tau_finite) if tau_finite else float("nan")
    peak_idx = tau.index(peak_tau) if tau_finite and peak_tau in tau else -1
    final_tau = _safe_float(final.get("tau_MPa"))
    peak_eps = _safe_float(rows[peak_idx].get("eps_total")) if peak_idx >= 0 else float("nan")
    post_peak_drop = peak_tau - final_tau if math.isfinite(peak_tau) and math.isfinite(final_tau) else float("nan")

    n_tau_decrease = sum(1 for v in d_tau if math.isfinite(v) and v < 0.0)
    n_tau_increase = sum(1 for v in d_tau if math.isfinite(v) and v > 0.0)
    n_plastic_over_imposed = sum(
        1 for dp, de in zip(d_eps_p, d_eps_total)
        if math.isfinite(dp) and math.isfinite(de) and dp > de
    )
    total_depin = sum(v for v in n_depin if math.isfinite(v))
    total_capture = sum(v for v in n_capture if math.isfinite(v))
    epsp_per_depin = eps_p / total_depin if total_depin > 0 and math.isfinite(eps_p) else float("nan")

    summary = {
        "run_dir": str(outdir),
        "model": params.get("model"),
        "T_K": params.get("temperature_K"),
        "rho_m2": params.get("forest_rho_m2"),
        "strain_rate_s": params.get("strain_rate"),
        "dt_s": params.get("dt"),
        "target_strain": target_strain,
        "valid_history_rows": n,
        "eps_total_final": eps_total,
        "eps_plastic_final": eps_p,
        "epsp_over_epstotal_final": eps_p / eps_total if eps_total else float("nan"),
        "tau_final_MPa": final_tau,
        "tau_tail_median_MPa": median(tau_tail),
        "tau_peak_MPa": peak_tau,
        "tau_peak_eps_total": peak_eps,
        "tau_post_peak_drop_MPa": post_peak_drop,
        "n_tau_increase_steps": n_tau_increase,
        "n_tau_decrease_steps": n_tau_decrease,
        "frac_tau_decrease_steps": n_tau_decrease / max(n_tau_decrease + n_tau_increase, 1),
        "n_steps_d_eps_p_gt_d_eps_total": n_plastic_over_imposed,
        "total_depin": total_depin,
        "total_capture": total_capture,
        "epsp_per_depin": epsp_per_depin,
        "n_depin_steps": len(depin_steps),
        "frac_steps_with_depin": len(depin_steps) / max(n, 1),
        "mean_d_eps_p_on_depin_steps": sum(finite(d_eps_p_depin)) / max(len(finite(d_eps_p_depin)), 1),
        "mean_d_eps_p_on_no_depin_steps": sum(finite(d_eps_p_no_depin)) / max(len(finite(d_eps_p_no_depin)), 1),
        "mean_d_tau_on_depin_steps_MPa": sum(finite(d_tau_depin)) / max(len(finite(d_tau_depin)), 1),
        "mean_d_tau_on_no_depin_steps_MPa": sum(finite(d_tau_no_depin)) / max(len(finite(d_tau_no_depin)), 1),
        "n_live_pins_tail_median": median(n_live_tail),
        "tau_local_median_tail_MPa": median(tau_local_tail),
        "frac_tau_local_capped_tail_median": median(cap_tail),
        "frac_tau_local_capped_max": max(finite(cap_all)) if finite(cap_all) else float("nan"),
        "max_crossing_rate_dt": max(finite(rdt_vals)) if finite(rdt_vals) else float("nan"),
        "tail_crossing_rate_dt_median": median(rdt_vals[-tail_n:]) if rdt_vals else float("nan"),
        "preflight_cross_R0_dt": params.get("preflight_cross_R0_dt"),
        "preflight_cross_Rfloor_dt": params.get("preflight_cross_Rfloor_dt"),
        "preflight_sweep_time_over_zero_force_lifetime": params.get("preflight_sweep_time_over_zero_force_lifetime"),
    }

    warnings = []
    if summary["max_crossing_rate_dt"] and math.isfinite(summary["max_crossing_rate_dt"]) and summary["max_crossing_rate_dt"] > 0.2:
        warnings.append("event_resolution: max(crossing_rate*dt) > 0.2; consider smaller dt")
    if summary["frac_tau_local_capped_tail_median"] and math.isfinite(summary["frac_tau_local_capped_tail_median"]) and summary["frac_tau_local_capped_tail_median"] > 0.05:
        warnings.append("local_stress_cap_active_in_tail")
    if summary["n_tau_decrease_steps"] == 0 and summary["total_depin"] > 0:
        warnings.append("stress_never_decreased_despite_depin_events")
    if summary["tau_post_peak_drop_MPa"] and math.isfinite(summary["tau_post_peak_drop_MPa"]) and summary["tau_post_peak_drop_MPa"] < 1e-6:
        warnings.append("no_post_peak_stress_drop")
    summary["warnings"] = warnings

    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    # One-row CSV for easy grep/aggregation.
    with open(outdir / "run_summary.csv", "w", newline="") as f:
        writer = _csv.DictWriter(f, fieldnames=list(summary.keys()))
        writer.writeheader()
        writer.writerow(summary)

    lines = []
    lines.append("Run summary")
    lines.append("===========")
    for key in [
        "T_K", "rho_m2", "strain_rate_s", "dt_s", "eps_total_final",
        "epsp_over_epstotal_final", "tau_tail_median_MPa", "tau_peak_MPa",
        "tau_post_peak_drop_MPa", "total_depin", "n_live_pins_tail_median",
        "tau_local_median_tail_MPa", "frac_tau_local_capped_tail_median",
        "max_crossing_rate_dt", "preflight_cross_Rfloor_dt",
        "preflight_sweep_time_over_zero_force_lifetime",
        "n_tau_decrease_steps", "n_steps_d_eps_p_gt_d_eps_total",
        "mean_d_tau_on_depin_steps_MPa", "mean_d_eps_p_on_depin_steps",
    ]:
        lines.append(f"{key}: {summary.get(key)}")
    if warnings:
        lines.append("warnings:")
        for w in warnings:
            lines.append(f"  - {w}")
    else:
        lines.append("warnings: none")
    (outdir / "run_summary.txt").write_text("\n".join(lines) + "\n")

'''
if 'def write_run_summary(' not in s:
    if helper_marker not in s:
        fail('Could not find main marker')
    s = s.replace(helper_marker, helper_code + helper_marker, 1)

# -----------------------------------------------------------------------------
# 2. Add CLI args.
# -----------------------------------------------------------------------------
old_args = '''    ap.add_argument("--pin-diagnostic-every", type=int, default=0)
    ap.add_argument("--pin-diagnostic-max-rows", type=int, default=200000)
'''
new_args = '''    ap.add_argument("--pin-diagnostic-every", type=int, default=0)
    ap.add_argument("--pin-diagnostic-max-rows", type=int, default=200000)

    ap.add_argument("--preflight-rdt-warn", type=float, default=0.1,
                    help="write a warning if the analytical high-drive/floor crossing rate times dt exceeds this value")
    ap.add_argument("--preflight-rdt-stop", type=float, default=0.0,
                    help="if >0, abort before the run when analytical high-drive/floor crossing rate times dt exceeds this value")
    ap.add_argument("--no-run-summary", action="store_true",
                    help="disable writing run_summary.{json,txt,csv} at the end of the run")
'''
if '--preflight-rdt-warn' not in s:
    if old_args not in s:
        fail('Could not find pin diagnostic args block')
    s = s.replace(old_args, new_args, 1)

# -----------------------------------------------------------------------------
# 3. Insert preflight diagnostics before params and add to params.
# -----------------------------------------------------------------------------
preflight_marker = '''    params = vars(args).copy()
    params.update({
'''
preflight_code = r'''    # v14 preflight: cheap kinetic/time-scale diagnostics before the expensive run.
    kT_pre = max(KB_EV * T, 1.0e-300)
    if args.crossing_drive_mode == "force_work":
        cross_DG0_eV = cross.DeltaG_force_eV(0.0, T, cross_force_scale_N)
    else:
        cross_DG0_eV = cross.DeltaG_eV(0.0, T)
    cross_floor_DG_eV = max(0.0, cross.floor_frac * cross.G0_eV(T) - KB_EV * T * cross.entropy_kB)
    cross_R0_s = _rate_from_barrier(cross_DG0_eV, T, args.attempt_frequency_s)
    cross_Rfloor_s = _rate_from_barrier(cross_floor_DG_eV, T, args.attempt_frequency_s)
    cross_R0_dt = cross_R0_s * args.dt if math.isfinite(cross_R0_s) else float("nan")
    cross_Rfloor_dt = cross_Rfloor_s * args.dt if math.isfinite(cross_Rfloor_s) else float("nan")
    cross_tau0_s = 1.0 / cross_R0_s if cross_R0_s and math.isfinite(cross_R0_s) and cross_R0_s > 0 else float("inf")
    free_sweep_speed_red_s = args.glide_jump_length_reduced * required_free_net_rate_s
    sweep_time_spacing_s = forest_spacing_red / max(free_sweep_speed_red_s, 1.0e-300)
    sweep_time_over_tau0 = sweep_time_spacing_s / max(cross_tau0_s, 1.0e-300)
    pin_survival_prob_one_spacing = math.exp(-sweep_time_over_tau0) if sweep_time_over_tau0 < 700 else 0.0

    preflight_warnings = []
    if math.isfinite(cross_Rfloor_dt) and cross_Rfloor_dt > args.preflight_rdt_warn:
        preflight_warnings.append(
            f"high-drive/floor crossing R*dt={cross_Rfloor_dt:.3g} exceeds warn threshold {args.preflight_rdt_warn:g}; reduce dt or avoid cap/floor regime"
        )
    if math.isfinite(cross_R0_dt) and cross_R0_dt > args.preflight_rdt_warn:
        preflight_warnings.append(
            f"zero-force crossing R0*dt={cross_R0_dt:.3g} exceeds warn threshold {args.preflight_rdt_warn:g}; thermal depinning unresolved"
        )
    if sweep_time_over_tau0 > 3.0:
        preflight_warnings.append(
            f"pins likely thermally transparent before loading: sweep_time/tau0={sweep_time_over_tau0:.3g}"
        )
    elif sweep_time_over_tau0 < 0.1:
        preflight_warnings.append(
            f"pins likely survive many spacings before thermal escape: sweep_time/tau0={sweep_time_over_tau0:.3g}; jamming/cap risk"
        )

    preflight = {
        "cross_DG0_eV": cross_DG0_eV,
        "cross_floor_DG_eV": cross_floor_DG_eV,
        "cross_R0_s": cross_R0_s,
        "cross_Rfloor_s": cross_Rfloor_s,
        "cross_R0_dt": cross_R0_dt,
        "cross_Rfloor_dt": cross_Rfloor_dt,
        "cross_zero_force_lifetime_s": cross_tau0_s,
        "forest_spacing_reduced": forest_spacing_red,
        "free_sweep_speed_reduced_s": free_sweep_speed_red_s,
        "sweep_time_one_spacing_s": sweep_time_spacing_s,
        "sweep_time_over_zero_force_lifetime": sweep_time_over_tau0,
        "pin_survival_prob_one_spacing": pin_survival_prob_one_spacing,
        "warnings": preflight_warnings,
    }
    (outdir / "preflight_diagnostics.json").write_text(json.dumps(preflight, indent=2, sort_keys=True))
    (outdir / "preflight_diagnostics.txt").write_text("\n".join([
        "Preflight diagnostics",
        "=====================",
        f"cross_DG0_eV: {cross_DG0_eV}",
        f"cross_floor_DG_eV: {cross_floor_DG_eV}",
        f"cross_R0_s: {cross_R0_s}",
        f"cross_Rfloor_s: {cross_Rfloor_s}",
        f"cross_R0_dt: {cross_R0_dt}",
        f"cross_Rfloor_dt: {cross_Rfloor_dt}",
        f"cross_zero_force_lifetime_s: {cross_tau0_s}",
        f"forest_spacing_reduced: {forest_spacing_red}",
        f"free_sweep_speed_reduced_s: {free_sweep_speed_red_s}",
        f"sweep_time_one_spacing_s: {sweep_time_spacing_s}",
        f"sweep_time_over_zero_force_lifetime: {sweep_time_over_tau0}",
        f"pin_survival_prob_one_spacing: {pin_survival_prob_one_spacing}",
        "warnings:",
        *(f"  - {w}" for w in preflight_warnings),
    ]) + "\n")
    if args.preflight_rdt_stop and args.preflight_rdt_stop > 0.0 and math.isfinite(cross_Rfloor_dt) and cross_Rfloor_dt > args.preflight_rdt_stop:
        raise SystemExit(f"preflight stop: cross_Rfloor_dt={cross_Rfloor_dt:.6g} > {args.preflight_rdt_stop:.6g}")

'''
if 'preflight = {' not in s:
    if preflight_marker not in s:
        fail('Could not find params marker for preflight')
    s = s.replace(preflight_marker, preflight_code + preflight_marker, 1)

old_params_piece = '''        "v8_cross_effective_barrier_zero_stress_eV": cross.DeltaG_eV(0.0, T),
    })
'''
new_params_piece = '''        "v8_cross_effective_barrier_zero_stress_eV": cross.DeltaG_eV(0.0, T),
        "preflight_cross_DG0_eV": cross_DG0_eV,
        "preflight_cross_floor_DG_eV": cross_floor_DG_eV,
        "preflight_cross_R0_s": cross_R0_s,
        "preflight_cross_Rfloor_s": cross_Rfloor_s,
        "preflight_cross_R0_dt": cross_R0_dt,
        "preflight_cross_Rfloor_dt": cross_Rfloor_dt,
        "preflight_cross_zero_force_lifetime_s": cross_tau0_s,
        "preflight_sweep_time_one_spacing_s": sweep_time_spacing_s,
        "preflight_sweep_time_over_zero_force_lifetime": sweep_time_over_tau0,
        "preflight_pin_survival_prob_one_spacing": pin_survival_prob_one_spacing,
        "preflight_warnings": preflight_warnings,
    })
'''
if 'preflight_cross_Rfloor_dt' not in s:
    if old_params_piece not in s:
        fail('Could not find params tail insertion point')
    s = s.replace(old_params_piece, new_params_piece, 1)

# -----------------------------------------------------------------------------
# 4. Add history columns.
# -----------------------------------------------------------------------------
old_hist_start = '''        "step", "time_s", "eps_total", "eps_plastic", "sigma_MPa", "tau_MPa",
'''
new_hist_start = '''        "step", "time_s", "eps_total", "eps_plastic", "sigma_MPa", "tau_MPa",
        "tau_before_step_MPa", "tau_after_step_MPa",
'''
if 'tau_before_step_MPa' not in s:
    if old_hist_start not in s:
        fail('Could not find history start')
    s = s.replace(old_hist_start, new_hist_start, 1)

old_hist_deps = '''        "d_eps_p_book_free_glide", "d_eps_p_book_total", "d_eps_p_actual",
'''
new_hist_deps = '''        "d_eps_p_book_free_glide", "d_eps_p_book_total", "d_eps_p_actual",
        "d_eps_total", "d_tau_step_MPa", "d_eps_p_over_d_eps_total", "d_eps_p_per_depin_step",
'''
if 'd_tau_step_MPa' not in s:
    if old_hist_deps not in s:
        fail('Could not find history d_eps insertion')
    s = s.replace(old_hist_deps, new_hist_deps, 1)

# -----------------------------------------------------------------------------
# 5. Track actual d_eps_total step.
# -----------------------------------------------------------------------------
old_step_eps = '''            time_s += args.dt
            eps_total = min(args.target_strain, eps_total + args.strain_rate * args.dt)
'''
new_step_eps = '''            time_s += args.dt
            eps_total_prev = eps_total
            eps_total = min(args.target_strain, eps_total + args.strain_rate * args.dt)
            d_eps_total_step = eps_total - eps_total_prev
'''
if 'd_eps_total_step' not in s:
    if old_step_eps not in s:
        fail('Could not find eps_total update block')
    s = s.replace(old_step_eps, new_step_eps, 1)

# -----------------------------------------------------------------------------
# 6. Compute after-step stress and strain/depin diagnostics.
# -----------------------------------------------------------------------------
old_after_epsp = '''            eps_p += d_eps_p

            # Approximate line length diagnostic.
'''
new_after_epsp = '''            eps_p += d_eps_p
            tau_after_step_MPa = args.elastic_modulus_MPa * (eps_total - eps_p)
            d_tau_step_MPa = tau_after_step_MPa - tau_MPa
            d_eps_p_over_d_eps_total = d_eps_p / max(d_eps_total_step, 1.0e-300)
            d_eps_p_per_depin_step = d_eps_p / max(float(n_depin), 1.0) if n_depin > 0 else 0.0

            # Approximate line length diagnostic.
'''
if 'tau_after_step_MPa' not in s:
    if old_after_epsp not in s:
        fail('Could not find eps_p update block')
    s = s.replace(old_after_epsp, new_after_epsp, 1)

# -----------------------------------------------------------------------------
# 7. Fix history row stress and add fields.
# -----------------------------------------------------------------------------
old_row_start = '''                "step": step, "time_s": time_s,
                "eps_total": eps_total, "eps_plastic": eps_p,
                "sigma_MPa": tau_MPa, "tau_MPa": tau_MPa,
'''
new_row_start = '''                "step": step, "time_s": time_s,
                "eps_total": eps_total, "eps_plastic": eps_p,
                "sigma_MPa": tau_after_step_MPa, "tau_MPa": tau_after_step_MPa,
                "tau_before_step_MPa": tau_MPa,
                "tau_after_step_MPa": tau_after_step_MPa,
'''
if '"tau_before_step_MPa": tau_MPa' not in s:
    if old_row_start not in s:
        fail('Could not find row stress start')
    s = s.replace(old_row_start, new_row_start, 1)

old_row_deps = '''                "d_eps_p_book_free_glide": d_eps_p_book_free_glide,
                "d_eps_p_book_total": d_eps_p_book_total,
                "d_eps_p_actual": d_eps_p_actual,
'''
new_row_deps = '''                "d_eps_p_book_free_glide": d_eps_p_book_free_glide,
                "d_eps_p_book_total": d_eps_p_book_total,
                "d_eps_p_actual": d_eps_p_actual,
                "d_eps_total": d_eps_total_step,
                "d_tau_step_MPa": d_tau_step_MPa,
                "d_eps_p_over_d_eps_total": d_eps_p_over_d_eps_total,
                "d_eps_p_per_depin_step": d_eps_p_per_depin_step,
'''
if '"d_tau_step_MPa": d_tau_step_MPa' not in s:
    if old_row_deps not in s:
        fail('Could not find row d_eps insertion')
    s = s.replace(old_row_deps, new_row_deps, 1)

# -----------------------------------------------------------------------------
# 8. Write run summary before run.finished.
# -----------------------------------------------------------------------------
old_finish = '''
    (outdir / "run.finished").write_text("")
    return 0
'''
new_finish = '''
    if not args.no_run_summary:
        try:
            write_run_summary(outdir, params, args.target_strain)
        except Exception as exc:
            (outdir / "run_summary.error.txt").write_text(str(exc) + "\\n")

    (outdir / "run.finished").write_text("")
    return 0
'''
if 'write_run_summary(outdir, params, args.target_strain)' not in s:
    if old_finish not in s:
        fail('Could not find finish block')
    s = s.replace(old_finish, new_finish, 1)

# Model name update.
s = s.replace(
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v13_actual_swept_strain"',
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v14_diagnostics_summary"',
    1,
)

# Write and compile.
dst.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(dst)], check=True)

checks = [
    "preflight_diagnostics.json",
    "preflight_cross_Rfloor_dt",
    "tau_before_step_MPa",
    "tau_after_step_MPa",
    "d_tau_step_MPa",
    "write_run_summary",
    "run_summary.json",
    "clean_arrhenius_taylor_explicit_obstacles_v14_diagnostics_summary",
]
for c in checks:
    ok = c in s
    print(("OK   " if ok else "FAIL ") + c)
    if not ok:
        raise SystemExit("source audit failed")

print("PATCH COMPLETE")
print("Run with DRIVER=clean_arrhenius_taylor_ddd_v14.py")
