#!/usr/bin/env python3
"""Patch clean_arrhenius_taylor_ddd_v14.py to v15 branch-specific crossing floor diagnostics.

v15 adds:
  --expfit-cross-floor-frac and --expfit-peierls-floor-frac
  branch-specific EXP floors, so high-drive crossing-rate resolution can be
  controlled without changing the zero-force barrier or Peierls branch.
  run_summary now includes preflight values correctly.
"""
from pathlib import Path
import shutil, subprocess, sys, time

src = Path("clean_arrhenius_taylor_ddd_v14.py")
if not src.exists():
    raise SystemExit("Could not find clean_arrhenius_taylor_ddd_v14.py")
backup = src.with_suffix(src.suffix + f".bak_before_v15_branch_floor_{int(time.time())}")
shutil.copy2(src, backup)
print(f"backup: {backup}")
s = src.read_text()

def must_replace(old, new, label):
    global s
    if new in s:
        print(f"already patched: {label}")
        return
    if old not in s:
        raise SystemExit(f"ERROR: could not find block for {label}")
    s = s.replace(old, new, 1)
    print(f"patched: {label}")

must_replace(
'    ap.add_argument("--expfit-floor-frac", type=float, default=0.0)\n    ap.add_argument("--expfit-cross-scale", type=float, default=0.40)',
'    ap.add_argument("--expfit-floor-frac", type=float, default=0.0)\n    ap.add_argument("--expfit-cross-floor-frac", type=float, default=-1.0,\n                    help="v15: forest/crossing EXP-floor fraction. >=0 overrides --expfit-floor-frac for crossing only. Use this to bound high-drive crossing rates without changing Peierls or the zero-force barrier.")\n    ap.add_argument("--expfit-peierls-floor-frac", type=float, default=-1.0,\n                    help="v15: Peierls EXP-floor fraction. >=0 overrides --expfit-floor-frac for Peierls only.")\n    ap.add_argument("--expfit-cross-scale", type=float, default=0.40)',
"branch floor args")

must_replace(
'def make_barrier(args, scale: float, entropy: float, sigc0_MPa: float) -> ExpFitBarrier:\n    return ExpFitBarrier(\n        scale=scale,\n        entropy_kB=entropy,\n        floor_frac=args.expfit_floor_frac,\n        T0_K=args.expfit_T0_K,',
'def make_barrier(args, scale: float, entropy: float, sigc0_MPa: float, floor_frac: float | None = None) -> ExpFitBarrier:\n    floor = args.expfit_floor_frac if floor_frac is None else float(floor_frac)\n    return ExpFitBarrier(\n        scale=scale,\n        entropy_kB=entropy,\n        floor_frac=floor,\n        T0_K=args.expfit_T0_K,',
"make_barrier floor parameter")

must_replace(
'    cross_sigc0 = args.expfit_cross_sigc0_MPa if args.expfit_cross_sigc0_MPa > 0 else args.expfit_sigc0_MPa\n    peierls_sigc0 = args.expfit_peierls_sigc0_MPa if args.expfit_peierls_sigc0_MPa > 0 else args.expfit_sigc0_MPa\n    cross = make_barrier(args, args.expfit_cross_scale, args.expfit_cross_entropy_kB, cross_sigc0)\n    peierls = make_barrier(args, args.expfit_peierls_scale, args.expfit_peierls_entropy_kB, peierls_sigc0)',
'    cross_sigc0 = args.expfit_cross_sigc0_MPa if args.expfit_cross_sigc0_MPa > 0 else args.expfit_sigc0_MPa\n    peierls_sigc0 = args.expfit_peierls_sigc0_MPa if args.expfit_peierls_sigc0_MPa > 0 else args.expfit_sigc0_MPa\n    # v15: branch-specific EXP-floor fractions.  The floor fraction changes the\n    # high-drive minimum barrier but leaves the zero-force barrier H(0)=G0 unchanged.\n    # This is the clean way to bound max crossing rates without retuning pin survival.\n    cross_floor_frac = args.expfit_cross_floor_frac if args.expfit_cross_floor_frac >= 0.0 else args.expfit_floor_frac\n    peierls_floor_frac = args.expfit_peierls_floor_frac if args.expfit_peierls_floor_frac >= 0.0 else args.expfit_floor_frac\n    cross = make_barrier(args, args.expfit_cross_scale, args.expfit_cross_entropy_kB, cross_sigc0, cross_floor_frac)\n    peierls = make_barrier(args, args.expfit_peierls_scale, args.expfit_peierls_entropy_kB, peierls_sigc0, peierls_floor_frac)',
"branch floor barrier creation")

s = s.replace('"model": "clean_arrhenius_taylor_explicit_obstacles_v14_diagnostics_summary"', '"model": "clean_arrhenius_taylor_explicit_obstacles_v15_branch_floor_preflight"', 1)

must_replace(
'        "v11_cross_work_scale_eV": cross_work_scale_eV,\n        "v3_plastic_strain_source": args.plastic_strain_source,',
'        "v11_cross_work_scale_eV": cross_work_scale_eV,\n        "v15_cross_floor_frac_used": cross_floor_frac,\n        "v15_peierls_floor_frac_used": peierls_floor_frac,\n        "v15_cross_floor_rate_is_upper_bound": True,\n        "preflight_cross_R0_dt": cross_R0_dt,\n        "preflight_cross_Rfloor_dt": cross_Rfloor_dt,\n        "preflight_sweep_time_over_zero_force_lifetime": sweep_time_over_tau0,\n        "preflight_pin_survival_prob_one_spacing": pin_survival_prob_one_spacing,\n        "preflight_cross_R0_s": cross_R0_s,\n        "preflight_cross_Rfloor_s": cross_Rfloor_s,\n        "preflight_cross_DG0_eV": cross_DG0_eV,\n        "preflight_cross_floor_DG_eV": cross_floor_DG_eV,\n        "v3_plastic_strain_source": args.plastic_strain_source,',
"params preflight/floor metadata")

must_replace(
'        "max_crossing_rate_dt", "preflight_cross_Rfloor_dt",\n        "preflight_sweep_time_over_zero_force_lifetime",',
'        "max_crossing_rate_dt", "preflight_cross_Rfloor_dt",\n        "preflight_sweep_time_over_zero_force_lifetime",\n        "preflight_pin_survival_prob_one_spacing",\n        "v15_cross_floor_frac_used",',
"run summary displayed preflight/floor fields")

src.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(src)], check=True)
print("PATCH COMPLETE: v15 branch-specific floor/preflight diagnostics")
print("Run with DRIVER=clean_arrhenius_taylor_ddd_v14.py (now patched to v15), or copy to clean_arrhenius_taylor_ddd_v15.py")
