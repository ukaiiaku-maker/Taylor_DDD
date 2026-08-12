#!/usr/bin/env python3
"""
patch_v11_to_v12_pinned_line_backstress.py

Create clean_arrhenius_taylor_ddd_v12.py from clean_arrhenius_taylor_ddd_v11.py.

Core change:
  The v10/v11 external_drive projection is kept only for completely unpinned
  periodic lines.  For lines with live pins, raw local tau_app - tau_back
  mobility is preserved so curvature/line-length work can reduce net swept
  area and produce Taylor hardening.

This is not a segment-window/artificial gating patch.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import time

src = Path("clean_arrhenius_taylor_ddd_v11.py")
if not src.exists():
    alt = Path("clean_arrhenius_taylor_ddd_v11(1).py")
    if alt.exists():
        src = alt
    else:
        raise SystemExit("Could not find clean_arrhenius_taylor_ddd_v11.py")

dst = Path("clean_arrhenius_taylor_ddd_v12.py")
s = src.read_text()
backup = src.with_suffix(src.suffix + f".bak_before_v12_{int(time.time())}")
shutil.copy2(src, backup)
print(f"backup: {backup}")
print(f"source: {src}")
print(f"dest:   {dst}")

def fail(msg):
    raise SystemExit("ERROR: " + msg)

old_arg = """    ap.add_argument("--backstress-com-projection", choices=["external_drive", "none"], default="external_drive",
                    help="v10: with backstress mobility on, project the internal line-tension \\
                         self-force so it cannot create a spurious net center-of-mass friction \\
                         on a free periodic line. external_drive preserves the line-mean \\
                         glide set by the applied Peierls stress while retaining curvature-driven \\
                         relative node motion. none reproduces the raw v9 local tau_app-tau_back law.")
"""
new_arg = old_arg + """    ap.add_argument("--project-backstress-on-pinned-lines", action="store_true",
                    help="v12 diagnostic/physics switch. By default the external-drive COM projection "
                         "is applied only to lines with no live pins. A pinned line keeps the raw "
                         "local tau_app - tau_back mobility, so curvature/line-length work can "
                         "reduce net swept area. Set this flag to recover old v10/v11 behavior.")
"""
if "--project-backstress-on-pinned-lines" not in s:
    if old_arg not in s:
        fail("Could not find backstress-com-projection argument block.")
    s = s.replace(old_arg, new_arg, 1)

old_proj = """            dx_all = dx_all_raw.copy()
            if args.backstress_mobility == "on" and args.backstress_com_projection == "external_drive":
                G_app_fwd = peierls.DeltaG_eV(max(tau_MPa, 0.0), T)
                G_app_rev = peierls.DeltaG_eV(max(-tau_MPa, 0.0), T)
                net_rate_app = pref_glide * (math.exp(-G_app_fwd / kT) - math.exp(-G_app_rev / kT))
                dx_app = float(np.clip(args.glide_jump_length_reduced * net_rate_app * args.dt,
                                       -args.max_free_dx_reduced, args.max_free_dx_reduced))
                for li_proj in range(nline):
                    fmask = free[li_proj]
                    if np.any(fmask):
                        mean_raw = float(np.mean(dx_all_raw[li_proj, fmask]))
                        dx_all[li_proj, fmask] = dx_all_raw[li_proj, fmask] - mean_raw + dx_app
                dx_all = np.clip(dx_all, -args.max_free_dx_reduced, args.max_free_dx_reduced)
"""
new_proj = """            dx_all = dx_all_raw.copy()
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
"""
if "n_pinned_raw_lines" not in s:
    if old_proj not in s:
        fail("Could not find v10 external_drive projection block.")
    s = s.replace(old_proj, new_proj, 1)

old_params = """        "v10_backstress_com_projection": args.backstress_com_projection,
        "v4_line_tension_smooth_frac": args.line_tension_smooth_frac,
"""
new_params = """        "v10_backstress_com_projection": args.backstress_com_projection,
        "v12_project_backstress_on_pinned_lines": args.project_backstress_on_pinned_lines,
        "v12_backstress_projection_rule": "project_only_unpinned_lines_unless_flag_set",
        "v4_line_tension_smooth_frac": args.line_tension_smooth_frac,
"""
if "v12_backstress_projection_rule" not in s:
    if old_params not in s:
        fail("Could not find params projection block.")
    s = s.replace(old_params, new_params, 1)

old_hist = """        "n_free_nodes", "free_node_fraction", "pinned_node_fraction",
        "blocked_pair_count", "tau_back_abs_mean_MPa", "tau_back_abs_p90_MPa",
"""
new_hist = """        "n_free_nodes", "free_node_fraction", "pinned_node_fraction",
        "mean_dx_raw_free_reduced", "mean_dx_after_projection_free_reduced", "dx_app_reduced",
        "n_com_projected_lines", "n_pinned_raw_lines",
        "blocked_pair_count", "tau_back_abs_mean_MPa", "tau_back_abs_p90_MPa",
"""
if "mean_dx_raw_free_reduced" not in s:
    if old_hist not in s:
        fail("Could not find history column insertion point.")
    s = s.replace(old_hist, new_hist, 1)

old_row = """                "n_free_nodes": n_free_nodes,
                "free_node_fraction": n_free_nodes / max(n_total_nodes, 1),
                "pinned_node_fraction": 1.0 - n_free_nodes / max(n_total_nodes, 1),
                "blocked_pair_count": len(blocked_until_far),
"""
new_row = """                "n_free_nodes": n_free_nodes,
                "free_node_fraction": n_free_nodes / max(n_total_nodes, 1),
                "pinned_node_fraction": 1.0 - n_free_nodes / max(n_total_nodes, 1),
                "mean_dx_raw_free_reduced": mean_dx_raw_free,
                "mean_dx_after_projection_free_reduced": mean_dx_after_projection_free,
                "dx_app_reduced": dx_app,
                "n_com_projected_lines": n_com_projected_lines,
                "n_pinned_raw_lines": n_pinned_raw_lines,
                "blocked_pair_count": len(blocked_until_far),
"""
if '"n_pinned_raw_lines": n_pinned_raw_lines' not in s:
    if old_row not in s:
        fail("Could not find history row insertion point.")
    s = s.replace(old_row, new_row, 1)

s = s.replace(
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v11_force_work"',
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v12_pinned_line_raw_backstress"',
    1,
)

dst.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(dst)], check=True)

checks = [
    "--project-backstress-on-pinned-lines",
    "n_pinned_raw_lines",
    "mean_dx_raw_free_reduced",
    "v12_backstress_projection_rule",
]
for c in checks:
    print(("OK   " if c in s else "FAIL ") + c)
    if c not in s:
        raise SystemExit("source audit failed")

print("PATCH COMPLETE")
print("Run with DRIVER=clean_arrhenius_taylor_ddd_v12.py")
