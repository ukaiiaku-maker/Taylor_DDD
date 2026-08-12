#!/usr/bin/env python3
"""
patch_v12_to_v13_actual_swept_strain.py

Create clean_arrhenius_taylor_ddd_v13.py from clean_arrhenius_taylor_ddd_v12.py.

Core change:
  Plastic strain/stress feedback is computed from the actual step-wise swept
  line motion after glide, relaxation, capture snapping, and depin state changes,
  not from the pre-capture/pre-depin dx_free bookkeeping array.

New behavior:
  --plastic-strain-source actual  (default)
      d_eps_p = swept area from dx_actual = min_image(x_lines_end - x_lines_start)
      using the final coordinates after all step updates.
  --plastic-strain-source free_glide
      old dx_free feedback retained, but actual swept strain is still logged.
  --plastic-strain-source total
      old dx_total feedback retained, but actual swept strain is still logged.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import time

src = Path("clean_arrhenius_taylor_ddd_v12.py")
if not src.exists():
    alt = Path("clean_arrhenius_taylor_ddd_v12(1).py")
    if alt.exists():
        src = alt
    else:
        raise SystemExit("Could not find clean_arrhenius_taylor_ddd_v12.py")

dst = Path("clean_arrhenius_taylor_ddd_v13.py")
s = src.read_text()
backup = src.with_suffix(src.suffix + f".bak_before_v13_actual_swept_{int(time.time())}")
shutil.copy2(src, backup)
print(f"backup: {backup}")
print(f"source: {src}")
print(f"dest:   {dst}")

def replace_once(old: str, new: str, label: str):
    global s
    if old not in s:
        raise SystemExit(f"ERROR: could not find block for {label}")
    s = s.replace(old, new, 1)
    print(f"OK   {label}")

old_arg = """    ap.add_argument("--plastic-strain-source", choices=["free_glide", "total"], default="free_glide",
                    help="free_glide: only driven forward glide sweeps plastic area "
                         "(line-tension relaxation does not contaminate eps_p).")
"""
new_arg = """    ap.add_argument("--plastic-strain-source", choices=["actual", "free_glide", "total"], default="actual",
                    help="actual: use measured step-wise line-coordinate change after glide, "
                         "relaxation, capture/depin, and capture snapping for plastic-strain "
                         "feedback. free_glide and total retain old bookkeeping modes as "
                         "diagnostics.")
"""
if 'choices=["actual", "free_glide", "total"]' not in s:
    replace_once(old_arg, new_arg, "plastic-strain-source actual/default")

old_params = """        "v3_plastic_strain_source": args.plastic_strain_source,
        "v5_signed_mobility": True,
"""
new_params = """        "v3_plastic_strain_source": args.plastic_strain_source,
        "v13_actual_swept_strain_feedback": (args.plastic_strain_source == "actual"),
        "v13_strain_feedback_rule": "actual uses min_image(x_end - x_start) after all step updates",
        "v5_signed_mobility": True,
"""
if '"v13_actual_swept_strain_feedback"' not in s:
    replace_once(old_params, new_params, "params v13 actual swept strain metadata")

old_hist = """        "forest_rho_actual_m2", "n_obstacles_active", "n_pinned_nodes", "n_crossed_total",
        "d_eps_p", "d_eps_p_swept_area", "line_length_reduced",
        "mean_dx_free_reduced", "mean_dx_relax_reduced", "max_dx_total_reduced",
"""
new_hist = """        "forest_rho_actual_m2", "n_obstacles_active", "n_pinned_nodes", "n_crossed_total",
        "d_eps_p", "d_eps_p_swept_area",
        "d_eps_p_book_free_glide", "d_eps_p_book_total", "d_eps_p_actual",
        "eps_plastic_book_free_glide", "eps_plastic_book_total", "eps_plastic_actual",
        "line_length_reduced",
        "mean_dx_free_reduced", "mean_dx_relax_reduced", "mean_dx_actual_reduced",
        "max_dx_total_reduced", "max_abs_dx_actual_reduced",
"""
if '"d_eps_p_book_free_glide"' not in s:
    replace_once(old_hist, new_hist, "history columns for actual/book strain diagnostics")

old_state = """    eps_total = 0.0
    eps_p = 0.0
    time_s = 0.0
    n_crossed_total = 0
    pin_rows_written = 0
"""
new_state = """    eps_total = 0.0
    # eps_p is the strain variable used for stress feedback.  In v13 the default
    # is the actual swept-area strain; old bookkeeping variants remain available
    # through --plastic-strain-source for A/B diagnostics.
    eps_p = 0.0
    eps_p_book_free_glide = 0.0
    eps_p_book_total = 0.0
    eps_p_actual = 0.0
    # Track unwrapped x positions for final-state diagnostics.  Stepwise
    # min-image increments are safe because max_free_dx and the relaxation
    # regularizer keep per-step motion far below Lx/2.
    x_unwrapped_lines = x_lines.copy()
    time_s = 0.0
    n_crossed_total = 0
    pin_rows_written = 0
"""
if 'eps_p_book_free_glide = 0.0' not in s:
    replace_once(old_state, new_state, "strain accumulators and unwrapped positions")

old_plastic = """            # v3: plastic strain from the DRIVEN forward glide only by default, so the
            # line-tension relaxation (shape response) does not leak into eps_p.
            dx_plastic = dx_free if args.plastic_strain_source == "free_glide" else dx_total
            d_eps_p = swept_strain_increment(dx_plastic, Lx, Lz, b, s_out)
            eps_p += d_eps_p
"""
new_plastic = """            # v13: plastic strain from actual swept line area by default.
            #
            # Old v12 used dx_free by default, which is a pre-capture/pre-depin
            # bookkeeping array.  That can decouple n_depin from eps_p when the
            # high-density line advances through many capture/depin events, because
            # the selected dx array does not necessarily equal the final coordinate
            # change of the line after relaxation, capture snapping, and state updates.
            #
            # Keep the old book values as diagnostics, but use the actual final-step
            # displacement for stress feedback when --plastic-strain-source=actual.
            dx_book_free_glide = dx_free
            dx_book_total = dx_total
            dx_actual = minimum_image_delta(x_lines - x_before, Lx)

            d_eps_p_book_free_glide = swept_strain_increment(dx_book_free_glide, Lx, Lz, b, s_out)
            d_eps_p_book_total = swept_strain_increment(dx_book_total, Lx, Lz, b, s_out)
            d_eps_p_actual = swept_strain_increment(dx_actual, Lx, Lz, b, s_out)

            eps_p_book_free_glide += d_eps_p_book_free_glide
            eps_p_book_total += d_eps_p_book_total
            eps_p_actual += d_eps_p_actual
            x_unwrapped_lines += dx_actual

            if args.plastic_strain_source == "free_glide":
                d_eps_p = d_eps_p_book_free_glide
            elif args.plastic_strain_source == "total":
                d_eps_p = d_eps_p_book_total
            else:
                d_eps_p = d_eps_p_actual
            eps_p += d_eps_p
"""
if 'd_eps_p_book_free_glide = swept_strain_increment' not in s:
    replace_once(old_plastic, new_plastic, "actual swept-area strain feedback block")

old_row = """                "n_crossed_total": n_crossed_total,
                "d_eps_p": d_eps_p,
                "d_eps_p_swept_area": d_eps_p,
                "line_length_reduced": line_len,
                "mean_dx_free_reduced": float(np.mean(dx_free)),
                "mean_dx_relax_reduced": float(np.mean(dx_relax)),
                "max_dx_total_reduced": float(np.max(dx_total)),
"""
new_row = """                "n_crossed_total": n_crossed_total,
                "d_eps_p": d_eps_p,
                "d_eps_p_swept_area": d_eps_p_actual,
                "d_eps_p_book_free_glide": d_eps_p_book_free_glide,
                "d_eps_p_book_total": d_eps_p_book_total,
                "d_eps_p_actual": d_eps_p_actual,
                "eps_plastic_book_free_glide": eps_p_book_free_glide,
                "eps_plastic_book_total": eps_p_book_total,
                "eps_plastic_actual": eps_p_actual,
                "line_length_reduced": line_len,
                "mean_dx_free_reduced": float(np.mean(dx_free)),
                "mean_dx_relax_reduced": float(np.mean(dx_relax)),
                "mean_dx_actual_reduced": float(np.mean(dx_actual)),
                "max_dx_total_reduced": float(np.max(dx_total)),
                "max_abs_dx_actual_reduced": float(np.max(np.abs(dx_actual))),
"""
if '"d_eps_p_book_free_glide": d_eps_p_book_free_glide' not in s:
    replace_once(old_row, new_row, "history row actual/book strain diagnostics")

old_final = """    with open(outdir / "single_glider_final_nodes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["line_id", "node_id", "x_reduced", "z_reduced", "pinned", "pinned_obs"])
        for li in range(nline):
            for j in range(nn):
                w.writerow([li, j, x_lines[li, j], z_lines[li, j], int(pinned[li, j]), int(pinned_obs[li, j])])
"""
new_final = """    with open(outdir / "single_glider_final_nodes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["line_id", "node_id", "x_reduced", "x_unwrapped_reduced", "z_reduced", "pinned", "pinned_obs"])
        for li in range(nline):
            for j in range(nn):
                w.writerow([
                    li, j,
                    x_lines[li, j],
                    x_unwrapped_lines[li, j],
                    z_lines[li, j],
                    int(pinned[li, j]),
                    int(pinned_obs[li, j]),
                ])
"""
if '"x_unwrapped_reduced"' not in s:
    replace_once(old_final, new_final, "final node unwrapped x output")

s = s.replace(
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v12_pinned_line_raw_backstress"',
    '"model": "clean_arrhenius_taylor_explicit_obstacles_v13_actual_swept_strain"',
    1,
)

dst.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(dst)], check=True)

checks = [
    'choices=["actual", "free_glide", "total"]',
    "v13_actual_swept_strain_feedback",
    "d_eps_p_book_free_glide",
    "d_eps_p_actual",
    "x_unwrapped_reduced",
]
for c in checks:
    print(("OK   " if c in s else "FAIL ") + c)
    if c not in s:
        raise SystemExit("source audit failed")

print("PATCH COMPLETE")
print("Run with DRIVER=clean_arrhenius_taylor_ddd_v13.py")
