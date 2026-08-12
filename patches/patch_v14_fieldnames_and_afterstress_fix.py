#!/usr/bin/env python3
"""
patch_v14_fieldnames_and_afterstress_fix.py

Fixes v14 diagnostics patch issues:
  1) defines tau_after_step_MPa / d_tau_step_MPa / strain-ratio diagnostics
     after eps_p is updated;
  2) adds the new diagnostics to hist_cols so csv.DictWriter accepts them;
  3) compiles the patched driver.

Run from the OpenDiS directory:
    python3 patch_v14_fieldnames_and_afterstress_fix.py
    python3 -m py_compile clean_arrhenius_taylor_ddd_v14.py
"""
from pathlib import Path
import shutil
import subprocess
import sys
import time

src = Path("clean_arrhenius_taylor_ddd_v14.py")
if not src.exists():
    raise SystemExit("Could not find clean_arrhenius_taylor_ddd_v14.py")

s = src.read_text()
backup = src.with_suffix(src.suffix + f".bak_before_v14_fieldnames_fix_{int(time.time())}")
shutil.copy2(src, backup)

# 1) Define after-step stress diagnostics if missing.
old = """            eps_p += d_eps_p

            # Approximate line length diagnostic.
"""
new = """            eps_p += d_eps_p

            # v14: after-step stress diagnostics.  tau_MPa above is the
            # beginning-of-step elastic stress used for the mobility solve.
            # History should report the after-step stress consistent with the
            # updated plastic strain, and retain the before-step value.
            tau_before_step_MPa = tau_MPa
            tau_after_step_MPa = args.elastic_modulus_MPa * (eps_total - eps_p)
            d_tau_step_MPa = tau_after_step_MPa - tau_before_step_MPa
            d_eps_p_over_d_eps_total = (
                d_eps_p / d_eps_total_step if abs(d_eps_total_step) > 1.0e-300 else 0.0
            )
            d_eps_p_per_depin_step = (
                d_eps_p / float(n_depin) if n_depin > 0 else 0.0
            )

            # Approximate line length diagnostic.
"""
if "tau_after_step_MPa = args.elastic_modulus_MPa * (eps_total - eps_p)" not in s:
    if old not in s:
        raise SystemExit("Could not find insertion point after eps_p += d_eps_p")
    s = s.replace(old, new, 1)

# 2) Ensure tau_before row uses explicit before-step variable, not whatever tau_MPa is later.
s = s.replace('"tau_before_step_MPa": tau_MPa,', '"tau_before_step_MPa": tau_before_step_MPa,')

# 3) Add row fields to hist_cols if they are not already listed.
old_cols = """        "d_eps_p", "d_eps_p_swept_area",
        "d_eps_p_book_free_glide", "d_eps_p_book_total", "d_eps_p_actual",
"""
new_cols = """        "d_eps_p", "d_eps_p_swept_area",
        "d_eps_total", "d_tau_step_MPa", "d_eps_p_over_d_eps_total", "d_eps_p_per_depin_step",
        "d_eps_p_book_free_glide", "d_eps_p_book_total", "d_eps_p_actual",
"""
if '"d_tau_step_MPa"' not in s[s.find('hist_cols = ['):s.find('event_cols = [')]:
    if old_cols not in s:
        raise SystemExit("Could not find hist_cols insertion point for v14 diagnostics")
    s = s.replace(old_cols, new_cols, 1)

src.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(src)], check=True)

print(f"backup: {backup}")
print("PATCH COMPLETE: v14 fieldnames and after-step stress diagnostics fixed")
print("Re-run: ./run_v14_sanity_onecase.sh")
