#!/usr/bin/env python3
"""
patch_v14_fix_after_step_stress.py

Fix clean_arrhenius_taylor_ddd_v14.py NameError:
    tau_after_step_MPa is not defined

The original v14 patch changed the history row to report after-step stress but
forgot to compute tau_after_step_MPa, d_tau_step_MPa, and related diagnostics
after eps_p is updated.  This patch inserts those definitions immediately after
the plastic strain increment and updates the history row to use the explicit
before-step stress variable.
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
backup = src.with_suffix(src.suffix + f".bak_before_v14_afterstress_fix_{int(time.time())}")
shutil.copy2(src, backup)

old = """            eps_p += d_eps_p

            # Approximate line length diagnostic.
"""
new = """            eps_p += d_eps_p

            # v14: after-step stress diagnostics.  The trial/mobility stress tau_MPa
            # above is the beginning-of-step elastic stress, before this step's plastic
            # swept-area increment.  For a consistent history row, report the stress
            # after eps_p has been updated, and retain the before-step value explicitly.
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

s = s.replace('"tau_before_step_MPa": tau_MPa,', '"tau_before_step_MPa": tau_before_step_MPa,', 1)

src.write_text(s)
subprocess.run([sys.executable, "-m", "py_compile", str(src)], check=True)

print(f"backup: {backup}")
print("PATCH COMPLETE: v14 after-step stress diagnostics now defined")
print("Re-run: ./run_v14_sanity_onecase.sh")
