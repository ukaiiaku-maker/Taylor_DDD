#!/usr/bin/env python3
from pathlib import Path
import sys

p = Path("arrhenius_single_glider_fixed_forest_test.py")
s = p.read_text()

checks = {
    "has expfit Peierls barrier function": "arrhenius_unpinned_glide_barrier_ev" in s,
    "has expfit Peierls velocity function": "arrhenius_unpinned_glide_velocity_reduced_s" in s,
    "main glide update calls expfit velocity": "vg = arrhenius_unpinned_glide_velocity_reduced_s(" in s,
    "main glide update converts MPa to Pa": "tau_mpa * 1.0e6" in s,
    "legacy glide retained only as fallback": "else:\n                    vg = net_glide_speed_reduced_per_s(" in s,
}

bad = False
print("PEIERLS EXPFIT CALLSITE AUDIT")
for name, ok in checks.items():
    print(("OK   " if ok else "FAIL ") + name)
    bad |= not ok

if bad:
    print("\nDo not run. Peierls expfit is not wired into the actual glide update.")
    sys.exit(2)

print("\nAUDIT PASSED")
