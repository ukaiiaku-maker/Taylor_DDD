#!/usr/bin/env python3
from pathlib import Path
import re
import sys

p = Path("arrhenius_single_glider_fixed_forest_test.py")
s = p.read_text()

print("PEIERLS / GLIDE SOURCE AUDIT")
print("="*72)

patterns = [
    "expfit_peierls_scale",
    "glide_barrier_frac_of_nuc",
    "glide_activation_volume_b3",
    "mobile_glide_prefactor",
    "Peierls",
    "peierls",
    "glide",
]

for pat in patterns:
    hits = [(m.start(), m.group(0)) for m in re.finditer(re.escape(pat), s)]
    print(f"{pat}: {len(hits)} hit(s)")

print("\nRelevant source snippets:")
for key in ["expfit_peierls_scale", "glide_barrier_frac_of_nuc", "glide_activation_volume_b3"]:
    for m in re.finditer(re.escape(key), s):
        a = max(0, m.start() - 900)
        b = min(len(s), m.end() + 900)
        print("\n" + "-"*72)
        print(f"SNIPPET AROUND: {key}")
        print("-"*72)
        print(s[a:b])

# Hard fail conditions.
failures = []

if "expfit_peierls_scale" not in s:
    failures.append("driver does not contain expfit_peierls_scale")

if "glide_barrier_frac_of_nuc" in s:
    print("\nWARNING: driver still contains glide_barrier_frac_of_nuc.")
    print("This may be okay only if it is not used when barrier-family=expfit_floor.")

# Check whether expfit_peierls_scale appears inside more than parser/help/run-command writing.
# This is heuristic but useful.
uses = []
for m in re.finditer("expfit_peierls_scale", s):
    context = s[max(0,m.start()-250):min(len(s),m.end()+250)]
    if "add_argument" not in context and "run_command" not in context:
        uses.append(context)

print(f"\nNon-parser uses of expfit_peierls_scale: {len(uses)}")
for i, u in enumerate(uses[:5], 1):
    print("\n" + "-"*72)
    print(f"USE {i}")
    print("-"*72)
    print(u)

if len(uses) == 0:
    failures.append("expfit_peierls_scale appears only in parser/command plumbing; Peierls may not use it")

if failures:
    print("\nAUDIT FAILED:")
    for f in failures:
        print(" -", f)
    sys.exit(2)

print("\nAUDIT PASSED at heuristic level.")
print("Next: run Peierls scale sensitivity smoke test to verify behavior.")
