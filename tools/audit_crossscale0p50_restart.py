#!/usr/bin/env python3
from pathlib import Path
import csv
import math

ROOT = Path("results/multiglider_crossscale0p50_fixedPeierls0p1_N4_T700_1100_density11_rate5e2_60k")
STEPS_EXPECTED = 60000

temps = [700, 800, 900, 1000, 1100]
rhos = [1e12, 5e12, 1e13, 5e13, 1e14, 5e14, 1e15, 5e15, 1e16, 5e16, 1e17]
seed = 11

def rtag(rho):
    return f"{rho:.1e}".replace("+", "").replace(".", "p")

def final_step(hist):
    if not hist.exists() or hist.stat().st_size == 0:
        return None
    last = None
    try:
        with hist.open() as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                last = row
        if last is None:
            return None
        return int(float(last.get("step", "nan")))
    except Exception:
        return None

rows = []
missing_or_incomplete = []

for T in temps:
    for rho in rhos:
        name = f"T{T}_rho{rtag(rho)}_N4_seed{seed}"
        d = ROOT / name
        hist = d / "single_glider_history.csv"
        step = final_step(hist)

        if step is None:
            status = "MISSING"
            missing_or_incomplete.append(d)
        elif step >= STEPS_EXPECTED:
            status = "COMPLETE"
        else:
            status = f"INCOMPLETE_step{step}"
            missing_or_incomplete.append(d)

        rows.append((T, rho, name, status, step))

print(f"ROOT = {ROOT}")
print(f"Expected cases = {len(rows)}")
print()
for T, rho, name, status, step in rows:
    print(f"{status:22s} T={T:<4g} rho={rho:.1e}  {name}")

print()
print(f"Missing/incomplete cases = {len(missing_or_incomplete)}")
for d in missing_or_incomplete:
    print(d)

(Path("restart_missing_or_incomplete_cases.txt")).write_text(
    "\n".join(str(d) for d in missing_or_incomplete) + "\n"
)
print("\nwrote: restart_missing_or_incomplete_cases.txt")
