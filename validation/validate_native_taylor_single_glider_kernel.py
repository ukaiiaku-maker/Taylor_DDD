#!/usr/bin/env python3
"""Regression of the native Taylor kernel against published v17 force work."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


KB_EV = 8.617333262145e-5


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    import pyexadis

    temperature = 1100.0
    b_m = 2.48e-10
    mu_pa = 80.0e9
    alpha = 0.5
    force_scale_factor = 0.25
    fc = force_scale_factor * alpha * mu_pa * b_m**2
    sigma_c_gpa = 1.497042242375928
    vstar_b3 = fc * b_m / (sigma_c_gpa * 1.0e9 * b_m**3)
    H = 0.40 * 1.908192
    entropy = -9.25
    floor = 0.50
    shape_a = 2.2056211004282904
    shape_n = 2.5207319790155385
    eta0 = 1.0e12
    dt = 1.0e-8

    cases = []
    max_barrier_error = 0.0
    max_rate_relative_error = 0.0
    spacing_invariance = True
    for ratio in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        force = ratio * fc
        reference_H = floor * H + (1.0 - floor) * H * math.exp(
            -shape_a * ratio**shape_n
        )
        reference_G = max(0.0, reference_H - KB_EV * temperature * entropy)
        reference_rate = eta0 * math.exp(-reference_G / (KB_EV * temperature))
        spacing_results = []
        for spacing_b in (100.0, 1000.0, 10000.0):
            native = dict(pyexadis.evaluate_taylor_line_tension_interaction(
                force, spacing_b * b_m, spacing_b * b_m,
                100.0e6, dt, temperature, H, entropy, sigma_c_gpa,
                floor, shape_a, shape_n, eta0, vstar_b3, 1.0,
                b_m, mu_pa, alpha,
            ))
            max_barrier_error = max(max_barrier_error, abs(native["G_used_eV"] - reference_G))
            max_rate_relative_error = max(
                max_rate_relative_error,
                abs(native["rate_s"] - reference_rate) / max(reference_rate, 1.0e-300),
            )
            spacing_results.append(native)
        spacing_invariance &= max(r["G_used_eV"] for r in spacing_results) == min(
            r["G_used_eV"] for r in spacing_results
        )
        cases.append({
            "force_ratio_F_over_Fc": ratio,
            "reference_G_used_eV": reference_G,
            "native_G_used_eV": spacing_results[0]["G_used_eV"],
            "native_tau_eff_Pa": spacing_results[0]["tau_eff_Pa"],
            "native_phi_geom_by_spacing": [r["phi_geom_L_over_b"] for r in spacing_results],
        })

    passed = (
        max_barrier_error <= 1.0e-12 and
        max_rate_relative_error <= 1.0e-12 and
        spacing_invariance
    )
    result = {
        "status": "passed" if passed else "failed",
        "reference": "main:clean_arrhenius_taylor_ddd_v17.py force_work mode",
        "force_scale_N": fc,
        "vstar_b3_for_exact_force_scale_mapping": vstar_b3,
        "max_barrier_absolute_error_eV": max_barrier_error,
        "max_rate_relative_error": max_rate_relative_error,
        "barrier_invariant_to_L_eff_at_fixed_force": spacing_invariance,
        "L_eff_over_b_is_diagnostic_only": spacing_invariance,
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
