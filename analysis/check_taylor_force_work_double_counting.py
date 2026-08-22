#!/usr/bin/env python3
"""Fail if audited Taylor interaction work is missing or multiplied by L_eff/b."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--identity-rtol", type=float, default=1.0e-10)
    parser.add_argument("--require-order-one", action="store_true")
    args = parser.parse_args()

    rows = []
    for path in args.audit:
        with path.open() as handle:
            for line in handle:
                row = json.loads(line)
                if row.get("interaction_class") and row.get("kinetically_eligible") == 1:
                    rows.append(row)
    if not rows:
        raise SystemExit("no audited Taylor D-D interaction rows")

    identity_errors = []
    ratios = []
    phi_geom = []
    for row in rows:
        force = float(row["F_event_used_N"])
        x_dagger = float(row["x_dagger_m"])
        vstar = float(row["v_star_m3"])
        tau = float(row["tau_eff_Pa"])
        expected = force * x_dagger / vstar
        scale = max(abs(expected), abs(tau), 1.0)
        identity_errors.append(abs(tau - expected) / scale)

        tau_app = float(row.get("tau_app_resolved_Pa", 0.0))
        geom = float(row.get("phi_geom_L_over_b", 0.0))
        tau_phi = tau_app * geom
        if abs(tau_phi) > 0.0 and math.isfinite(tau_phi) and math.isfinite(tau):
            ratios.append(abs(tau / tau_phi))
            phi_geom.append(abs(geom))

    max_error = max(identity_errors)
    identity_ok = max_error <= args.identity_rtol
    median_ratio = float(np.median(ratios)) if ratios else None
    order_one = median_ratio is not None and 0.1 <= median_ratio <= 10.0
    slope = None
    if len(ratios) >= 3:
        x = np.log(np.maximum(np.asarray(phi_geom), 1.0e-300))
        y = np.log(np.maximum(np.asarray(ratios), 1.0e-300))
        if float(np.ptp(x)) > 1.0e-12:
            slope = float(np.polyfit(x, y, 1)[0])

    result = {
        "status": "passed" if identity_ok and (order_one or not args.require_order_one) else "failed",
        "interaction_rows": len(rows),
        "force_work_identity": "tau_eff = F_event_used_N * x_dagger_m / v_star_m3",
        "max_relative_identity_error": max_error,
        "identity_passed": identity_ok,
        "force_work_is_not_multiplied_by_L_eff_over_b": identity_ok,
        "force_work_over_phi_app_median": median_ratio,
        "force_work_over_phi_app_log_slope_vs_L_eff_over_b": slope,
        "clean_geometry_order_one_required": args.require_order_one,
        "clean_geometry_order_one_passed": order_one if args.require_order_one else None,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
