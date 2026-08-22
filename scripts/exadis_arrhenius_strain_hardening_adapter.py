#!/usr/bin/env python3
"""Prototype adapter for Arrhenius/TST ExaDiS strain-hardening work.

This script has two modes:

1. `--dry-run`: no ExaDiS dependency.  It writes the mechanism parameter file,
   analytical Taylor peak table, and implementation checklist.

2. `--run-stock-exadis`: attempts to import `pyexadis` and run a stock ExaDiS
   strain-hardening example structure.  This is intentionally audit-only at this
   stage.  It does not claim to replace native ExaDiS mobility/topology kernels.

The purpose is to provide a stable entry point for the branch while the native
Arrhenius mobility/topology hooks are developed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

# Allow use from repository root without installing as a package.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from arrhenius_tst_laws import (  # noqa: E402
    ArrheniusHazard,
    ExpFloorBarrier,
    analytical_taylor_exp_floor_stress_MPa,
    scan_peak_by_temperature,
)


DEFAULT_MECHANISMS = {
    "peierls_glide": {
        "barrier_family": "exp_floor",
        "H_eV": 0.05,
        "S_kB": 0.0,
        "sigma_c_GPa": 14.5,
        "floor_fraction": 0.0,
        "a": 6.65607,
        "n": 2.15276,
        "eta0_s": 1.0e12,
        "stress_convention": "resolved_glide_stress",
        "rate_convention": "signed_forward_minus_reverse",
    },
    "forest_depinning": {
        "barrier_family": "exp_floor",
        "H_eV": 0.50,
        "S_kB": -9.0,
        "sigma_c_GPa": 14.5,
        "floor_fraction": 0.20,
        "a": 6.65607,
        "n": 2.15276,
        "eta0_s": 1.0e12,
        "stress_convention": "force_work_tau_eff",
        "rate_convention": "hazard_release",
        "vstar_b3": 10.0,
        "x_dagger_rule": "vstar_over_b_squared",
    },
    "junction_unzip": {
        "barrier_family": "exp_floor",
        "H_eV": 0.50,
        "S_kB": -9.0,
        "sigma_c_GPa": 14.5,
        "floor_fraction": 0.20,
        "a": 6.65607,
        "n": 2.15276,
        "eta0_s": 1.0e12,
        "stress_convention": "reaction_force_work_tau_eff",
        "rate_convention": "hazard_unzip",
    },
    "cross_slip": {
        "barrier_family": "exp_floor",
        "H_eV": 0.50,
        "S_kB": -9.0,
        "sigma_c_GPa": 14.5,
        "floor_fraction": 0.20,
        "a": 6.65607,
        "n": 2.15276,
        "eta0_s": 1.0e12,
        "stress_convention": "primary_minus_cross_slip_resolved_stress",
        "rate_convention": "competing_hazard",
    },
}


def mechanism_to_law(mech: dict) -> ArrheniusHazard:
    barrier = ExpFloorBarrier(
        H_eV=mech["H_eV"],
        S_kB=mech["S_kB"],
        sigma_c_Pa=mech["sigma_c_GPa"] * 1.0e9,
        f=mech["floor_fraction"],
        a=mech["a"],
        n=mech["n"],
    )
    return ArrheniusHazard(barrier=barrier, eta0_s=mech.get("eta0_s", 1.0e12))


def write_default_config(outdir: Path) -> Path:
    config = {
        "name": "arrhenius_exadis_fcc_strain_hardening_seed_config",
        "description": "Initial TST parameter schema for ExaDiS/ParaDiS strain-hardening reparameterization.",
        "mechanisms": DEFAULT_MECHANISMS,
        "global": {
            "burgers_m": 2.55e-10,
            "strain_rate_s": 1.0e3,
            "temperatures_K": [850, 900, 950, 1000, 1050],
            "densities_m2_for_audit": [1e15, 2e15, 3e15, 5e15, 7e15, 1e16, 1.5e16, 2e16, 3e16],
        },
    }
    path = outdir / "arrhenius_mechanism_config.json"
    path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    return path


def write_analytical_peak_table(outdir: Path, args) -> Path:
    law = mechanism_to_law(DEFAULT_MECHANISMS["forest_depinning"])
    densities = [float(x) for x in args.densities]
    temps = [float(x) for x in args.temperatures]

    rows = []
    for T in temps:
        for rho in densities:
            tau, regime, G_req, tau_eff = analytical_taylor_exp_floor_stress_MPa(
                rho_m2=rho,
                T_K=T,
                strain_rate_s=args.strain_rate,
                b_m=args.burgers,
                law=law,
            )
            rows.append({
                "T_K": T,
                "rho_m2": rho,
                "tau_app_MPa": tau,
                "tau_eff_GPa": tau_eff,
                "regime": regime,
                "G_req_eV": G_req,
            })

    peak_rows = scan_peak_by_temperature(
        temperatures_K=temps,
        densities_m2=densities,
        strain_rate_s=args.strain_rate,
        b_m=args.burgers,
        law=law,
    )

    path = outdir / "analytical_taylor_peak_table.csv"
    with path.open("w") as handle:
        handle.write("T_K,rho_m2,tau_app_MPa,tau_eff_GPa,regime,G_req_eV\n")
        for row in rows:
            handle.write(
                f"{row['T_K']},{row['rho_m2']},{row['tau_app_MPa']},{row['tau_eff_GPa']},{row['regime']},{row['G_req_eV']}\n"
            )

    peak_path = outdir / "analytical_taylor_peak_summary.json"
    peak_path.write_text(json.dumps(peak_rows, indent=2, sort_keys=True) + "\n")
    return path


def write_checklist(outdir: Path) -> Path:
    text = """# Native ExaDiS implementation checklist

1. Build ExaDiS with Python bindings and verify the stock FCC strain-hardening example runs.
2. Add audit columns for every candidate event family:
   - local stress or force-work coordinate
   - barrier_eV
   - rate_s
   - probability_per_step
   - selected/not selected
3. Implement ArrheniusPeierls mobility using signed forward-minus-reverse rates.
4. Implement ArrheniusTopology for junction zip/unzip and forest depinning hazards.
5. Keep remeshing numerical only.
6. Compare stock ExaDiS, audit-only ExaDiS, Arrhenius mobility only, and full Arrhenius topology.
7. Reject any implementation that hides a deterministic Taylor stress or release threshold in a stress cap.
"""
    path = outdir / "native_exadis_checklist.md"
    path.write_text(text)
    return path


def run_stock_exadis_audit(args) -> None:
    try:
        import pyexadis  # type: ignore
        from pyexadis_base import (  # type: ignore
            ExaDisNet,
            DisNetManager,
            SimulateNetworkPerf,
            CalForce,
            MobilityLaw,
            TimeIntegration,
            Collision,
            Topology,
            Remesh,
        )
    except ImportError as exc:
        raise SystemExit(
            "Cannot import pyexadis.  Use --dry-run on machines without ExaDiS, "
            "or build ExaDiS with -DEXADIS_PYTHON_BINDING=On."
        ) from exc

    pyexadis.initialize()
    try:
        state = {
            "crystal": "fcc",
            "burgmag": args.burgers,
            "mu": 54.6e9,
            "nu": 0.324,
            "a": 6.0,
            "maxseg": 2000.0,
            "minseg": 300.0,
            "rtol": 10.0,
            "rann": 10.0,
            "nextdt": 1e-10,
            "maxdt": 1e-9,
        }
        data = Path(args.exadis_data)
        if not data.exists():
            raise SystemExit(f"ExaDiS data file not found: {data}")

        G = ExaDisNet()
        G.read_paradis(str(data))
        net = DisNetManager(G)

        calforce = CalForce(force_mode="SUBCYCLING_MODEL", state=state, Ngrid=64, cell=net.cell)
        mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
        timeint = TimeIntegration(integrator="Subcycling", rgroups=[0.0, 100.0, 600.0, 1600.0], state=state, force=calforce, mobility=mobility)
        collision = Collision(collision_mode="Retroactive", state=state)
        topology = Topology(topology_mode="TopologyParallel", state=state, force=calforce, mobility=mobility)
        remesh = Remesh(remesh_rule="LengthBased", state=state)

        sim = SimulateNetworkPerf(
            calforce=calforce,
            mobility=mobility,
            timeint=timeint,
            collision=collision,
            topology=topology,
            remesh=remesh,
            cross_slip=None,
            vis=None,
            loading_mode="strain_rate",
            erate=args.strain_rate,
            edir=np.array([0.0, 0.0, 1.0]),
            max_strain=args.max_strain,
            burgmag=state["burgmag"],
            state=state,
            print_freq=1,
            plot_freq=10,
            write_freq=100,
            write_dir=str(args.outdir / "stock_exadis_output"),
            restart=None,
        )
        sim.run(net, state)
    finally:
        pyexadis.finalize()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results/arrhenius_exadis_adapter"))
    parser.add_argument("--dry-run", action="store_true", help="write configs and analytical tables without importing ExaDiS")
    parser.add_argument("--run-stock-exadis", action="store_true", help="attempt stock ExaDiS strain-hardening run")
    parser.add_argument("--exadis-data", type=Path, default=Path("core/exadis/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data"))
    parser.add_argument("--max-strain", type=float, default=1.0e-3)
    parser.add_argument("--strain-rate", type=float, default=1.0e3)
    parser.add_argument("--burgers", type=float, default=2.55e-10)
    parser.add_argument("--temperatures", nargs="+", default=[850, 900, 950, 1000, 1050])
    parser.add_argument("--densities", nargs="+", default=[1e15, 2e15, 3e15, 5e15, 7e15, 1e16, 1.5e16, 2e16, 3e16])
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    config_path = write_default_config(args.outdir)
    peak_table = write_analytical_peak_table(args.outdir, args)
    checklist = write_checklist(args.outdir)

    print(f"wrote {config_path}")
    print(f"wrote {peak_table}")
    print(f"wrote {checklist}")

    if args.run_stock_exadis:
        run_stock_exadis_audit(args)
    elif not args.dry_run:
        print("No ExaDiS run requested.  Use --run-stock-exadis, or --dry-run to suppress this message.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
