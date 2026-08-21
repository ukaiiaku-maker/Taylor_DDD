#!/usr/bin/env python3
"""Run the stock high-performance ExaDiS driver with native audit off/on.

This is an instrumentation-only gate.  It never installs an Arrhenius law and
fails if the runtime-enabled native recorder changes the stock trajectory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


def _clean(value: Any) -> Any:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    return value


def _git_sha(path: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _network_digest(system: Any) -> str:
    digest = hashlib.sha256()
    nodes = np.asarray(system.get_nodes_array(), dtype=np.float64)
    segs = np.asarray(system.get_segs_array(), dtype=np.float64)
    digest.update(nodes.tobytes(order="C"))
    digest.update(segs.tobytes(order="C"))
    return digest.hexdigest()


def _file_digest(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _line_count(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def _normalized_restart_digest(path: Path) -> str | None:
    """Hash restart content after removing ExaDiS's wall-clock timestamp."""
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b"date_and_time "):
                continue
            digest.update(line)
    return digest.hexdigest()


def _write_comparison_svg(disabled_path: Path, enabled_path: Path, output: Path) -> None:
    """Write a dependency-free stress/strain and density comparison plot."""
    disabled = np.atleast_2d(np.loadtxt(disabled_path, comments="#"))
    enabled = np.atleast_2d(np.loadtxt(enabled_path, comments="#"))
    width, height = 960, 420
    left, right, top, bottom = 70, 25, 35, 55
    gap = 70
    panel_width = (width - left - right - gap) / 2
    panel_height = height - top - bottom

    def panel_points(data: np.ndarray, xcol: int, ycol: int, x0: float,
                     xmin: float, xmax: float, ymin: float, ymax: float) -> str:
        xr = max(xmax - xmin, np.finfo(float).eps)
        yr = max(ymax - ymin, np.finfo(float).eps)
        points = []
        for row in data:
            x = x0 + panel_width * (row[xcol] - xmin) / xr
            y = top + panel_height * (1.0 - (row[ycol] - ymin) / yr)
            points.append(f"{x:.2f},{y:.2f}")
        return " ".join(points)

    panels = []
    for index, (column, title, ylabel) in enumerate((
        (2, "Stress–strain", "Stress (Pa)"),
        (3, "Density–strain", "Density (m⁻²)"),
    )):
        x0 = left + index * (panel_width + gap)
        xmin = float(min(disabled[:, 1].min(), enabled[:, 1].min()))
        xmax = float(max(disabled[:, 1].max(), enabled[:, 1].max()))
        ymin = float(min(disabled[:, column].min(), enabled[:, column].min()))
        ymax = float(max(disabled[:, column].max(), enabled[:, column].max()))
        if ymin == ymax:
            ymin -= max(abs(ymin), 1.0) * 0.01
            ymax += max(abs(ymax), 1.0) * 0.01
        off = panel_points(disabled, 1, column, x0, xmin, xmax, ymin, ymax)
        on = panel_points(enabled, 1, column, x0, xmin, xmax, ymin, ymax)
        panels.append(f'''<g>
  <rect x="{x0:.1f}" y="{top}" width="{panel_width:.1f}" height="{panel_height}" fill="white" stroke="#9ca3af"/>
  <polyline points="{off}" fill="none" stroke="#2563eb" stroke-width="2.2"/>
  <polyline points="{on}" fill="none" stroke="#dc2626" stroke-width="1.6" stroke-dasharray="6 4"/>
  <text x="{x0 + panel_width/2:.1f}" y="22" text-anchor="middle" font-size="16">{title}</text>
  <text x="{x0 + panel_width/2:.1f}" y="{height - 14}" text-anchor="middle" font-size="13">Strain</text>
  <text x="{x0 - 48:.1f}" y="{top + panel_height/2:.1f}" text-anchor="middle" font-size="13" transform="rotate(-90 {x0 - 48:.1f} {top + panel_height/2:.1f})">{ylabel}</text>
  <text x="{x0:.1f}" y="{height - bottom + 18}" font-size="11">{xmin:.3g}</text>
  <text x="{x0 + panel_width:.1f}" y="{height - bottom + 18}" text-anchor="end" font-size="11">{xmax:.3g}</text>
  <text x="{x0 - 8:.1f}" y="{top + 4}" text-anchor="end" font-size="11">{ymax:.3g}</text>
  <text x="{x0 - 8:.1f}" y="{top + panel_height}" text-anchor="end" font-size="11">{ymin:.3g}</text>
</g>''')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#f8fafc"/>
<g font-family="system-ui, sans-serif" fill="#111827">{''.join(panels)}
<line x1="350" y1="402" x2="375" y2="402" stroke="#2563eb" stroke-width="2.2"/><text x="380" y="406" font-size="12">audit off</text>
<line x1="470" y1="402" x2="495" y2="402" stroke="#dc2626" stroke-width="1.6" stroke-dasharray="6 4"/><text x="500" y="406" font-size="12">audit on</text>
</g></svg>\n'''
    output.write_text(svg)


def _audit_counts(path: Path) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    mechanisms: Counter[str] = Counter()
    labels: defaultdict[str, Counter[str]] = defaultdict(Counter)
    if not path.exists():
        return {}, {}
    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            mechanism = str(row.get("mechanism", "unknown"))
            mechanisms[mechanism] += 1
            accepted = row.get("accepted_stock", -1)
            if accepted in (0, 1):
                labels[mechanism]["accepted" if accepted else "rejected"] += 1
    return dict(sorted(mechanisms.items())), {
        key: dict(sorted(counts.items())) for key, counts in sorted(labels.items())
    }


def _build_modules(state: dict[str, Any], net: Any, cross_slip: bool):
    from pyexadis_base import (
        CalForce,
        Collision,
        CrossSlip,
        MobilityLaw,
        Remesh,
        TimeIntegration,
        Topology,
    )

    force = CalForce(force_mode="SUBCYCLING_MODEL", state=state, Ngrid=64, cell=net.cell)
    mobility = MobilityLaw(
        mobility_law="FCC_0", state=state,
        Medge=64103.0, Mscrew=64103.0, vmax=4000.0,
    )
    integrator = TimeIntegration(
        integrator="Subcycling", rgroups=[0.0, 100.0, 600.0, 1600.0],
        state=state, force=force, mobility=mobility,
    )
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(
        topology_mode="TopologyParallel", state=state,
        force=force, mobility=mobility,
    )
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    cross = (
        CrossSlip(cross_slip_mode="ForceBasedParallel", state=state, force=force)
        if cross_slip else None
    )
    return force, mobility, integrator, collision, topology, remesh, cross


def _run_case(args: argparse.Namespace, outdir: Path, audit_enabled: bool) -> dict[str, Any]:
    import pyexadis
    from pyexadis_base import (
        DisNetManager,
        ExaDisNet,
        SimulateNetworkPerf,
        get_exadis_params,
    )

    state: dict[str, Any] = {
        "crystal": "fcc",
        "burgmag": args.burgers,
        "mu": 54.6e9,
        "nu": 0.324,
        "a": 6.0,
        "maxseg": 2000.0,
        "minseg": 300.0,
        "rtol": 10.0,
        "rann": 10.0,
        "nextdt": 1.0e-10,
        "maxdt": 1.0e-9,
    }

    graph = ExaDisNet()
    graph.read_paradis(str(args.exadis_data), verbose=False)
    net = DisNetManager(graph)
    modules = _build_modules(state, net, args.cross_slip)
    force, mobility, integrator, collision, topology, remesh, cross = modules

    sim = SimulateNetworkPerf(
        calforce=force,
        mobility=mobility,
        timeint=integrator,
        collision=collision,
        topology=topology,
        remesh=remesh,
        cross_slip=cross,
        vis=None,
        loading_mode="strain_rate",
        erate=args.strain_rate,
        edir=np.array([0.0, 0.0, 1.0]),
        max_strain=args.max_strain,
        burgmag=state["burgmag"],
        state=state,
        print_freq=args.print_freq,
        plot_freq=None,
        write_freq=args.output_frequency,
        write_dir=str(outdir),
        out_props=["Step", "Strain", "Stress", "Density", "Nnodes", "Nsegs", "DT", "Time"],
    )

    params = get_exadis_params(state)
    system = pyexadis.System(net.get_disnet(ExaDisNet).net, params)
    system.set_neighbor_cutoff(force.force.neighbor_cutoff)
    driver = pyexadis.Driver(system)
    native_modules = [
        force.force, mobility.mobility, integrator.integrator,
        collision.collision, topology.topology, remesh.remesh,
    ]
    if cross is not None:
        native_modules.append(cross.cross_slip)
    driver.set_modules(*native_modules)
    driver.outputdir = str(outdir)
    driver.set_simulation("")

    audit_path = outdir / "event_audit.jsonl"
    if audit_enabled:
        driver.enable_audit(str(audit_path), args.audit_stride)
        if not driver.audit_enabled():
            raise RuntimeError("driver.enable_audit() did not enable the native recorder")

    ctrl = sim.get_exadis_ctrl(state)
    driver.initialize(ctrl)
    stepper = pyexadis.Driver.MAX_STRAIN(args.max_strain)
    while stepper.iterate(driver):
        driver.step(ctrl)

    state = driver.update_state(state)
    if audit_enabled:
        driver.disable_audit()

    mechanisms, labels = _audit_counts(audit_path)
    summary = {
        "audit_enabled": audit_enabled,
        "strain": float(state["strain"]),
        "stress_Pa": float(state["stress"]),
        "pstrain": float(state["pstrain"]),
        "density_m2": float(state["density"]),
        "Nnodes": int(system.number_of_nodes()),
        "Nsegs": int(system.number_of_segs()),
        "dt_s": float(state["dt"]),
        "time_s": float(state["time"]),
        "istep": int(state["istep"]),
        "network_sha256": _network_digest(system),
        "stress_strain_rows": _line_count(outdir / "stress_strain_dens.dat"),
        "stress_strain_sha256": _file_digest(outdir / "stress_strain_dens.dat"),
        "config_0_sha256": _file_digest(outdir / "config.0.data"),
        "restart_0_sha256": _file_digest(outdir / "restart.0.exadis"),
        "restart_0_normalized_sha256": _normalized_restart_digest(outdir / "restart.0.exadis"),
        "audit_path": str(audit_path) if audit_enabled else None,
        "audit_rows_by_mechanism": mechanisms,
        "candidate_labels_by_mechanism": labels,
    }
    (outdir / "final_summary.json").write_text(
        json.dumps(_clean(summary), indent=2, sort_keys=True) + "\n"
    )
    return summary


def _compare(disabled: dict[str, Any], enabled: dict[str, Any], tolerance: float) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    passed = True
    for key in ("strain", "stress_Pa", "pstrain", "density_m2", "dt_s", "time_s"):
        left = float(disabled[key])
        right = float(enabled[key])
        absolute = abs(right - left)
        scale = max(1.0, abs(left), abs(right))
        ok = absolute <= tolerance * scale
        passed = passed and ok
        checks[key] = {
            "audit_disabled": left,
            "audit_enabled": right,
            "abs_diff": absolute,
            "rel_diff": absolute / scale,
            "ok": ok,
        }
    for key in (
        "Nnodes", "Nsegs", "istep", "network_sha256", "stress_strain_rows",
        "stress_strain_sha256", "config_0_sha256", "restart_0_normalized_sha256",
    ):
        left = disabled[key]
        right = enabled[key]
        ok = left == right
        passed = passed and ok
        checks[key] = {"audit_disabled": left, "audit_enabled": right, "ok": ok}
    return {"passed": passed, "tolerance": tolerance, "checks": checks}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prove runtime native audit on/off invariance for stock FCC ExaDiS."
    )
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--exadis-root", type=Path, default=Path("core/exadis"))
    parser.add_argument(
        "--exadis-data", type=Path,
        default=Path("core/exadis/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("results/exadis_native_audit"))
    parser.add_argument("--max-strain", type=float, default=1.0e-5)
    parser.add_argument("--strain-rate", type=float, default=1.0e3)
    parser.add_argument("--burgers", type=float, default=2.55e-10)
    parser.add_argument("--audit-stride", type=int, default=1)
    parser.add_argument("--print-freq", type=int, default=1)
    parser.add_argument("--output-frequency", type=int, default=2_000_000_000)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    parser.add_argument("--cross-slip", action="store_true")
    parser.add_argument(
        "--require-candidate-labels",
        action="append",
        choices=("cross_slip", "collision", "topology_split"),
        default=[],
        help="fail unless this mechanism has both accepted and rejected stock candidates",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    args.repo_root = args.repo_root.resolve()
    args.exadis_root = (args.repo_root / args.exadis_root).resolve() if not args.exadis_root.is_absolute() else args.exadis_root
    args.exadis_data = (args.repo_root / args.exadis_data).resolve() if not args.exadis_data.is_absolute() else args.exadis_data
    args.outdir = (args.repo_root / args.outdir).resolve() if not args.outdir.is_absolute() else args.outdir

    if not args.exadis_data.exists():
        raise SystemExit(f"ExaDiS input data not found: {args.exadis_data}")

    import pyexadis
    if not getattr(pyexadis, "NATIVE_EVENT_AUDIT_COMPILED", False):
        raise SystemExit("pyexadis was built without EXADIS_ENABLE_EVENT_AUDIT=ON")
    if not hasattr(pyexadis.Driver(), "enable_audit"):
        raise SystemExit("native audit pybind smoke test failed: Driver.enable_audit missing")

    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True)

    threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
    pyexadis.initialize(num_threads=threads, verbose=True)
    try:
        disabled_dir = args.outdir / "audit_disabled"
        enabled_dir = args.outdir / "audit_enabled"
        disabled_dir.mkdir()
        enabled_dir.mkdir()
        disabled = _run_case(args, disabled_dir, audit_enabled=False)
        enabled = _run_case(args, enabled_dir, audit_enabled=True)
    finally:
        pyexadis.finalize()

    comparison = _compare(disabled, enabled, args.tolerance)
    comparison_plot = args.outdir / "audit_comparison.svg"
    _write_comparison_svg(
        disabled_dir / "stress_strain_dens.dat",
        enabled_dir / "stress_strain_dens.dat",
        comparison_plot,
    )
    report = {
        "status": "passed" if comparison["passed"] else "failed",
        "instrumentation_only": True,
        "arrhenius_replacements_connected": False,
        "repository_sha": _git_sha(args.repo_root),
        "exadis_source_sha": _git_sha(args.exadis_root),
        "python_executable": sys.executable,
        "pyexadis_module": pyexadis.__file__,
        "kokkos_backend": "OpenMP+Serial",
        "openmp_threads": threads,
        "mpi_used": False,
        "comparison_plot": str(comparison_plot),
        "audit_disabled": disabled,
        "audit_enabled": enabled,
        "invariance": comparison,
    }
    (args.outdir / "invariance_summary.json").write_text(
        json.dumps(_clean(report), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(_clean(comparison), indent=2, sort_keys=True))
    if not comparison["passed"]:
        raise SystemExit("native audit-on/audit-off invariance failed")
    if enabled["audit_rows_by_mechanism"].get("mobility_fcc0", 0) == 0:
        raise SystemExit("native audit produced no mobility_fcc0 rows")
    for mechanism in args.require_candidate_labels:
        counts = enabled["candidate_labels_by_mechanism"].get(mechanism, {})
        if counts.get("accepted", 0) == 0 or counts.get("rejected", 0) == 0:
            raise SystemExit(
                f"native audit did not produce accepted and rejected labels for {mechanism}: {counts}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
