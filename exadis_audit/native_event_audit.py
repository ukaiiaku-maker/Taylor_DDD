#!/usr/bin/env python3
"""Run stock or gate-passed Arrhenius ExaDiS with native audit off/on."""

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


def _network_line_diagnostics(system: Any, burgers_m: float) -> dict[str, float]:
    nodes = np.asarray(system.get_nodes_array(), dtype=float)
    segs = np.asarray(system.get_segs_array(), dtype=float)
    cell = system.get_cell()
    volume_m3 = float(cell.volume()) * burgers_m**3
    if segs.size == 0:
        return {
            "volume_m3": volume_m3,
            "mobile_line_density_m2": 0.0,
            "forest_intersecting_line_density_m2": 0.0,
            "junction_density_m3": 0.0,
            "mean_segment_length_m": 0.0,
        }
    nodeids = segs[:, :2].astype(int)
    positions = nodes[:, 2:5]
    r1 = positions[nodeids[:, 0]]
    r2 = np.asarray(cell.closest_image(Rref=r1, R=positions[nodeids[:, 1]]))
    lengths_m = np.linalg.norm(r2 - r1, axis=1) * burgers_m
    degrees = np.bincount(nodeids.ravel(), minlength=len(nodes))
    constraints = nodes[:, 5].astype(int)
    mobile = (constraints[nodeids[:, 0]] == 0) & (constraints[nodeids[:, 1]] == 0)
    junction_seg = np.linalg.norm(segs[:, 2:5], axis=1) > 1.01
    intersects = (
        (degrees[nodeids[:, 0]] >= 3) |
        (degrees[nodeids[:, 1]] >= 3) |
        junction_seg
    )
    return {
        "volume_m3": volume_m3,
        "mobile_line_density_m2": float(lengths_m[mobile].sum() / volume_m3),
        "forest_intersecting_line_density_m2": float(lengths_m[intersects].sum() / volume_m3),
        "junction_density_m3": float(np.count_nonzero(degrees >= 3) / volume_m3),
        "mean_segment_length_m": float(np.mean(lengths_m)),
    }


def _network_density_m2(graph: Any, burgers_m: float) -> float:
    """Compute line density for an ExaDisNet, including periodic images."""
    data = graph.export_data()
    cell_data = data["cell"]
    nodes = np.asarray(data["nodes"]["positions"], dtype=float)
    nodeids = np.asarray(data["segs"]["nodeids"], dtype=int)
    if nodeids.size == 0:
        return 0.0
    import pyexadis

    cell = pyexadis.Cell(
        h=cell_data["h"], origin=cell_data["origin"],
        is_periodic=cell_data["is_periodic"],
    )
    r1 = nodes[nodeids[:, 0]]
    r2 = np.asarray(cell.closest_image(Rref=r1, R=nodes[nodeids[:, 1]]))
    line_internal = float(np.linalg.norm(r2 - r1, axis=1).sum())
    volume_internal = float(cell.volume())
    return line_internal / (volume_internal * burgers_m * burgers_m)


def _scale_initial_density(graph: Any, density_factor: float) -> tuple[float, float]:
    """Uniformly scale cell and nodes so rho_new/rho_old=density_factor."""
    if not np.isfinite(density_factor) or density_factor <= 0.0:
        raise ValueError("--density-factor must be finite and positive")
    before = _network_density_m2(graph, 1.0)
    if density_factor == 1.0:
        return before, before
    data = graph.export_data()
    scale = density_factor ** -0.5
    h = np.asarray(data["cell"]["h"], dtype=float)
    origin = np.asarray(data["cell"]["origin"], dtype=float)
    center = origin + 0.5 * np.sum(h, axis=0)
    data["cell"]["h"] = scale * h
    data["cell"]["origin"] = center + scale * (origin - center)
    positions = np.asarray(data["nodes"]["positions"], dtype=float)
    data["nodes"]["positions"] = center + scale * (positions - center)
    graph.import_data(data)
    after = _network_density_m2(graph, 1.0)
    ratio = after / before
    if not np.isclose(ratio, density_factor, rtol=1.0e-10, atol=0.0):
        raise RuntimeError(
            f"initial-density scaling mismatch: requested {density_factor}, got {ratio}"
        )
    return before, after


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


def _discrete_audit_metrics(path: Path) -> dict[str, dict[str, Any]]:
    metrics: defaultdict[str, Counter[str]] = defaultdict(Counter)
    rdt_max: defaultdict[str, float] = defaultdict(float)
    if not path.exists():
        return {}
    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            mechanism = str(row.get("mechanism", ""))
            if mechanism not in {"topology_split", "cross_slip"}:
                continue
            metrics[mechanism]["rows"] += 1
            if row.get("geometry_admissible") == 1:
                metrics[mechanism]["geometry_admissible"] += 1
            if row.get("kinetically_eligible") == 1 and float(row.get("H_eV", 0.0)) > 0.0:
                metrics[mechanism]["kinetic_rows"] += 1
            if row.get("accepted_arrhenius") == 1:
                metrics[mechanism]["accepted_arrhenius"] += 1
            if row.get("deterministic_high_hazard") == 1:
                metrics[mechanism]["deterministic_high_hazard"] += 1
            rdt_max[mechanism] = max(rdt_max[mechanism], float(row.get("Rdt", 0.0)))
    return {
        mechanism: {**dict(counts), "Rdt_max": rdt_max[mechanism]}
        for mechanism, counts in sorted(metrics.items())
    }


def _arrhenius_mobility_params(args: argparse.Namespace) -> dict[str, Any]:
    if args.arrhenius_mobility == "off":
        return {}
    if args.arrhenius_config_data is None:
        raise ValueError("--arrhenius-mobility requires --arrhenius-config")
    block_name = "mobility_peierls"
    block = args.arrhenius_config_data.get(block_name)
    if not isinstance(block, dict):
        raise ValueError(f"Arrhenius config lacks object {block_name!r}")
    if block.get("replacement_eligible") is not True:
        raise ValueError(
            "Arrhenius mobility parameter block is not marked replacement_eligible=true"
        )
    required = ("H_eV", "S_kB", "sigma_c_GPa", "f", "a", "n", "jump_b", "vstar_b3")
    missing = [name for name in required if block.get(name) is None]
    if missing:
        raise ValueError(f"Arrhenius mobility config is incomplete: {', '.join(missing)}")
    coupling = block.get("anisotropic_coupling", {})
    params = {
        "mode": args.arrhenius_mobility,
        "temperature_K": args.arrhenius_temperature_K,
        "eta0_s": args.arrhenius_eta0_default,
        "vmax": float(block.get("vmax", 4000.0)),
        **{name: float(block[name]) for name in required},
        "a_nn": float(coupling.get("a_nn", 0.0)),
        "a_mm": float(coupling.get("a_mm", 0.0)),
        "a_np": float(coupling.get("a_np", 0.0)),
        "a_other": float(coupling.get("a_other", 0.0)),
    }
    site_multiplicity = float(block.get("site_multiplicity", 1.0))
    stress_concentration = float(block.get("stress_concentration_phi", 1.0))
    intrinsic_vstar = float(
        block.get("intrinsic_activation_volume_b3", block["vstar_b3"])
    )
    if site_multiplicity <= 0.0:
        raise ValueError("Arrhenius mobility site_multiplicity must be positive")
    if not np.isclose(
        stress_concentration * intrinsic_vstar, params["vstar_b3"], rtol=1e-12
    ):
        raise ValueError(
            "Arrhenius mobility requires vstar_b3 = "
            "stress_concentration_phi * intrinsic_activation_volume_b3"
        )
    prefactor_exponent = float(
        block.get("attempt_frequency_temperature_exponent", 0.0)
    )
    prefactor_reference = float(
        block.get("attempt_frequency_reference_temperature_K", 900.0)
    )
    if prefactor_reference <= 0.0:
        raise ValueError("Arrhenius mobility prefactor reference must be positive")
    params["eta0_s"] = (
        site_multiplicity
        * args.arrhenius_eta0_default
        * (args.arrhenius_temperature_K / prefactor_reference) ** prefactor_exponent
    )
    if args.arrhenius_mobility == "full":
        for name in ("H_screw_eV", "sigma_c_screw_GPa", "vstar_screw_b3"):
            if block.get(name) is None:
                raise ValueError(f"full Arrhenius mobility config lacks {name}")
            params[name] = float(block[name])
    return params


def _arrhenius_discrete_params(
    args: argparse.Namespace, mechanism: str, block: dict[str, Any], seed_offset: int
) -> Any:
    import pyexadis

    required = ("H_eV", "S_kB", "sigma_c_GPa", "f", "a", "n", "vstar_b3")
    missing = [name for name in required if block.get(name) is None]
    if missing:
        raise ValueError(
            f"Arrhenius {mechanism} config is incomplete: {', '.join(missing)}"
        )
    coupling = block.get("anisotropic_coupling", {})
    site_multiplicity = float(block.get("site_multiplicity", 1.0))
    if site_multiplicity <= 0.0:
        raise ValueError(f"Arrhenius {mechanism} site_multiplicity must be positive")
    if "stress_concentration_phi" in block:
        raise ValueError(
            f"Arrhenius {mechanism} must not impose stress_concentration_phi; "
            "the native Taylor force-work kernel audits phi_eff"
        )
    if block.get("stress_concentration_mode") != "single_glider_line_tension_force_work":
        raise ValueError(
            f"Arrhenius {mechanism} requires "
            "stress_concentration_mode=single_glider_line_tension_force_work"
        )
    if block.get("L_eff_mode") != "harmonic_adjacent_arms":
        raise ValueError(
            f"Arrhenius {mechanism} requires L_eff_mode=harmonic_adjacent_arms"
        )
    force_source_name = block.get(
        "force_work_source", "native_trial_force_preferred"
    )
    force_source_values = {
        "native_trial_force_preferred": 0,
        "line_tension_reconstruction": 1,
    }
    if force_source_name not in force_source_values:
        raise ValueError(
            f"Arrhenius {mechanism} has unsupported force_work_source={force_source_name}"
        )
    x_dagger_b = float(block.get("x_dagger_b", 0.0))
    line_tension_alpha = float(block.get("line_tension_alpha", 0.0))
    if x_dagger_b <= 0.0 or line_tension_alpha <= 0.0:
        raise ValueError(
            f"Arrhenius {mechanism} requires positive x_dagger_b and line_tension_alpha"
        )
    prefactor_exponent = float(
        block.get("attempt_frequency_temperature_exponent", 0.0)
    )
    prefactor_reference = float(
        block.get("attempt_frequency_reference_temperature_K", 900.0)
    )
    if prefactor_reference <= 0.0:
        raise ValueError(f"Arrhenius {mechanism} prefactor reference must be positive")
    eta0_effective = (
        site_multiplicity
        * float(block.get("eta0_s", args.arrhenius_eta0_default))
        * (args.arrhenius_temperature_K / prefactor_reference) ** prefactor_exponent
    )
    return pyexadis.Arrhenius_DiscreteEvent_Params(
        True,
        float(args.arrhenius_temperature_K),
        float(block["H_eV"]),
        float(block["S_kB"]),
        float(block["sigma_c_GPa"]),
        float(block["f"]),
        float(block["a"]),
        float(block["n"]),
        eta0_effective,
        float(block["vstar_b3"]),
        float(args.burgers),
        float(coupling.get("a_nn", 0.0)),
        float(coupling.get("a_mm", 0.0)),
        float(coupling.get("a_np", 0.0)),
        float(coupling.get("a_tc", coupling.get("a_other", 0.0))),
        float(block.get("high_hazard_Rdt", 20.0)),
        int(block.get("seed", 1469598103934665603 + seed_offset)),
        x_dagger_b,
        line_tension_alpha,
        force_source_values[force_source_name],
    )


def _arrhenius_event_family(
    args: argparse.Namespace, family: str, names: tuple[str, ...]
) -> list[Any]:
    if args.arrhenius_config_data is None:
        raise ValueError(f"--arrhenius-{family.replace('_', '-')} requires --arrhenius-config")
    family_block = args.arrhenius_config_data.get(family)
    if not isinstance(family_block, dict) or family_block.get("replacement_eligible") is not True:
        raise ValueError(
            f"Arrhenius {family} block is not marked replacement_eligible=true"
        )
    mechanisms = family_block.get("mechanisms")
    if not isinstance(mechanisms, dict):
        raise ValueError(f"Arrhenius {family} block lacks mechanisms")
    missing = [name for name in names if not isinstance(mechanisms.get(name), dict)]
    if missing:
        raise ValueError(f"Arrhenius {family} lacks mechanism blocks: {', '.join(missing)}")
    parameter_sets = args.arrhenius_config_data.get("interaction_parameter_sets", {})
    resolved = []
    for index, name in enumerate(names):
        block = mechanisms[name]
        parameter_set = args.interaction_parameter_set or block.get("parameter_set")
        if parameter_set is not None:
            base = parameter_sets.get(parameter_set)
            if not isinstance(base, dict):
                raise ValueError(
                    f"Arrhenius {name} references missing interaction parameter set "
                    f"{parameter_set}"
                )
            block = ({**block, **base} if args.interaction_parameter_set
                     else {**base, **block})
        resolved.append(_arrhenius_discrete_params(args, name, block, index))
    return resolved


def _build_modules(state: dict[str, Any], net: Any, args: argparse.Namespace):
    import pyexadis
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
    if args.arrhenius_mobility == "off":
        mobility = MobilityLaw(
            mobility_law="FCC_0", state=state,
            Medge=64103.0, Mscrew=64103.0, vmax=4000.0,
        )
    else:
        mobility = MobilityLaw(
            mobility_law="FCC_0_ARRHENIUS", state=state,
            **_arrhenius_mobility_params(args),
        )
    integrator = TimeIntegration(
        integrator="Subcycling", rgroups=[0.0, 100.0, 600.0, 1600.0],
        state=state, force=force, mobility=mobility,
    )
    collision = Collision(collision_mode="Retroactive", state=state)
    topology_kwargs: dict[str, Any] = {"force": force, "mobility": mobility}
    if args.arrhenius_topology == "on":
        topology_kwargs["arrhenius_events"] = _arrhenius_event_family(
            args,
            "topology",
            (
                "junction_zip",
                "junction_unzip",
                "junction_destruction",
                "junction_reconfiguration",
                "forest_depinning_like_release",
            ),
        )
    topology = Topology(
        topology_mode="TopologyParallel", state=state, **topology_kwargs
    )
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    cross = None
    if args.cross_slip or args.arrhenius_cross_slip == "on":
        cross_kwargs: dict[str, Any] = {"force": force}
        if args.arrhenius_cross_slip == "on":
            cross_kwargs["arrhenius_events"] = _arrhenius_event_family(
                args,
                "cross_slip",
                ("plane_change", "zipper_propagation"),
            )
        cross = CrossSlip(
            cross_slip_mode="ForceBasedParallel", state=state, **cross_kwargs
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
        "nextdt": min(
            1.0e-10,
            args.max_strain / (args.minimum_steps * args.strain_rate)
            if args.minimum_steps > 0 else 1.0e-10,
        ),
        "maxdt": min(
            1.0e-9,
            args.max_strain / (args.minimum_steps * args.strain_rate)
            if args.minimum_steps > 0 else 1.0e-9,
        ),
    }

    graph = ExaDisNet()
    graph.read_paradis(str(args.exadis_data), verbose=False)
    _base_density_internal, scaled_density_internal = _scale_initial_density(
        graph, args.density_factor
    )
    initial_density_m2 = scaled_density_internal / (args.burgers * args.burgers)
    net = DisNetManager(graph)
    modules = _build_modules(state, net, args)
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
    discrete_metrics = _discrete_audit_metrics(audit_path)
    line_diagnostics = _network_line_diagnostics(system, args.burgers)
    summary = {
        "audit_enabled": audit_enabled,
        "density_factor": float(args.density_factor),
        "initial_density_m2": float(initial_density_m2),
        "strain": float(state["strain"]),
        "stress_Pa": float(state["stress"]),
        "pstrain": float(state["pstrain"]),
        "density_m2": float(state["density"]),
        "Nnodes": int(system.number_of_nodes()),
        "Nsegs": int(system.number_of_segs()),
        "network_sane": bool(system.is_sane()),
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
        "arrhenius_discrete_audit": discrete_metrics,
        **line_diagnostics,
    }
    (outdir / "final_summary.json").write_text(
        json.dumps(_clean(summary), indent=2, sort_keys=True) + "\n"
    )
    return summary


def _compare(disabled: dict[str, Any], enabled: dict[str, Any], tolerance: float) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    passed = True
    for key in (
        "strain", "stress_Pa", "pstrain", "initial_density_m2", "density_m2",
        "density_factor", "dt_s", "time_s",
    ):
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
        description="Run stock or gate-passed native Arrhenius FCC ExaDiS."
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
    parser.add_argument(
        "--minimum-steps", type=int, default=0,
        help="cap maxdt so the target strain uses at least this many steps",
    )
    parser.add_argument("--burgers", type=float, default=2.55e-10)
    parser.add_argument(
        "--density-factor", type=float, default=1.0,
        help="target initial line-density factor via uniform geometry scaling",
    )
    parser.add_argument("--audit-stride", type=int, default=1)
    parser.add_argument("--print-freq", type=int, default=1)
    parser.add_argument("--output-frequency", type=int, default=2_000_000_000)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    parser.add_argument(
        "--audit-enabled-only", action="store_true",
        help="run one audited case for calibration/validation; do not claim invariance",
    )
    parser.add_argument("--cross-slip", action="store_true")
    parser.add_argument(
        "--arrhenius-mobility", choices=("off", "peierls", "full"), default="off"
    )
    parser.add_argument("--arrhenius-topology", choices=("off", "on"), default="off")
    parser.add_argument("--arrhenius-cross-slip", choices=("off", "on"), default="off")
    parser.add_argument(
        "--arrhenius-collision", choices=("off", "activated-only"), default="off"
    )
    parser.add_argument("--arrhenius-temperature-K", type=float, default=None)
    parser.add_argument("--arrhenius-eta0-default", type=float, default=1.0e12)
    parser.add_argument("--arrhenius-config", type=Path, default=None)
    parser.add_argument(
        "--interaction-parameter-set", default=None,
        help="override every D-D mechanism with a named interaction_parameter_sets block",
    )
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
    args.arrhenius_config_data = None
    if args.arrhenius_config is not None:
        config_path = (
            args.repo_root / args.arrhenius_config
            if not args.arrhenius_config.is_absolute() else args.arrhenius_config
        ).resolve()
        args.arrhenius_config_data = json.loads(config_path.read_text())
        if args.arrhenius_temperature_K is None:
            args.arrhenius_temperature_K = float(args.arrhenius_config_data.get("temperature_K", 0.0))
    if args.arrhenius_mobility != "off" and not (args.arrhenius_temperature_K and args.arrhenius_temperature_K > 0.0):
        raise SystemExit("requested Arrhenius mobility requires a positive temperature")
    if args.minimum_steps < 0:
        raise SystemExit("--minimum-steps must be nonnegative")
    if args.strain_rate <= 0.0 or args.max_strain <= 0.0:
        raise SystemExit("--strain-rate and --max-strain must be positive")
    unsupported = []
    if args.arrhenius_collision == "activated-only": unsupported.append("collision")
    if unsupported:
        raise SystemExit(
            "requested Arrhenius module has not passed its replacement gate: " + ", ".join(unsupported)
        )

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
        enabled_dir = args.outdir / "audit_enabled"
        enabled_dir.mkdir()
        if args.audit_enabled_only:
            disabled = None
        else:
            disabled_dir = args.outdir / "audit_disabled"
            disabled_dir.mkdir()
            disabled = _run_case(args, disabled_dir, audit_enabled=False)
        enabled = _run_case(args, enabled_dir, audit_enabled=True)
    finally:
        pyexadis.finalize()

    comparison = None if disabled is None else _compare(disabled, enabled, args.tolerance)
    comparison_plot = None
    if disabled is not None:
        comparison_plot = args.outdir / "audit_comparison.svg"
        _write_comparison_svg(
            disabled_dir / "stress_strain_dens.dat",
            enabled_dir / "stress_strain_dens.dat",
            comparison_plot,
        )
    report = {
        "status": (
            "calibration_trace_complete" if comparison is None else
            ("passed" if comparison["passed"] else "failed")
        ),
        "instrumentation_only": (
            args.arrhenius_mobility == "off"
            and args.arrhenius_topology == "off"
            and args.arrhenius_cross_slip == "off"
        ),
        "arrhenius_replacements_connected": (
            args.arrhenius_mobility != "off"
            or args.arrhenius_topology == "on"
            or args.arrhenius_cross_slip == "on"
        ),
        "arrhenius_stage": (
            "A3_arrhenius_mobility_topology_cross_slip"
            if args.arrhenius_cross_slip == "on" else
            "A2_arrhenius_mobility_topology"
            if args.arrhenius_topology == "on" else
            "A0_stock" if args.arrhenius_mobility == "off" else
            "A1_arrhenius_peierls_only"
        ),
        "repository_sha": _git_sha(args.repo_root),
        "exadis_source_sha": _git_sha(args.exadis_root),
        "python_executable": sys.executable,
        "pyexadis_module": pyexadis.__file__,
        "kokkos_backend": "OpenMP+Serial",
        "openmp_threads": threads,
        "mpi_used": False,
        "comparison_plot": str(comparison_plot) if comparison_plot else None,
        "audit_enabled_only": args.audit_enabled_only,
        "audit_disabled": disabled,
        "audit_enabled": enabled,
        "invariance": comparison,
    }
    (args.outdir / "invariance_summary.json").write_text(
        json.dumps(_clean(report), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(_clean(comparison), indent=2, sort_keys=True))
    if comparison is not None and not comparison["passed"]:
        raise SystemExit("native audit-on/audit-off invariance failed")
    expected_mobility = (
        "mobility_fcc0" if args.arrhenius_mobility == "off" else
        "mobility_fcc0_arrhenius"
    )
    if enabled["audit_rows_by_mechanism"].get(expected_mobility, 0) == 0:
        raise SystemExit(f"native audit produced no {expected_mobility} rows")
    for mechanism in args.require_candidate_labels:
        counts = enabled["candidate_labels_by_mechanism"].get(mechanism, {})
        if counts.get("accepted", 0) == 0 or counts.get("rejected", 0) == 0:
            raise SystemExit(
                f"native audit did not produce accepted and rejected labels for {mechanism}: {counts}"
            )
    for requested, mechanism in (
        (args.arrhenius_topology == "on", "topology_split"),
        (args.arrhenius_cross_slip == "on", "cross_slip"),
    ):
        if requested:
            metrics = enabled["arrhenius_discrete_audit"].get(mechanism, {})
            if metrics.get("kinetic_rows", 0) == 0:
                raise SystemExit(
                    f"requested Arrhenius {mechanism} produced no audited kinetic rows: {metrics}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
