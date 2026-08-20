#!/usr/bin/env python3
"""Binding-level ExaDiS event-conjugate audit.

This module instruments the stock ExaDiS Python stepping path without changing
accepted mechanics.  It is intentionally audit-only.  It does not connect
Arrhenius hazards or alter mobility, topology, collision, cross-slip, or remesh
selection rules.

The performance driver, ``SimulateNetworkPerf``, executes the stock C++ driver
without Python per-module callbacks.  This file therefore uses the stock Python
``SimulateNetwork`` stepping path, whose hook methods expose the sequence

    force -> mobility -> integration -> cross slip -> collision -> topology -> remesh -> response

using the same ExaDiS modules.  This is the binding-level instrumentation path.
A later native C++ patch is still needed to expose hidden device-side trial
candidate internals in ``TopologyParallel``, ``CrossSlipParallel``, and
``CollisionRetroactive``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import numpy as np


# -----------------------------------------------------------------------------
# Audit recorder
# -----------------------------------------------------------------------------

@dataclass
class EventAuditRecorder:
    """JSONL event audit recorder.

    The recorder is disabled by default.  Passing a path enables it.  This is the
    binding-level equivalent of ``enable_audit(path, stride)``.
    """

    path: Optional[Path] = None
    stride: int = 1
    enabled: bool = False

    def __post_init__(self) -> None:
        if self.path is not None:
            self.path = Path(self.path)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._handle = self.path.open("w")
            self.enabled = True
        else:
            self._handle = None
            self.enabled = False
        self.rows = 0

    def should_record(self, step: int) -> bool:
        return self.enabled and self.stride > 0 and (int(step) % int(self.stride) == 0)

    def write(self, row: Dict[str, Any]) -> None:
        if not self.enabled:
            return
        clean = {}
        for key, value in row.items():
            clean[key] = _json_clean(value)
        self._handle.write(json.dumps(clean, sort_keys=True) + "\n")
        self.rows += 1

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> "EventAuditRecorder":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def enable_audit(path: Optional[os.PathLike[str] | str], stride: int = 1) -> EventAuditRecorder:
    """Return an enabled or disabled audit recorder.

    ``path=None`` leaves auditing disabled.  This mirrors the requested API while
    keeping audit disabled by default.
    """
    return EventAuditRecorder(path=None if path is None else Path(path), stride=stride)


def _json_clean(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, float)):
        val = float(value)
        if math.isnan(val):
            return None
        if math.isinf(val):
            return "inf" if val > 0 else "-inf"
        return val
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (list, tuple)):
        return [_json_clean(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _json_clean(v) for k, v in value.items()}
    return value


# -----------------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------------

def norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def unit(v: np.ndarray) -> np.ndarray:
    n = norm(v)
    if n <= 0.0:
        return np.zeros(3)
    return np.asarray(v, dtype=float) / n


def state_stress_vector_to_matrix(stress6: Iterable[float]) -> np.ndarray:
    s = np.asarray(stress6, dtype=float).ravel()
    if s.size != 6:
        return np.zeros((3, 3))
    # ExaDiS Python state order is xx, yy, zz, yz, xz, xy.
    return np.array(
        [
            [s[0], s[5], s[4]],
            [s[5], s[1], s[3]],
            [s[4], s[3], s[2]],
        ],
        dtype=float,
    )


def pbc_delta(cell: Any, r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
    """Return closest-image displacement when pyexadis exposes it.

    Falls back to direct displacement.  Rows using the fallback remain useful for
    local audits but should be labeled as not PBC-corrected if the simulation is
    periodic.
    """
    try:
        rpbc = np.asarray(cell.closest_image(Rref=np.asarray([r1]), R=np.asarray([r2])))[0]
        return rpbc - r1
    except Exception:
        return r2 - r1


def get_connectivity_from_segments(num_nodes: int, seg_nodeids: np.ndarray) -> list[list[int]]:
    conn: list[list[int]] = [[] for _ in range(num_nodes)]
    for sid, (n1, n2) in enumerate(seg_nodeids.astype(int)):
        if 0 <= n1 < num_nodes:
            conn[n1].append(sid)
        if 0 <= n2 < num_nodes:
            conn[n2].append(sid)
    return conn


def network_summary(N: Any, state: Dict[str, Any], label: str) -> Dict[str, Any]:
    G = N.get_disnet()
    try:
        nodes = G.get_nodes_data()
        segs = G.get_segs_data()
        n_nodes = int(len(nodes["positions"]))
        n_segs = int(len(segs["nodeids"]))
    except Exception:
        n_nodes = int(N.num_nodes())
        n_segs = int(N.num_segments())

    return {
        "record_type": "module_state",
        "label": label,
        "step": int(state.get("istep", -1)),
        "time_s": float(state.get("time", 0.0)),
        "dt_s": float(state.get("dt", 0.0)),
        "strain": float(state.get("strain", 0.0)),
        "stress_Pa": float(state.get("stress", 0.0)),
        "density": float(state.get("density", 0.0)),
        "Nnodes": n_nodes,
        "Nsegs": n_segs,
        "dEp": np.asarray(state.get("dEp", np.zeros(6))).ravel(),
        "dWp": np.asarray(state.get("dWp", np.zeros(3))).ravel(),
    }


# -----------------------------------------------------------------------------
# Mobility and force audit
# -----------------------------------------------------------------------------

def audit_mobility_force_records(N: Any, state: Dict[str, Any], recorder: EventAuditRecorder, phase: str) -> None:
    """Record segment/node force quantities exposed through Python bindings.

    Important convention: nodal forces from ExaDiS are total nodal forces.  The
    segment rows below project that total nodal force onto each connected arm.
    They are therefore not a true per-segment force decomposition unless a later
    native force-kernel patch provides per-segment force contributions.  This is
    explicitly recorded in ``force_decomposition``.
    """
    step = int(state.get("istep", -1))
    if not recorder.should_record(step):
        return

    G = N.get_disnet()
    nodes = G.get_nodes_data()
    segs = G.get_segs_data()
    tags = np.asarray(nodes["tags"], dtype=int)
    pos = np.asarray(nodes["positions"], dtype=float)
    constraints = np.asarray(nodes["constraints"], dtype=int).ravel()
    nodeids = np.asarray(segs["nodeids"], dtype=int)
    burgs = np.asarray(segs["burgers"], dtype=float)
    planes = np.asarray(segs["planes"], dtype=float)
    cell = G.cell

    forces = np.asarray(state.get("nodeforces", np.empty((0, 3))), dtype=float)
    force_tags = np.asarray(state.get("nodeforcetags", np.empty((0, 2))), dtype=int)
    vels = np.asarray(state.get("nodevels", np.empty((0, 3))), dtype=float)
    vel_tags = np.asarray(state.get("nodeveltags", np.empty((0, 2))), dtype=int)
    sigma = state_stress_vector_to_matrix(state.get("applied_stress", np.zeros(6)))

    force_by_tag = {tuple(t): forces[i] for i, t in enumerate(force_tags) if i < len(forces)}
    vel_by_tag = {tuple(t): vels[i] for i, t in enumerate(vel_tags) if i < len(vels)}

    for sid, (n1, n2) in enumerate(nodeids):
        if n1 < 0 or n2 < 0 or n1 >= len(pos) or n2 >= len(pos):
            continue
        r1 = pos[n1]
        r2 = pos[n2]
        dr = pbc_delta(cell, r1, r2)
        L = norm(dr)
        if L <= 0.0:
            continue
        line = unit(dr)
        burg = burgs[sid]
        bmag = norm(burg)
        plane = unit(planes[sid])
        glide_dir = unit(np.cross(plane, line))
        if norm(glide_dir) <= 0.0:
            # Degenerate or missing plane.  Use Burgers-projected direction as a
            # fallback and mark it.
            glide_dir = unit(burg - np.dot(burg, line) * line)
            glide_dir_source = "fallback_burgers_perp_line"
        else:
            glide_dir_source = "cross_plane_line"

        character_screw_abs = abs(float(np.dot(unit(burg), line))) if bmag > 0.0 else np.nan
        fpk_per_length = np.cross(sigma @ burg, line) if bmag > 0.0 else np.zeros(3)
        fpk_half_segment = 0.5 * L * fpk_per_length
        fpk_glide_N = float(np.dot(fpk_half_segment, glide_dir))
        tau_external_pk_Pa = fpk_glide_N / (bmag * 0.5 * L) if bmag > 0 and L > 0 else np.nan

        for local_node, endpoint in [(int(n1), "n1"), (int(n2), "n2")]:
            tag = tuple(tags[local_node])
            f_node = np.asarray(force_by_tag.get(tag, np.zeros(3)), dtype=float)
            v_node = np.asarray(vel_by_tag.get(tag, np.zeros(3)), dtype=float)
            f_glide_N = float(np.dot(f_node, glide_dir))
            tau_total_nodal_projected_Pa = f_glide_N / (bmag * 0.5 * L) if bmag > 0 and L > 0 else np.nan

            recorder.write({
                "record_type": "mobility_force_audit",
                "phase": phase,
                "step": step,
                "time_s": float(state.get("time", 0.0)),
                "dt_s": float(state.get("dt", 0.0)),
                "node_id": local_node,
                "node_tag": tag,
                "node_constraint": int(constraints[local_node]) if local_node < len(constraints) else None,
                "endpoint": endpoint,
                "seg_id": int(sid),
                "n1": int(n1),
                "n2": int(n2),
                "L_m": L,
                "half_L_m": 0.5 * L,
                "burg_x": float(burg[0]),
                "burg_y": float(burg[1]),
                "burg_z": float(burg[2]),
                "bmag_m": bmag,
                "plane_x": float(plane[0]),
                "plane_y": float(plane[1]),
                "plane_z": float(plane[2]),
                "line_x": float(line[0]),
                "line_y": float(line[1]),
                "line_z": float(line[2]),
                "character_screw_abs": character_screw_abs,
                "glide_dir_source": glide_dir_source,
                "glide_x": float(glide_dir[0]),
                "glide_y": float(glide_dir[1]),
                "glide_z": float(glide_dir[2]),
                "force_x_N": float(f_node[0]),
                "force_y_N": float(f_node[1]),
                "force_z_N": float(f_node[2]),
                "force_glide_from_total_nodal_force_N": f_glide_N,
                "tau_from_total_nodal_force_Pa": tau_total_nodal_projected_Pa,
                "external_PK_half_segment_glide_force_N": fpk_glide_N,
                "tau_external_PK_Pa": tau_external_pk_Pa,
                "vel_x": float(v_node[0]),
                "vel_y": float(v_node[1]),
                "vel_z": float(v_node[2]),
                "power_total_nodal_W": float(np.dot(f_node, v_node)),
                "power_glide_projected_W": float(f_glide_N * np.dot(v_node, glide_dir)),
                "force_decomposition": "total_nodal_force_projected_to_each_connected_arm_not_native_per_segment",
            })


# -----------------------------------------------------------------------------
# Candidate audits exposed from network geometry
# -----------------------------------------------------------------------------

def audit_topology_candidates(N: Any, state: Dict[str, Any], recorder: EventAuditRecorder, phase: str) -> None:
    step = int(state.get("istep", -1))
    if not recorder.should_record(step):
        return
    G = N.get_disnet()
    nodes = G.get_nodes_data()
    segs = G.get_segs_data()
    nodeids = np.asarray(segs["nodeids"], dtype=int)
    conn = get_connectivity_from_segments(len(nodes["positions"]), nodeids)
    for node_id, arms in enumerate(conn):
        if len(arms) >= 3:
            recorder.write({
                "record_type": "topology_candidate_audit",
                "phase": phase,
                "step": step,
                "time_s": float(state.get("time", 0.0)),
                "node_id": int(node_id),
                "degree": int(len(arms)),
                "connected_seg_ids": [int(x) for x in arms],
                "candidate_type": "multi_node_split_candidate_by_degree",
                "before_after_power_available": False,
                "native_status": "binding_candidate_only_requires_TopologyParallel_trial_power_hook",
            })


def audit_cross_slip_candidates(N: Any, state: Dict[str, Any], recorder: EventAuditRecorder, phase: str) -> None:
    step = int(state.get("istep", -1))
    if not recorder.should_record(step):
        return
    G = N.get_disnet()
    nodes = G.get_nodes_data()
    segs = G.get_segs_data()
    pos = np.asarray(nodes["positions"], dtype=float)
    constraints = np.asarray(nodes["constraints"], dtype=int).ravel()
    nodeids = np.asarray(segs["nodeids"], dtype=int)
    burgs = np.asarray(segs["burgers"], dtype=float)
    planes = np.asarray(segs["planes"], dtype=float)
    conn = get_connectivity_from_segments(len(pos), nodeids)
    cell = G.cell

    for node_id, arms in enumerate(conn):
        if len(arms) != 2:
            continue
        if int(constraints[node_id]) != 0:
            continue
        s1, s2 = arms
        # Require similar Burgers vectors and screw-like arms.
        b1 = burgs[s1]
        b2 = burgs[s2]
        b = unit(b1 + b2) if norm(b1 + b2) > 0 else unit(b1)
        chars = []
        lengths = []
        for sid in [s1, s2]:
            n1, n2 = nodeids[sid]
            other = n2 if n1 == node_id else n1
            if other < 0 or other >= len(pos):
                continue
            line = unit(pbc_delta(cell, pos[node_id], pos[other]))
            chars.append(abs(float(np.dot(b, line))))
            lengths.append(norm(pbc_delta(cell, pos[node_id], pos[other])))
        if not chars:
            continue
        screw_metric = min(chars)
        if screw_metric < 0.98:
            continue
        recorder.write({
            "record_type": "cross_slip_candidate_audit",
            "phase": phase,
            "step": step,
            "time_s": float(state.get("time", 0.0)),
            "node_id": int(node_id),
            "seg_ids": [int(s1), int(s2)],
            "screw_metric_min_abs_dot_b_line": screw_metric,
            "segment_lengths_m": lengths,
            "primary_plane_1": planes[s1],
            "primary_plane_2": planes[s2],
            "tau_primary_available": False,
            "tau_cross_available": False,
            "native_status": "binding_candidate_only_requires_CrossSlipParallel_force_projection_hook",
        })


def audit_module_delta(before: Dict[str, Any], after: Dict[str, Any], recorder: EventAuditRecorder, mechanism: str) -> None:
    if not recorder.enabled:
        return
    step = int(after.get("step", before.get("step", -1)))
    if not recorder.should_record(step):
        return
    recorder.write({
        "record_type": "module_delta_audit",
        "mechanism": mechanism,
        "step": step,
        "time_s": after.get("time_s", 0.0),
        "dt_s": after.get("dt_s", 0.0),
        "delta_Nnodes": int(after.get("Nnodes", 0)) - int(before.get("Nnodes", 0)),
        "delta_Nsegs": int(after.get("Nsegs", 0)) - int(before.get("Nsegs", 0)),
        "delta_strain": float(after.get("strain", 0.0)) - float(before.get("strain", 0.0)),
        "delta_stress_Pa": float(after.get("stress_Pa", 0.0)) - float(before.get("stress_Pa", 0.0)),
        "delta_density": float(after.get("density", 0.0)) - float(before.get("density", 0.0)),
        "collision_classification": (
            "topology_or_geometry_changed" if mechanism == "collision" and (
                int(after.get("Nnodes", 0)) != int(before.get("Nnodes", 0)) or
                int(after.get("Nsegs", 0)) != int(before.get("Nsegs", 0))
            ) else "no_count_change"
        ),
        "deterministic_geometry_only": mechanism in {"collision", "remesh"},
    })


# -----------------------------------------------------------------------------
# Audited Python stepping driver
# -----------------------------------------------------------------------------

class AuditedSimulateNetworkMixin:
    """Mixin overriding SimulateNetwork step hooks with audit records."""

    def __init__(self, *args, audit_recorder: Optional[EventAuditRecorder] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.audit_recorder = audit_recorder or EventAuditRecorder()

    def _record_state(self, N: Any, state: Dict[str, Any], label: str) -> Dict[str, Any]:
        row = network_summary(N, state, label)
        if self.audit_recorder.should_record(int(state.get("istep", -1))):
            self.audit_recorder.write(row)
        return row

    def step_integrate(self, N: Any, state: Dict[str, Any]):
        self._record_state(N, state, "before_integrate")
        self.save_old_nodes(N, state)

        self._record_state(N, state, "before_calforce")
        self.calforce.NodeForce(N, state)
        self._record_state(N, state, "after_calforce")
        audit_mobility_force_records(N, state, self.audit_recorder, "after_calforce_before_mobility")

        self._record_state(N, state, "before_mobility")
        self.mobility.Mobility(N, state)
        self._record_state(N, state, "after_mobility")
        audit_mobility_force_records(N, state, self.audit_recorder, "after_mobility")

        self._record_state(N, state, "before_time_integration")
        self.timeint.Update(N, state)
        self.plastic_strain(N, state)
        self._record_state(N, state, "after_time_integration")
        return state

    def step_topological_operations(self, N: Any, state: Dict[str, Any]):
        self._record_state(N, state, "before_topological_operations")
        audit_topology_candidates(N, state, self.audit_recorder, "before_topological_operations")
        audit_cross_slip_candidates(N, state, self.audit_recorder, "before_topological_operations")

        if self.cross_slip is not None:
            before = self._record_state(N, state, "before_cross_slip")
            self.cross_slip.Handle(N, state)
            after = self._record_state(N, state, "after_cross_slip")
            audit_module_delta(before, after, self.audit_recorder, "cross_slip")

        if self.collision is not None:
            before = self._record_state(N, state, "before_collision")
            self.collision.HandleCol(N, state)
            after = self._record_state(N, state, "after_collision")
            audit_module_delta(before, after, self.audit_recorder, "collision")

        if self.topology is not None:
            before = self._record_state(N, state, "before_topology")
            self.topology.Handle(N, state)
            after = self._record_state(N, state, "after_topology")
            audit_module_delta(before, after, self.audit_recorder, "topology")

        if self.remesh is not None:
            before = self._record_state(N, state, "before_remesh")
            self.remesh.Remesh(N, state)
            after = self._record_state(N, state, "after_remesh")
            audit_module_delta(before, after, self.audit_recorder, "remesh")

        self._record_state(N, state, "after_topological_operations")
        return state

    def step_update_response(self, N: Any, state: Dict[str, Any]):
        before = self._record_state(N, state, "before_update_response")
        out = super().step_update_response(N, state)
        after = self._record_state(N, state, "after_update_response")
        audit_module_delta(before, after, self.audit_recorder, "update_response")
        return out

    def run(self, N: Any, state: Dict[str, Any]):
        try:
            return super().run(N, state)
        finally:
            self.audit_recorder.close()


# The actual subclass is created after importing pyexadis_base in run_stock_case().


# -----------------------------------------------------------------------------
# Stock FCC example runner
# -----------------------------------------------------------------------------

def add_exadis_python_path(exadis_root: Path) -> None:
    p = exadis_root / "python"
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def build_stock_modules(state: Dict[str, Any], net: Any, use_cross_slip: bool = False):
    from pyexadis_base import CalForce, Collision, CrossSlip, MobilityLaw, Remesh, TimeIntegration, Topology

    calforce = CalForce(force_mode="SUBCYCLING_MODEL", state=state, Ngrid=64, cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(
        integrator="Subcycling",
        rgroups=[0.0, 100.0, 600.0, 1600.0],
        state=state,
        force=calforce,
        mobility=mobility,
    )
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(topology_mode="TopologyParallel", state=state, force=calforce, mobility=mobility)
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    cross_slip = CrossSlip(cross_slip_mode="ForceBasedParallel", state=state, force=calforce) if use_cross_slip else None
    return calforce, mobility, timeint, collision, topology, remesh, cross_slip


def run_stock_case(args: argparse.Namespace, audit_path: Optional[Path], outdir: Path) -> Dict[str, Any]:
    add_exadis_python_path(args.exadis_root)
    import pyexadis
    from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetwork

    class AuditedSimulateNetwork(AuditedSimulateNetworkMixin, SimulateNetwork):
        pass

    pyexadis.initialize()
    try:
        state: Dict[str, Any] = {
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

        G = ExaDisNet()
        G.read_paradis(str(args.exadis_data))
        net = DisNetManager(G)
        calforce, mobility, timeint, collision, topology, remesh, cross_slip = build_stock_modules(
            state, net, use_cross_slip=args.cross_slip
        )

        recorder = enable_audit(audit_path, stride=args.audit_stride) if audit_path else EventAuditRecorder()
        sim_cls = AuditedSimulateNetwork if audit_path else SimulateNetwork
        sim = sim_cls(
            calforce=calforce,
            mobility=mobility,
            timeint=timeint,
            collision=collision,
            topology=topology,
            remesh=remesh,
            cross_slip=cross_slip,
            vis=None,
            loading_mode="strain_rate",
            erate=args.strain_rate,
            edir=np.array([0.0, 0.0, 1.0]),
            max_strain=args.max_strain,
            burgmag=state["burgmag"],
            state=state,
            print_freq=args.print_freq,
            plot_freq=None,
            write_freq=args.write_freq,
            write_dir=str(outdir),
            audit_recorder=recorder,
            exadis_plastic_strain=True,
        ) if audit_path else sim_cls(
            calforce=calforce,
            mobility=mobility,
            timeint=timeint,
            collision=collision,
            topology=topology,
            remesh=remesh,
            cross_slip=cross_slip,
            vis=None,
            loading_mode="strain_rate",
            erate=args.strain_rate,
            edir=np.array([0.0, 0.0, 1.0]),
            max_strain=args.max_strain,
            burgmag=state["burgmag"],
            state=state,
            print_freq=args.print_freq,
            plot_freq=None,
            write_freq=args.write_freq,
            write_dir=str(outdir),
            exadis_plastic_strain=True,
        )
        final_state = sim.run(net, state)
        summary = network_summary(net, final_state, "final")
        summary["audit_rows"] = recorder.rows if audit_path else 0
        summary["audit_path"] = str(audit_path) if audit_path else None
        (outdir / "final_summary.json").write_text(json.dumps(_json_clean(summary), indent=2, sort_keys=True) + "\n")
        return summary
    finally:
        pyexadis.finalize()


def compare_summaries(disabled: Dict[str, Any], enabled: Dict[str, Any], tolerance: float) -> Dict[str, Any]:
    keys_float = ["strain", "stress_Pa", "density"]
    keys_int = ["Nnodes", "Nsegs"]
    diffs: Dict[str, Any] = {}
    passed = True
    for k in keys_float:
        d = float(disabled.get(k, 0.0))
        e = float(enabled.get(k, 0.0))
        diff = abs(e - d)
        scale = max(1.0, abs(d), abs(e))
        ok = diff <= tolerance * scale
        passed = passed and ok
        diffs[k] = {"disabled": d, "enabled": e, "abs_diff": diff, "rel_diff": diff / scale, "ok": ok}
    for k in keys_int:
        d = int(disabled.get(k, -1))
        e = int(enabled.get(k, -2))
        ok = d == e
        passed = passed and ok
        diffs[k] = {"disabled": d, "enabled": e, "ok": ok}
    return {"passed": passed, "checks": diffs}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run stock ExaDiS FCC strain-hardening with binding-level audit on/off.")
    p.add_argument("--exadis-root", type=Path, default=Path("core/exadis"))
    p.add_argument("--exadis-data", type=Path, default=Path("core/exadis/examples/22_fcc_Cu_15um_1e3/180chains_16.10e.data"))
    p.add_argument("--outdir", type=Path, default=Path("results/exadis_binding_event_audit"))
    p.add_argument("--max-strain", type=float, default=1.0e-5)
    p.add_argument("--strain-rate", type=float, default=1.0e3)
    p.add_argument("--burgers", type=float, default=2.55e-10)
    p.add_argument("--audit-stride", type=int, default=1)
    p.add_argument("--print-freq", type=int, default=1)
    p.add_argument("--write-freq", type=int, default=0)
    p.add_argument("--cross-slip", action="store_true", help="enable stock ForceBasedParallel cross-slip for candidate audit")
    p.add_argument("--tolerance", type=float, default=1.0e-10)
    p.add_argument("--mode", choices=["disabled", "enabled", "both"], default="both")
    return p


def main() -> int:
    args = build_parser().parse_args()
    args.exadis_root = args.exadis_root.resolve()
    args.exadis_data = args.exadis_data.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.exadis_root.exists():
        raise SystemExit(f"ExaDiS root not found: {args.exadis_root}")
    if not args.exadis_data.exists():
        raise SystemExit(f"ExaDiS FCC data file not found: {args.exadis_data}")

    summaries: Dict[str, Any] = {}
    if args.mode in {"disabled", "both"}:
        out = args.outdir / "audit_disabled"
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        summaries["disabled"] = run_stock_case(args, audit_path=None, outdir=out)

    if args.mode in {"enabled", "both"}:
        out = args.outdir / "audit_enabled"
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        summaries["enabled"] = run_stock_case(args, audit_path=out / "event_audit.jsonl", outdir=out)

    if args.mode == "both":
        comparison = compare_summaries(summaries["disabled"], summaries["enabled"], args.tolerance)
        summaries["comparison"] = comparison
        (args.outdir / "audit_invariance_summary.json").write_text(
            json.dumps(_json_clean(summaries), indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(_json_clean(summaries["comparison"]), indent=2, sort_keys=True))
        if not comparison["passed"]:
            raise SystemExit("Audit-enabled run is not invariant relative to audit-disabled run")
    else:
        (args.outdir / "audit_run_summary.json").write_text(json.dumps(_json_clean(summaries), indent=2, sort_keys=True) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
