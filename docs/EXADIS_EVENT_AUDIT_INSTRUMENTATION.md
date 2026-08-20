# ExaDiS event-conjugate audit instrumentation

This document records the instrumentation layer added before any Arrhenius/TST
kinetic replacement is connected.

## What is implemented in this branch

The branch now contains a binding-level audit driver:

```text
exadis_audit/binding_event_audit.py
scripts/run_exadis_binding_event_audit.sh
```

The audit driver uses the stock ExaDiS Python module stack but runs through the
Python `SimulateNetwork` stepping path rather than `SimulateNetworkPerf`.  This
is deliberate.  The performance driver executes the stock C++ driver without
Python callbacks.  The Python stepping path exposes the sequence

```text
force -> mobility -> integration -> cross slip -> collision -> topology -> remesh -> response
```

and permits audit records before and after each module call without changing
accepted stock mechanics.

## EventAuditRecorder

`EventAuditRecorder` is disabled by default.  The binding-level API is

```python
from exadis_audit.binding_event_audit import enable_audit
recorder = enable_audit("event_audit.jsonl", stride=1)
```

A `None` path leaves the recorder disabled.  Audit rows are JSONL records.  This
is intentionally separate from the Arrhenius hazard law.  No barrier, rate, or
probability is used to accept or reject a stock ExaDiS event.

## Mobility-force audit

The audit records one row per segment endpoint after force calculation and after
mobility calculation.  Rows include

```text
node id and tag
segment id and endpoints
segment length
Burgers vector
plane normal
line direction
screw-character metric
nodal force
projected glide force
resolved stress inferred from total nodal force
external PK half-segment glide force
resolved stress inferred from applied stress
node velocity
nodal power and projected glide power
```

The important label is

```text
force_decomposition = total_nodal_force_projected_to_each_connected_arm_not_native_per_segment
```

This is a deliberate safety label.  ExaDiS exposes total nodal force through the
Python bindings.  The binding-level audit can project that total force onto each
connected arm, but this is not a native per-segment force decomposition.  A later
native force-kernel audit is still required before using these forces as
mechanism-specific `F_event` values in a TST barrier.

## Driver-level before/after module audit

The audited stepping class records state before and after

```text
calforce
mobility
time integration
cross slip
collision
topology
remesh
mechanical response update
```

State rows include step, time, dt, strain, stress, density, node count, segment
count, plastic strain increment, and plastic spin increment.  Module-delta rows
record changes in those quantities across each module.

## Topology candidate audit

The binding-level topology audit currently identifies multi-node split
candidates from network degree.  It records nodes with degree greater than or
equal to three and their connected segment ids.

It does not yet expose the native `TopologyParallel` trial power table.  The
source code location for the native hook is `src/topology_types/topology_parallel.h`,
where `SplitDisNet` trial configurations and force/mobility evaluations are used
to compare split candidates.  The next native hook must export the before/after
trial force, velocity, power, and accepted stock decision for each candidate.

## Cross-slip candidate audit

The binding-level cross-slip audit identifies screw-like two-arm nodes and logs
the current segment planes and screw metric.  It does not yet expose the native
candidate-plane force projections.

The native source location is `src/cross_slip_types/cross_slip_parallel.h`, where
ExaDiS identifies screw candidates, candidate planes, nodal force projection, and
stock cross-slip selection.  The next native hook must export `tau_primary`,
`tau_cross`, selected plane, rejected planes, and accepted stock event.

## Collision classification audit

The binding-level collision audit records before/after node and segment counts
around the collision module and labels count-changing events.  It cannot yet
distinguish core-overlap cleanup from activated collision candidates unless that
classification is exported from `src/collision_types/collision_retroactive.cpp`.

The native hook must classify each collision as

```text
deterministic core overlap
numerical cleanup
activated candidate reaction
```

Only the activated category is eligible for Arrhenius replacement.

## Validation run

Run the audit locally from the Taylor_DDD branch:

```bash
cd /Users/sdillon/Taylor_DDD
git switch arrhenius-exadis-strain-hardening

PYTHON_BIN=/Users/sdillon/Taylor_DDD/.venv-opendis/bin/python \
EXADIS_ROOT=core/exadis \
MAX_STRAIN=1.0e-5 \
AUDIT_STRIDE=1 \
bash scripts/run_exadis_binding_event_audit.sh
```

Expected outputs:

```text
results/exadis_binding_event_audit/audit_disabled/final_summary.json
results/exadis_binding_event_audit/audit_enabled/final_summary.json
results/exadis_binding_event_audit/audit_enabled/event_audit.jsonl
results/exadis_binding_event_audit/audit_invariance_summary.json
```

The validation criterion is that audit-enabled and audit-disabled runs give the
same final strain, stress, density, node count, and segment count within the
specified tolerance.

## What is not yet implemented

This branch still does not connect Arrhenius hazards to stock ExaDiS.  It also
does not yet provide native C++ device-side candidate tables for topology,
cross-slip, and collision.  Those are required before the Arrhenius model can be
connected without double counting force, line length, stress amplification, or
trial-configuration energetics.

The next commit should patch the native ExaDiS source tree so the hidden trial
candidate data are exported to a host-side audit buffer.  Only after that audit
is numerically invariant can `ArrheniusMobilityLaw`, `ArrheniusTopology`,
`ArrheniusCrossSlip`, and activated `ArrheniusCollision` be connected.
