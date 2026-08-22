# Native ExaDiS event-audit instrumentation

## Status

This branch now patches the high-performance C++ ExaDiS driver itself. The
earlier Python stepping audit remains available for historical comparison, but
it is not the evidence used by the native invariance gate.

The native patch is instrumentation only:

- `EXADIS_ENABLE_EVENT_AUDIT` defaults to `OFF` at compile time.
- an audit-enabled build still records nothing until `Driver.enable_audit()` is
  called;
- the recorder does not evaluate a barrier, draw a random number, change a
  candidate list, or replace a stock decision;
- `invariance_summary.json` explicitly records
  `arrhenius_replacements_connected: false`.

No trained Arrhenius hazard is connected in this branch.

## Reproducible source and patch

The patch is pinned to LLNL ExaDiS commit
`20ea2e82cdb919581c0611c338a6e46f6ad3f008` and is stored at
`exadis_native_patches/0001-native-event-audit.patch`.

```bash
git clone https://github.com/LLNL/ExaDiS core/exadis
git -C core/exadis checkout 20ea2e82cdb919581c0611c338a6e46f6ad3f008
bash scripts/apply_exadis_native_audit_patch.sh
bash scripts/build_exadis_native_audit.sh
```

The build helper accepts `PYTHON_BIN`, `BUILD_DIR`, `CMAKE_CXX_COMPILER`,
`FFTW_INC_DIR`, `FFTW_LIB_DIR`, and `BUILD_JOBS` overrides.

## Native hooks

### Driver stages

The C++ performance driver records before/after state around:

```text
force -> mobility -> integration -> plastic strain -> glide-plane reset
      -> optional cross slip -> collision -> topology -> remesh -> response
```

State rows include step, time, real time step, strain, stress, plastic strain,
density, node/segment counts, applied stress, and plastic strain/spin tensors.

### FCC_0 mobility

The stock subcycling driver has a zero-force outer mobility call, so the recorder
arms at integration entry and claims the first native FCC_0 evaluation inside
that integrator. It mirrors the device network once per audited outer step and
writes one row for each node/connected-arm pair. Rows contain physical
segment length, Burgers vector, plane, line direction, screw character, total
nodal force, native nodal velocity, their glide projections, power, and the
full applied stress tensor.

The force label is intentionally conservative:

```text
force_decomposition=total_nodal_force_projected_to_connected_arm
```

`tau_local_Pa = F_glide/(b L)` is a diagnostic arm projection. It is not claimed
to be a universal event-conjugate activation stress. Remaining internal
subcycling evaluations are not audited, which keeps the call structure intact
and prevents duplicate multi-gigabyte traces.

### TopologyParallel

Every native split trial is recorded after the stock trial-power table and
winner are known. Rows contain node and split IDs, arm-set mask, trial positions,
before/after/delta power in watts, and the stock accepted/rejected decision,
including neighbor-conflict and power-threshold rejection reasons.

### ForceBasedParallel cross slip

Every screw candidate is recorded after stock conflict resolution. Rows contain
both segment IDs, neighboring nodes, current and candidate planes, Burgers and
line directions, segment lengths, nodal force, primary/cross-plane force
projections, force threshold, diagnostic resolved stresses, proposed event type,
and the final stock accepted/rejected label.

### CollisionRetroactive

Segment and hinge collision candidates are recorded with closest-point
fractions, distance, segment IDs, rejection reason, and stock outcome. The
current native collision implementation is geometric core-overlap handling, so
every such row carries:

```text
deterministic_geometry_only = true
classification = core_overlap_geometry
not_eligible_for_arrhenius_fit
```

No activated collision law has been inferred from deterministic cleanup.

## Invariance gates

Run the short candidate-coverage gate:

```bash
PYTHON_BIN=/path/to/python \
BUILD_DIR=core/exadis/build-audit \
bash scripts/run_exadis_native_candidate_smoke.sh
```

This requires both accepted and rejected labels for `cross_slip`, `collision`,
and `topology_split`, and compares the audit-off/on final state and full network
digest.

Run the full stock strain-hardening gate:

```bash
PYTHON_BIN=/path/to/python \
BUILD_DIR=core/exadis/build-audit \
OMP_NUM_THREADS=1 \
MAX_STRAIN=1.0e-5 \
bash scripts/run_exadis_native_audit.sh
```

Outputs are:

```text
results/exadis_native_audit/audit_disabled/final_summary.json
results/exadis_native_audit/audit_enabled/final_summary.json
results/exadis_native_audit/audit_enabled/event_audit.jsonl
results/exadis_native_audit/invariance_summary.json
```

The exact reproducibility proof uses one OpenMP worker. Stock ExaDiS repeated
runs with four workers produced different node counts on this machine because
parallel candidate ordering is not deterministic; that behavior also occurs
with audit disabled. It must not be confused with an audit effect. Multi-thread
statistical equivalence is a separate validation problem.

The current recorder is restricted to serial/rank-0 output. Enabling it on a
nonzero MPI rank fails explicitly rather than interleaving JSONL records.

### Recorded validation evidence

Using ExaDiS `20ea2e82cdb919581c0611c338a6e46f6ad3f008`, the OpenMP+Serial
Kokkos backend, one OpenMP worker, and no MPI, the full `1.0e-5` gate passed at
step 23. Audit-off/on both ended at strain `1.021951400751762e-05`, stress
`1307390.4867648063 Pa`, plastic strain `1.1768853406074455e-06`, density
`1.1485301879790017e12 m^-2`, 11,088 nodes, and 11,940 segments. The complete
node/segment digest was identical:

```text
c8a04cdfa8b9ee27f09ed910592d54facee2e65fae1915b7a5e3224518617d52
```

The enabled full trace contains 529,620 FCC_0 mobility rows, 24,798 topology
split rows, 3,764 deterministic collision rows, and 414 driver stage rows.

The two-step cross-slip coverage gate also passed exactly. It contains 359
cross-slip candidates (59 accepted, 300 rejected), 2,275 topology candidates
(122 accepted, 2,153 rejected), and 271 collision candidates (107 accepted,
164 rejected). Of 31,740 mobility rows, 27,931 have nonzero projected velocity.

## Gate before kinetic replacement

A candidate fit is only a calibrated mechanism surrogate. Real barriers can
depend on character, stress, temperature, junction character, and other state.
Native replacement remains blocked until the relevant campaign identifies and
validates activation enthalpy, activation entropy, effective mechanical-work
tensor (or reduced `phi V*`), site multiplicity, and event strain increment.

For anisotropic work the fitting path uses `W = sigma:A`, resolving Schmid and
non-Schmid components from the recorded stress tensor and slip geometry.
Independent alternatives may be summed in hazard space. Sequential obstacles
must instead use renewal/residence-time treatment; they must not be represented
as parallel hazard channels.
