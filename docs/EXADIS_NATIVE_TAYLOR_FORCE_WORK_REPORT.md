# Native ExaDiS Taylor force-work implementation report

## Outcome

The native implementation work is real and buildable, but the Taylor-model
mission has **not** passed its physics acceptance gates.

ExaDiS now has one shared `TaylorLineTensionInteraction` kernel. Native
topology split candidates classified as forest release, junction zip, junction
unzip, junction destruction, or junction reconfiguration call it. Native
cross-slip plane-change and zipper candidates call the same kernel. Core-overlap
collision and numerical cleanup remain deterministic because the current
collision path exposes no separately identified activated reaction class.

The interaction bias used by the EXP-floor law is exactly

```text
tau_eff = F_event * x_dagger / vstar
```

Native trial force work is preferred. A line-tension reconstruction using
`T = alpha * mu * b^2` and the two most event-aligned adjacent arms is the
fallback. The harmonic adjacent-arm length is recorded, but `L_eff/b` is not
multiplied into the force-work bias. `phi_eff = tau_eff/tau_app,resolved` and
`phi_geom = L_eff/b` are diagnostics only. The runtime configuration rejects a
fixed `stress_concentration_phi` in every D-D topology and cross-slip block.

The high-barrier 3-D density sweep does not show Taylor strengthening. The
required success statement therefore cannot be made.

## Native implementation

Patch `exadis_native_patches/0005-native-taylor-line-tension-force-work.patch`
adds the shared kernel and applies after patches 0001 through 0004 against
LLNL ExaDiS SHA `20ea2e82cdb919581c0611c338a6e46f6ad3f008`.

The implementation provides:

- a single C++ force-work evaluator in
  `src/arrhenius/taylor_line_tension_interaction.{h,cpp}`;
- native-trial-force and line-tension-reconstruction modes;
- harmonic adjacent-arm length plus minimum, sum, and geometric-amplification
  diagnostics;
- the five classified native topology hazards through one call site;
- cross-slip plane-change and zipper hazards through the same call site;
- persistent cumulative hazard and renewal-generation metadata inherited from
  the native A3 discrete-event implementation;
- barrier, rate, `Rdt`, force, work, length, effective-phi, swept-area, and
  event-selection audit fields;
- named single-glider regression, ExaDiS reference, and high-barrier scaling
  parameter sets;
- fail-closed parsing that prohibits a fixed D-D stress-concentration factor;
- density diagnostics and a fatal Taylor-slope validator which does not accept
  a plateau.

For selected events, the audit's swept area is calculated from the executable
native split/cross-slip trial coordinates. It is not an independent
post-topology reconstruction of every swept surface. The driver continues to
obtain its trajectory plastic strain from native ExaDiS network motion. This
distinction is retained explicitly rather than labeling the audit estimate as
a universally exact swept surface.

## Validation status

| Gate | Status | Evidence |
|---|---|---|
| Audit on/off invariance | Passed | Exact equality of all compared scalars, stress/density history hash, normalized restart hash, and complete network hash through the three-step A3 test |
| Native kernel versus v17 | Passed | Barrier absolute error 0 and rate relative error 0 for `F/Fc = 0...4`; fixed-force barrier unchanged for `L_eff/b = 100...10000` |
| Force-work double counting | Passed | `tau_eff = F*x_dagger/vstar` has zero relative identity error over 1,533 long-sweep interaction rows |
| Gate 1: dynamic single-glider ExaDiS reproduction | Incomplete | The pybind/native kernel regression is not a dynamic one-glider/fixed-forest ExaDiS density sweep |
| Gate 2: high-barrier 3-D density scaling | Failed | Stress is flat to weakly decreasing; no common-strain log slope is in 0.4–0.6 |
| Gate 3: full A3 campaign with Taylor kernel | Incomplete | A build, smoke run, and short exact-invariance run passed, but the required Taylor density response did not |
| Gate 4: peak search | Not run | Correctly withheld because Gate 2 failed |

The exact-invariance run ends at strain `3.2e-7`, stress
`44294.87957203192 Pa`, density `1.1529642947272441e12 m^-2`, 10,927 nodes,
11,795 segments, and network SHA-256
`3c987716570baa3cb37a26213ad1a938849ecc3575eb6ce2230451ac1661fe29`
with auditing both disabled and enabled.

The native smoke audit contains 409 Taylor interaction rows: 154
forest-depinning-like releases, 11 junction-zip candidates, 174 cross-slip
plane-change candidates, and 70 zipper-propagation candidates. A class can only
be validated when the starting network actually presents that geometry; this
smoke run does not exercise junction unzip, destruction, or reconfiguration.

## Failed density gate

The six-factor sweep (`0.25, 0.5, 1, 2, 4, 8`) reaches approximately
`2.05e-6` strain and fails at every common comparison strain. A longer bracket
sweep (`0.25, 1, 8`) reaches `1.08497e-5` strain and also fails:

| Strain | log(stress) / log(density) slope |
|---:|---:|
| `2e-7` | `+0.0209` |
| `1e-6` | `-0.00752` |
| `2e-6` | `-0.00987` |
| `5e-6` | `-0.0114` |
| `1e-5` | `-0.0175` |

At the final long-sweep state, stress decreases from 1.559 MPa at density
factor 0.25 to 1.464 MPa at factor 8. This is not Taylor strengthening.

The required density diagnostics distinguish the failure mode:

| Quantity | 0.25x | 1x | 8x |
|---|---:|---:|---:|
| total line density (m^-2) | `2.882e11` | `1.153e12` | `9.224e12` |
| mobile line density (m^-2) | `2.882e11` | `1.153e12` | `9.224e12` |
| forest/intersecting proxy (m^-2) | `5.742e10` | `3.949e11` | `5.529e12` |
| junction density (m^-3) | `5.322e16` | `4.261e17` | `9.560e18` |
| interaction-candidate density (m^-3) | `9.074e15` | `5.867e16` | `1.240e18` |
| accepted-interaction density (m^-3) | `0` | `0` | `0` |
| mean segment length (m) | `3.409e-7` | `3.297e-7` | `2.417e-7` |
| median `L_eff` (m) | `3.816e-7` | `3.407e-7` | `2.206e-7` |
| median `L_eff/b` | `1496` | `1336` | `865` |
| median `tau_eff` (Pa) | `1.530e8` | `1.642e8` | `4.873e8` |
| median `tau_eff/tau_app` on nonzero applied rows | `-28.7` | `91.3` | `575` |
| load-bearing-candidate density (m^-3) | `9.037e15` | `5.867e16` | `1.234e18` |
| transparent-candidate density (m^-3) | `0` | `2.963e14` | `0` |

The geometric contact population is therefore changing: candidate and
junction densities rise strongly and median `L_eff` falls. The shared kernel is
also demonstrably called. The failure is instead consistent with a
startup/Orowan transient dominated by the mobile network: final plastic strain
is only `6.75e-8`, `2.09e-7`, and `7.23e-7` for the three density factors, and
no high-barrier interaction candidate is accepted. The applied resolved stress
is also near zero for many audited candidates, so `phi_eff` is signed and
ill-conditioned there; it should not be used alone to infer clean
single-glider scaling.

The reported forest/intersecting density is a network-topology proxy (segments
touching degree-three nodes or junction-character segments), not persistent
contact-species metadata. Almost all line length is classified mobile. This is
a second reason not to interpret the present initial-network rescaling as the
required fixed-source single-glider forest-density experiment.

## Remaining work before a success claim

1. Construct the required dynamic ExaDiS one-glider/fixed-forest geometry and
   compare its stress, force work, event rates, and density scaling directly to
   v17. The existing exact kernel regression is necessary but not sufficient.
2. Add persistent D-D contact identities and participating-segment metadata so
   load-bearing forest contacts can be varied independently of mobile source
   density. The current topology-class inference is based on native node and
   Burgers/connectivity state.
3. Run to a post-startup flow regime with accepted interaction events while
   keeping `Rdt` resolved, then repeat all six density factors through at least
   `1e-5` common strain.
4. Reconstruct selected-event swept surfaces after topology execution if exact
   per-event plastic strain, rather than the executable native trial estimate,
   is required in the audit.
5. Only after the 0.4–0.6 high-barrier slope passes, rerun the full A3
   temperature/rate campaign and begin a transparency-peak search.

## Reproduction

```bash
EXADIS_ROOT=core/exadis bash scripts/apply_exadis_native_audit_patch.sh

EVENT_AUDIT=ON BUILD_DIR=core/exadis/build-audit-llvm \
  PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
  bash scripts/build_exadis_native_audit.sh

EVENT_AUDIT=OFF BUILD_DIR=core/exadis/build-noaudit-llvm \
  PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
  bash scripts/build_exadis_native_audit.sh

PYTHONPATH=core/exadis/build-audit-llvm/python \
  /opt/anaconda3/envs/opendis311/bin/python \
  validation/validate_native_taylor_single_glider_kernel.py \
  --output results/exadis_native_taylor_kernel_regression_v1/validation.json

OUTROOT=results/exadis_native_taylor_density_gate2 \
MAX_STRAIN=1e-5 \
INTERACTION_PARAMETER_SET=taylor_scaling_test_barrier \
OMP_NUM_THREADS=1 \
  bash scripts/run_exadis_native_density_sweep.sh
```

The density runner always writes the density diagnostics and double-counting
check even when the strict Taylor-slope validator returns failure.
