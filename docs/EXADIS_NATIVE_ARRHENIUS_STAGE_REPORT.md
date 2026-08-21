# Native ExaDiS Arrhenius staged-conversion report

## Decision

The native conversion reached and passed A1: stock ExaDiS remains the default,
and a selectable native `MobilityFCC0Arrhenius` executes the trained directional
EXP-floor Peierls surrogate. A0 stock audit invariance and A1 audit invariance
are exact in one-thread mode. The long A1 stress, density, and network trajectory
gates pass.

The model is **not** a completed all-mechanism Arrhenius ExaDiS conversion.
Topology and cross slip failed their replacement gates and are deliberately not
connected. Collision remains deterministic core-overlap cleanup. Therefore A2,
A3, A4, and the density/temperature/rate production matrix were gated out.

## Provenance

```text
Taylor_DDD branch: arrhenius-exadis-strain-hardening
LLNL ExaDiS:       https://github.com/LLNL/ExaDiS.git
ExaDiS base SHA:   20ea2e82cdb919581c0611c338a6e46f6ad3f008
Python runtime:    /opt/anaconda3/envs/opendis311/bin/python
Fit runtime:       /Users/sdillon/Taylor_DDD/.venv-opendis/bin/python
Kokkos:            OpenMP + Serial
Validation workers: 1
MPI:               no
```

The two ordered patches are
`exadis_native_patches/0001-native-event-audit.patch` and
`exadis_native_patches/0002-native-arrhenius-exp-floor.patch`. They apply from
a clean checkout of the ExaDiS SHA above and reproduce the tested native files.

## Native implementation

The shared kernel in `src/arrhenius/arrhenius_exp_floor.{h,cpp}` implements

```text
G = H [f + (1-f) exp(-a (max(tau_eff,0)/sigma_c)^n)] - kB T S
G_floor = f H - kB T S
G_used = max(G, G_floor, 0)
R = eta0 exp(-G_used/(kB T))
P = -expm1(-R dt)
```

It supports scalar Schmid/non-Schmid effective stress and the tensor path
`W = sigma:A`, `tau_eff = W/v*`. Parameter validation is fail-closed and checks
finite, physical values.

`MobilityFCC0Arrhenius` is a native Kokkos-compatible nonlinear mobility; it
does not overwrite `MobilityFCC0`. It reconstructs the constrained nodal
generalized stress from projected nodal force and half the attached line
length, evaluates forward and reverse EXP-floor rates, and uses

```text
v = jump_b [R(+tau_eff,T) - R(-tau_eff,T)] direction.
```

The continuous mobility uses the native nonlinear subcycling integrator as a
mean-rate/multi-hit representation. `P` is audited as the probability of at
least one jump, but no Bernoulli draw is made for mobility. This is why saturated
`P` and `Rdt > 1` are permitted only with the declared adaptive mean-rate path.

Every directional evaluation records the temperature, applied stress tensor,
activation tensor, effective stress, activation work, barrier parameters,
`G`, `G_floor`, `G_used`, rate, `Rdt`, `P`, force, velocity, and network context.
Driver stage rows retain strain, plastic strain, density, node/segment counts,
stress, and plastic distortion increments.

Runtime selection is exposed through the native runner as:

```text
--arrhenius-mobility off|peierls|full
--arrhenius-topology off|on
--arrhenius-cross-slip off|on
--arrhenius-collision off|activated-only
--arrhenius-temperature-K
--arrhenius-eta0-default
--arrhenius-config
```

Stock is the default. `full`, topology, cross-slip, and activated-collision
requests fail closed unless the necessary fitted block has passed its gate.

## Mobility calibration gate

The campaign uses two complete initial configurations: the LLNL example state
and a second valid ParaDiS configuration after one deterministic native
evolution step. State 1 is training data; state 2 is held out wholesale across
all nodes, four rates, evolved snapshots, and five evaluated temperatures. The
stock `FCC_0` law has no temperature input, so its measured response is
explicitly evaluated as the target at 300, 500, 700, 900, and 1100 K; these are
not misrepresented as independent temperature-sensitive stock simulations.

The exact production config used for A1 was re-evaluated on the strengthened
646,800-row nodal table:

| Gate metric | Training | Held-out state 2 | Required | Result |
|---|---:|---:|---:|---|
| RMS log10 velocity error | 0.37713 | 0.41983 | <0.50 / <0.75 | pass |
| Sign accuracy | 1.00000 | 0.99995 | >0.97 held out | pass |
| Median velocity ratio | 1.25191 | 1.26340 | 0.5 to 2.0 | pass |
| Parameters on bounds | 0 | 0 | 0 | pass |
| Initial network states | 2 | 2 | at least 2 | pass |
| Temperatures evaluated | 5 | 5 | multiple | pass |

Its median and 95th-percentile `Rdt` are 66.61 and 733.78. The `Rdt` gate passes
only because the connected mobility is a continuous forward-minus-reverse mean
rate under native nonlinear subcycling, not a one-draw-per-step sampler.

The tested production parameters are:

```text
H                    0.04819349051199606 eV
S                   -0.15164376817344127 kB
sigma_c              0.039583654511514035 GPa
f                    0.024310824237702533
a                   36.50895532789194
n                    0.9325720424960138
jump                 0.23154770572269381 b
eta0                  1.0e12 s^-1 (fixed)
characteristic phiV* 390.8009690090542 b^3
non-Schmid terms      0 (not identified by this loading set)
```

The characteristic activation volume is the fitted EXP-floor derivative scale
`H(1-f)an/sigma_c`; the local differential volume also contains the stress-
dependent power/exponential factors. It is not a universal constant.

## A0 and A1 validation

The fresh post-kernel full A0 audit-on/audit-off pair is exact and also matches
the established full run: identical 24-row curves, complete network hash,
normalized restart hash, step/count fields, and scalar state with zero numerical
difference. The final A0 state is

```text
strain     1.021951400751762e-5
stress     1,307,390.4867648063 Pa
pstrain    1.1768853406074455e-6
density    1.1485301879790017e12 m^-2
nodes      11,088
segments   11,940
steps      23
network    c8a04cdfa8b9ee27f09ed910592d54facee2e65fae1915b7a5e3224518617d52
```

The final source builds with `EXADIS_ENABLE_EVENT_AUDIT=ON` and `OFF`. The
disabled binding reports `NATIVE_EVENT_AUDIT_COMPILED=false` and exposes no
driver audit methods; the enabled binding reports true and ran a final native
A1 smoke with 9,206 directional Arrhenius rows.

The two-step A1 run also passed exact audit-off/audit-on invariance. Its short
trajectory differs from stock by 1.66% in final stress, 0.021% in density,
0.037% in node count, and 0.119% in segment count.

The audit-enabled long A1 run reached strain `1.0849660268953602e-5` in 19
steps. Compared over common strain with A0:

| Quantity | Difference |
|---|---:|
| Normalized stress-curve RMSE | 0.989% |
| Maximum stress-curve deviation / stock peak | 1.737% |
| Maximum density-curve relative deviation | 0.196% |
| Final stress (different terminal strain) | 7.513% |
| Final density | 0.183% |
| Final plastic strain | 4.185% |
| Final node count | 0.054% |
| Final segment count | 0.117% |

All `Rdt`, `P`, and `G_used` values in the 385,190 A1 directional hazard rows
are finite. `G_used` ranges from
0.0135513 to 0.0599544 eV and no row enters an unintended transparent-barrier
state. The long native `Rdt` median and 95th percentile are 238.18 and 598.86;
`P` is consequently saturated for many rows, as expected for the explicitly
audited multi-hit mean-rate mobility.

## Mechanisms rejected by gate

### Topology/junction evolution

The strengthened fit has 17,775 native source rows and 810 positive rows. Its
whole-network held-out AUC is 0.99547, accuracy 0.98862, and class-balance error
13.17%, but it is not eligible:

- `f` and `n` are pinned to optimizer bounds;
- 95th-percentile `Rdt` is 4.0603, above the discrete-event limit;
- event classes are not resolved beyond generic multi-node split trials.

The exact source blocker is
`src/topology_types/topology_parallel.h`, particularly
`TopologyParallel::SplitMultiNode::operator()` and
`split_multi_nodes_parallel()`. ExaDiS constructs trial splits and selects a
maximum-power configuration per node. The available native records contain arm
masks, trial positions, and before/after force, velocity, and power, but the
algorithm does not expose physically distinct `junction_zip`, `junction_unzip`,
destruction, reconfiguration, or forest-depinning candidate objects. It also has
no persistent residence state for sequential obstacles. Hazard-gating the
generic rows would invent class semantics and violate the independent-pathway
versus renewal/residence-time rule.

### Cross slip

The strengthened fit has 3,419 source rows and 455 positives. Its held-out AUC
is 0.82557, accuracy 0.83621, and class-balance error 8.40%, but it is not
eligible:

- `S`, `sigma_c`, and `n` are pinned to optimizer bounds;
- 95th-percentile `Rdt` is 1.8030;
- rejected candidates do not all contain executable alternate geometry.

The exact blocker is `src/cross_slip_types/cross_slip_parallel.h`, in
`FindCrossSlipEvents::operator(team)`. The code computes a candidate plane and
force components for screw-like nodes, but fills executable event type and
branch-specific node positions only after deterministic force, pinning, and
geometry checks. A rejected `type == -1` row therefore lacks the complete
native operation that an Arrhenius acceptance would need. Switching such a row
would create geometry that stock ExaDiS never constructed.

Converting this safely requires a source-level extension that constructs and
audits valid alternate configurations for rejected physical candidates, plus
barrier data or labels that separate physical cross slip from geometry
admissibility. The audit alone cannot supply those missing states.

### Collision

`src/collision_types/collision_retroactive.cpp` reports the observed class as
deterministic core-overlap geometry. It remains stock cleanup. No physically
activated collision/reaction candidate class was observed, so
`activated-only` fails closed instead of wrapping numerical cleanup in a
fictitious barrier.

## Consequences and uncertainty

Independent future physical candidates must compete by summed hazards;
sequential obstacles require residence-time/renewal state. No additive stress,
max/min stress switch, or linear activation barrier was introduced.

The fitted EXP-floor mobility is a calibrated mechanism surrogate. Its negative
fitted entropy partly compensates for the temperature-independent stock target;
the data do not identify non-Schmid coefficients, site multiplicity independent
of the fixed prefactor, or a universal activation tensor. Real barriers may
depend on line and junction character, temperature, local chemistry, image
forces, and other state variables.

Because topology and cross slip are blocked at native candidate construction
and validation, running the requested 280-condition density/temperature/rate/
cross-slip matrix would falsely label A1 as a complete Arrhenius strain-
hardening model. That campaign is intentionally not run.

## Reproduction

```bash
# Apply both pinned native patches.
EXADIS_ROOT=core/exadis bash scripts/apply_exadis_native_audit_patch.sh

# Build audit-enabled and audit-disabled bindings.
EVENT_AUDIT=ON  BUILD_DIR=core/exadis/build-audit-llvm \
  PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
  bash scripts/build_exadis_native_audit.sh
EVENT_AUDIT=OFF BUILD_DIR=core/exadis/build-noaudit-llvm \
  PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
  bash scripts/build_exadis_native_audit.sh

# Recreate two-state fit data and all mechanism summaries.
bash scripts/run_exadis_native_calibration_campaign.sh

# Stock A0 and connected A1.
bash scripts/run_exadis_native_audit.sh
bash scripts/run_exadis_native_arrhenius_A1.sh
bash scripts/validate_exadis_native_arrhenius_A1.sh
```

Machine-readable generated evidence is under
`results/exadis_native_audit_post_kernel`,
`results/exadis_native_calibration_campaign/fit`,
`results/exadis_native_discrete_hazard_fit`, and
`results/exadis_native_arrhenius_A1_long`.
