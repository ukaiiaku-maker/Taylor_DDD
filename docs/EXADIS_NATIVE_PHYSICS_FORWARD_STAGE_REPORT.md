# Native ExaDiS physics-forward Arrhenius stage report

## Outcome

The native path reaches A3: a directional EXP-floor FCC mobility, five
physically classified topology hazards, and geometry-first cross-slip hazards
execute inside ExaDiS. Stock behavior remains the default. Core-overlap
collision and numerical remeshing remain deterministic. The final one-thread
A3 audit-on/audit-off pair agrees bit for bit in every compared trajectory
field and complete network hash.

This is a calibrated mechanism model, not a claim that the listed barriers are
universal material constants. Real barriers can depend on dislocation and
junction character, temperature, local stress, chemistry, image forces, and
other state variables that the present loading set does not identify.

## Native event mechanics

The shared native EXP-floor kernel evaluates

```text
G = H [f + (1-f) exp(-a (max(tau_eff,0)/sigma_c)^n)] - kB T S
G_floor = f H - kB T S
G_used = max(G, G_floor, 0)
R = eta0 exp(-G_used/(kB T))
P = -expm1(-R dt)
```

It supports scalar Schmid/non-Schmid work and `W = sigma:A`,
`tau_eff = W/v*`. The configuration records activation enthalpy, activation
entropy, intrinsic activation volume, stress-concentration factor, site
multiplicity, and the native swept-area event strain increment. It verifies
`vstar = phi * intrinsic_activation_volume` before constructing a module.

Discrete topology and cross-slip candidates have stable event keys and sampled
exponential thresholds. Successive steps accumulate `R dt` until a threshold
is crossed. Independent admissible alternatives therefore compete in hazard
space, while an uncompleted sequential obstacle retains residence time and
cumulative hazard. Consumption or invalidation advances a renewal generation.
Large `R dt` uses a declared deterministic high-hazard limit rather than an
invalid one-Bernoulli-draw approximation.

Topology first constructs native split trials and separates numerical filters
from junction zip, junction unzip, destruction, reconfiguration, and
forest-depinning-like release candidates. Its conjugate work uses the native
before/after trial forces and displacement. The stock maximum-power decision is
retained only as an audit reference in Arrhenius mode.

Cross slip first constructs executable FCC plane-change or zipper geometry.
The native post-construction force-sign check remains a deterministic stability
condition; the stock force-magnitude preference is replaced by kinetics. The
Schmid term uses the total local native force difference between cross and
primary planes, so it includes applied, self, and interaction stress. External
normal, shear, and signed-normal terms supply the configured non-Schmid tensor
contribution. Every eligible decision audits the stress tensor, activation
tensor, `W`, barrier, hazard, residence state, threshold, and execution result.

## Parameter status

The physics-forward mobility is a reference-state-preserving transformation of
the trained A1 surrogate. Relative to the earlier `S = 0` mapping, a 0.04 eV
stress-independent enthalpy offset is paired with site multiplicity
`1.6749048555`; all directional rates remain exactly identical at 900 K for
every stress. The production values are `H = 0.0999544 eV`, `f = 0.529567`,
`sigma_c = 0.0395837 GPa`, `a = 36.509`, `n = 0.932572`, `jump = 0.231548 b`,
and `phi v* = 390.801 b^3`. The transition-state prefactor is
`eta0(T) = eta0(900 K) T/900 K`, with `eta0(900 K) = 1e12 s^-1` per site.

Topology uses class-specific mechanism priors (`H = 0.70–0.85 eV`,
`sigma_c = 0.20–0.25 GPa`, and `phi v* = 15–20 b^3`). Cross-slip plane change
and zipper propagation use calibrated surrogate barriers with
`H = 2.00/1.95 eV`, fixed high-drive floors `fH = 0.040/0.0375 eV`,
`sigma_c = 0.020 GPa`, `n = 1.2`, `phi v* = 15 b^3`, and site multiplicity
1.25. Raising the zero-drive enthalpy while preserving the high-drive floor and
the barrier near the reference driven candidates prevents undriven residence
hazards from masquerading as mechanically selected cross slip.

These topology and cross-slip values are trajectory-calibrated mechanism
surrogates. They have not been identified as atomistic material constants.

## Acceptance evidence

The final A3 run at 900 K and `2.2e-7` strain has exact audit invariance:
stress `30633.5359625 Pa`, density `1.153113920287591e12 m^-2`, 10,917 nodes,
11,791 segments, and network SHA-256
`14e67dd977dfd9722cd963249b7232f238d036f402b2c990efd517aea0ae6914`
match with auditing on and off.

The final A3 audit contains 2,160 topology rows, 125 kinetic topology rows, and
68 accepted topology events (`Rdt_max = 76.417`). It contains 337 cross-slip
rows, 128 geometrically stable kinetic candidates, and 15 accepted cross-slip
events (`Rdt_max = 89.557`). Both families exercise persistent residence state
and the deterministic high-hazard path.

Against the stock cross-slip reference, the A3 trajectory passes the 30% gate
with substantially smaller deviations:

| Quantity | Difference |
|---|---:|
| Stress maximum normalized by stock peak | 0.853% |
| Stress normalized RMSE | 0.606% |
| Final density | 0.000447% |
| Final node count | 0.119% |
| Final segment count | 0.0593% |

The physics-forward campaign covers all 20 combinations of
`T = 300, 500, 700, 900, 1100 K` and
`strain rate = 1e1, 1e2, 1e3, 1e4 s^-1` at a common `2.2e-7` comparison strain.
Every rate branch strengthens monotonically. Every temperature branch softens
within the declared 5% stochastic tolerance. The only local reversal is the
low-rate 700–900 K pair (20.10 to 21.00 kPa, 4.5%); the complete low-rate trend
still falls from 27.39 kPa at 300 K to 19.57 kPa at 1100 K. All networks are
sane, all audited hazards are finite, every Arrhenius decision is present, and
every large-`Rdt` row has cumulative-hazard or deterministic-limit metadata.

Initial-density factors 0.5, 1.0, and 2.0 are applied exactly and remain ordered
after evolution. At `2.05e-6` strain their stresses are 295.29, 294.52, and
292.10 kPa: a 1.1% athermal density plateau. This does **not** resolve
Taylor-like strengthening, so the validator reports
`taylor_like_strengthening_observed: false`; it is a documented remaining
physics limitation rather than a claimed pass hidden by tolerance.

## Scope and reproduction

The short reference trajectory observes forest-depinning-like release,
junction zip, destruction, and reconfiguration. Other classes require initial
states that actually present those geometries. The results do not validate
arbitrary crystals, loading tensors, or atomistically transferable barriers.

The ordered patch stack applies to LLNL ExaDiS SHA
`20ea2e82cdb919581c0611c338a6e46f6ad3f008` and exactly reproduces all 23 tested
native source files.

```bash
EXADIS_ROOT=core/exadis bash scripts/apply_exadis_native_audit_patch.sh

EVENT_AUDIT=ON BUILD_DIR=core/exadis/build-audit-llvm \
  PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
  bash scripts/build_exadis_native_audit.sh

OUTDIR=results/exadis_native_arrhenius_A3 \
  bash scripts/run_exadis_native_arrhenius_A3.sh

OUTROOT=results/exadis_native_arrhenius_physics_sweep \
  bash scripts/run_exadis_native_physics_sweep.sh

OUTROOT=results/exadis_native_arrhenius_density_sweep \
  bash scripts/run_exadis_native_density_sweep.sh
```
