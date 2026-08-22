# ExaDiS anisotropic Arrhenius hazard fitting

## Purpose

This document describes the calibration step that comes after the stock ExaDiS audit and before any Arrhenius law is wired into native ExaDiS kinetics.

The goal is to fit equivalent transition-state hazards to the stock ExaDiS response for the three mechanism families requested next:

1. mobility / FCC_0 glide response,
2. cross-slip candidate selection,
3. collision or annihilation candidate selection.

The fitted laws are **not** yet native replacement physics. They are calibrated
mechanism surrogates. A later implementation must show that inserting them into
ExaDiS preserves event semantics and does not double-count length, line tension,
stress amplification, or trial-configuration energetics.

## Theoretical structure

The attached Arrhenius-hazard paper writes each channel as a thermally activated strain-rate contribution with activation enthalpy, entropy, stress bias, site multiplicity, and strain increment. The useful decomposition for ExaDiS is

```text
rate_i = eta0_i exp[-G_i(stress,state,T)/(kB T)]
```

with the mechanical bias represented by an event-conjugate activation work. In the scalar reduced form this is controlled by the product `phi_i v*_i`. In the anisotropic form the scalar product is replaced by

```text
W_i = sigma : A_i
A_i = phi_i V*_i
```

For a slip system, the activation work can be decomposed into Schmid and non-Schmid terms,

```text
tau_eff,s = tau_s + a_nn sigma_nn,s + a_mm sigma_mm,s + a_np sigma_np,s
```

where `tau_s` is the resolved shear stress, `sigma_nn,s` is the normal stress on the slip plane, `sigma_mm,s` is the normal stress along the slip direction, and `sigma_np,s` is an optional secondary shear or non-planar component. If the audit does not expose non-glide components, the calibration falls back to Schmid-only coupling.

## Fitted models

### Mobility / FCC_0

The stock `FCC_0` mobility is overdamped. It is not a threshold law. The equivalent Arrhenius law is therefore fitted to the velocity response, not to a binary event.

The initial fitted form is the current reduced Peierls law,

```text
G_+(tau,T) = max(0, H_P - tau v*_P / eV) - kB T S_P
G_-(tau,T) = max(0, H_P + tau v*_P / eV) - kB T S_P
v_pred = jump * b * eta0 [exp(-G_+ / kBT) - exp(-G_- / kBT)]
```

The fitter estimates an equivalent prefactor, activation volume, jump distance, and available non-Schmid coefficients against the native velocity records.

### Cross slip

Cross-slip is treated as a competing candidate hazard once the audit exposes candidate screw segments and stock accepted/rejected labels. The fitted branch uses the corrected EXP-floor free energy,

```text
G(tau,T) = H [ f + (1-f) exp(-a (tau_eff/sigma_c)^n) ] - kB T S_kB
```

where the entropy is outside the floor function. Candidate geometry should come from ExaDiS. The fitted law should only replace kinetic acceptance, not candidate generation.

### Collision and annihilation

Collision rows labelled as deterministic core overlap or numerical cleanup are excluded from fitting. Only activated or ambiguous collision candidates are eligible for an Arrhenius hazard fit. This distinction is essential because core-overlap cleanup is geometric bookkeeping rather than an activated transition.

## Workflow

Build and run the native audit first.

```bash
PYTHON_BIN=/path/to/python bash scripts/build_exadis_native_audit.sh
PYTHON_BIN=/path/to/python OMP_NUM_THREADS=1 \
  bash scripts/run_exadis_native_candidate_smoke.sh
PYTHON_BIN=/path/to/python OMP_NUM_THREADS=1 MAX_STRAIN=1.0e-5 \
  bash scripts/run_exadis_native_audit.sh
```

Then fit the equivalent hazards.

```bash
PYTHON_BIN=/path/to/python \
ROOT=results/exadis_native_candidate_smoke/audit_enabled \
TEMPERATURE_K=900 \
STRAIN_RATE_S=1.0e3 \
bash scripts/fit_exadis_anisotropic_hazards.sh
```

Expected outputs:

```text
results/exadis_anisotropic_hazard_fit/anisotropic_hazard_fit_summary.json
results/exadis_anisotropic_hazard_fit/mobility_fit_observed_vs_predicted.csv
results/exadis_anisotropic_hazard_fit/cross_slip_fit_observed_vs_predicted.csv
results/exadis_anisotropic_hazard_fit/collision_fit_observed_vs_predicted.csv
```

Some files may be absent if the audit did not produce the relevant candidate rows.

## Acceptance criteria before native kinetic replacement

A fitted law is usable for native implementation only when all of the following are true:

1. Audit-enabled and audit-disabled stock ExaDiS runs have identical or near-identical stress, strain, plastic strain, density, node count, and segment count.
2. Mobility fit reproduces the stock velocity magnitude and sign over the force/stress range used by the strain-hardening example.
3. Cross-slip fit is trained on actual ExaDiS candidate rows, not only before/after network differences.
4. Collision fit excludes deterministic geometry cleanup and is trained only on activated or ambiguous candidate reactions.
5. The fitted stress coupling does not use both `phi` and an already-amplified force-work stress for the same mechanism.
6. The fitted law reports `R dt` and `P = 1-exp(-R dt)` for every candidate during validation.

## Current limitations and replacement gate

The native audit exposes actual accepted/rejected cross-slip, topology, and
collision candidates. FCC_0 still provides a total nodal force and velocity
projected onto each connected arm, not a unique per-arm barrier force.

The stock collision path observed here is deterministic core-overlap handling;
the fitter excludes it and reports `insufficient_candidate_labels`. A single
temperature trajectory also cannot identify activation enthalpy, activation
entropy, prefactor, and activation volume independently. Consequently fitted
EXP-floor parameters are reported with `replacement_eligible: false`. They must
not be promoted to universal constants or wired into native acceptance until a
multi-temperature/state campaign and held-out trajectory validation clear the
blockers recorded in the fit summary.

Independent event alternatives combine by summing hazards. Sequential barriers
require a renewal/residence-time model, so their waiting times are composed
rather than their hazards summed.

## Current native-audit fit result

The two-step candidate trace was fitted at the requested 900 K calibration
anchor. It did not clear the replacement gate:

- mobility used 27,793 nonzero rows but had `3.4157` decades RMS log-velocity
  error, despite `91.5%` sign agreement;
- cross slip used 359 accepted/rejected rows and reached `83.6%` thresholded
  classification accuracy, but `a`, `sigma_c`, and all three non-Schmid
  coefficients landed on optimization bounds;
- all 271 collision rows were deterministic geometry and were excluded, giving
  `status: insufficient_candidate_labels`;
- a single temperature cannot separately identify `H`, `S`, prefactor,
  effective activation volume, site multiplicity, and event strain increment.

Accordingly `results/exadis_native_hazard_fit/anisotropic_hazard_fit_summary.json`
records `native_replacement_authorized: false` and
`arrhenius_replacements_connected: false`. No native Arrhenius replacement was
created from these underdetermined surrogates.
