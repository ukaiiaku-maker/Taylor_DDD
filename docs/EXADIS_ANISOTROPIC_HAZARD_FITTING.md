# ExaDiS anisotropic Arrhenius hazard fitting

## Purpose

This document describes the calibration step that comes after the stock ExaDiS audit and before any Arrhenius law is wired into native ExaDiS kinetics.

The goal is to fit equivalent transition-state hazards to the stock ExaDiS response for the three mechanism families requested next:

1. mobility / FCC_0 glide response,
2. cross-slip candidate selection,
3. collision or annihilation candidate selection.

The fitted laws are **not** yet native replacement physics. They are calibration targets. A later implementation must show that inserting these laws into ExaDiS preserves the intended mechanism-level event semantics and does not double-count length, line tension, or stress amplification.

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

Run the audit first.

```bash
PYTHON_BIN=/Users/sdillon/Taylor_DDD/.venv-opendis/bin/python \
EXADIS_ROOT=core/exadis \
MAX_STRAIN=1.0e-5 \
AUDIT_STRIDE=1 \
bash scripts/run_exadis_binding_event_audit.sh
```

Then fit the equivalent hazards.

```bash
PYTHON_BIN=/Users/sdillon/Taylor_DDD/.venv-opendis/bin/python \
ROOT=results/exadis_binding_event_audit/audit_enabled \
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

## Current limitations

The binding audit can project exposed total nodal forces onto connected arms, but it does not yet expose native force-kernel per-segment force decomposition. Native C++ instrumentation is still required to separate total nodal force into event-conjugate force components for forest depinning, junction unzipping, and cross-slip.

The cross-slip and collision fits become meaningful only after the audit exposes candidate-level accepted/rejected rows. If the candidate columns are absent, the fitter reports `no_data` or `no_acceptance_column` rather than inventing parameters.
