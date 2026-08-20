# Native ExaDiS Arrhenius hook scaffold

This directory contains the first native-interface scaffold for moving the
Arrhenius/TST laws from the reduced Taylor drivers into ExaDiS/ParaDiS-style
strain-hardening simulations.

The current file is self-contained because this repository does not vendor the
ExaDiS source tree.  It is intended to be copied or mapped into the ExaDiS
source after the stock example audit identifies the local stress and force-work
coordinates exposed by the current mobility, topology, cross-slip, and collision
modules.

## Classes

- `ArrheniusMobilityLaw`
  - replacement target for Peierls glide mobility
  - uses signed forward-minus-reverse hazards
  - input coordinate is resolved glide stress

- `ArrheniusTopology`
  - replacement target for forest depinning and junction zip/unzip
  - geometry generates candidate reactions
  - hazards select kinetic events
  - force-work coordinate is `tau_eff = F_event x_dagger / v_star`

- `ArrheniusCrossSlip`
  - replacement target for force-based cross slip
  - crystallography generates candidate cross-slip systems
  - selection is a competing hazard

- `ArrheniusCollision`
  - only for activated collision/annihilation/reaction processes
  - core overlap and numerical collision cleanup remain deterministic and must be audited as `deterministic_geometry_only`

## Integration rule

Do not turn a stress cap, deterministic junction force threshold, or acritical
Taylor stress into hidden physics.  Every activated mechanism must output at
least:

```text
step,time,mechanism,tau_local,tau_eff,force,x_dagger,v_star,barrier_eV,rate_s,probability_dt,selected
```

Remeshing remains numerical geometry maintenance and should not receive an
Arrhenius hazard.
