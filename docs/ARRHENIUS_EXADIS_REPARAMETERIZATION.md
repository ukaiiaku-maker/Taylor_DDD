# Arrhenius / transition-state reparameterization of ExaDiS strain-hardening simulations

## Goal

The goal of this branch is to convert the full 3-D strain-hardening workflow used in ParaDiS/ExaDiS-style simulations into a transition-state-theory description.  The intended end state is not an athermal Taylor-strengthening law with an Arrhenius correction.  It is a DDD model in which every mechanism that permits dislocation motion or network evolution is represented as a thermally activated hazard with an explicit barrier, prefactor, stress mapping, and event increment.

The immediate target is the ExaDiS large-scale FCC strain-hardening example, because it already organizes the simulation around separable modules for force calculation, mobility, time integration, collision, topology, remeshing, and optional cross slip.  The stock ExaDiS example `examples/22_fcc_Cu_15um_1e3/example_fcc_Cu_15um_1e3.py` uses `CalForce`, `MobilityLaw`, `TimeIntegration`, `Collision`, `Topology`, `Remesh`, `CrossSlip`, and `SimulateNetworkPerf` with strain-rate loading.  That modular layout is the most useful insertion point for our Arrhenius laws.

## Existing model elements to replace or audit

The ExaDiS strain-hardening examples already contain rate-like ingredients.  The starting point should not be to discard the existing mechanics.  The starting point should be to identify the stress or force variable that each module uses and replace the deterministic response with a barrier-crossing hazard, or, where direct replacement is not yet exposed through Python, audit the equivalent barrier that the current model implies.

Mechanism-level mapping:

1. **Peierls glide / segment mobility**

   Current role: segment velocity from a mobility law such as `FCC_0`.

   TST replacement: signed forward-minus-reverse Arrhenius rate for glide increments,

   `v ∝ l_jump [R(+tau_eff) - R(-tau_eff)]`,

   with the barrier evaluated using the resolved glide stress or force-work equivalent.

2. **Forest / Taylor depinning**

   Current role: junction and forest constraints generate forces and topology changes; the network response produces an apparent Taylor hardening.

   TST replacement: explicit depinning hazard for load-bearing contacts or junctions.  The preferred stress convention is force-work based,

   `tau_eff = F_PK x_dagger / v_star`,

   rather than the diagnostic average stress `F_PK/(b L_eff)`.  This convention avoids double counting segment length while retaining the mechanical work that lowers the barrier.

3. **Binary junction formation and destruction**

   Current role: deterministic topology rules create and break reactions when geometric and force criteria are satisfied.

   TST replacement: junction zipping, unzipping, and destruction each receive an activation barrier.  The driving coordinate should be the local force-work change for the proposed reaction.  The deterministic geometry rule becomes a candidate generator and detailed-balance filter, not the kinetic law itself.

4. **Cross slip**

   Current role: optional force-based cross-slip mode.

   TST replacement: cross-slip is a competing event family with its own EXP-floor or KAA-like barrier, using the stress difference between primary and cross-slip planes as the local bias.

5. **Collision and annihilation**

   Current role: geometric collision handling and line annihilation.

   TST replacement: short-range annihilation can remain deterministic when it is purely geometric core overlap.  Any activated annihilation, climb-assisted escape, or junction bypass must become a hazard.  The distinction must be explicit in the code and outputs.

6. **Remeshing**

   Current role: numerical geometry maintenance.

   TST replacement: none.  Remeshing is numerical and should not become physics.  It must not change integrated hazards, accumulated plastic strain, or effective event counts.

## Common barrier family

The first implementation uses the same EXP-floor form used in the nanopillar and reduced Taylor studies,

```text
G(tau,T) = G0(T) [ f + (1-f) exp(-a (tau_eff/sigma_c)^n) ]
G0(T) = H - k_B T S
R = eta0 exp[-G(tau,T)/(k_B T)]
P(dt) = 1 - exp[-R dt]
```

The barrier input `tau_eff` must always be the local stress conjugate to the proposed event.  If the mechanism is naturally expressed by a force and activated displacement, the stress input is constructed from force work,

```text
tau_eff = F_event x_dagger / v_star.
```

This makes the barrier law portable across Peierls glide, forest depinning, junction unzipping, cross slip, and source activation.

## Stress-scale normalization across event families

The mechanism-specific stress scale `sigma_c` should not be treated as universal.  The same functional form can be used while allowing the stress scale to reflect the activated coordinate.  The practical normalization is:

```text
sigma_c,i = sigma_c,reference / phi_i
```

or equivalently,

```text
tau_eff,i = phi_i tau_local,i.
```

Only one convention should be active for a given mechanism.  The code must not multiply by `phi_i` twice.

## Parameter schema

Each event family should be described by a JSON-serializable object:

```json
{
  "mechanism": "forest_depinning",
  "barrier_family": "exp_floor",
  "H_eV": 0.50,
  "S_kB": -9.0,
  "sigma_c_GPa": 14.5,
  "floor_fraction": 0.20,
  "a": 6.65607,
  "n": 2.15276,
  "eta0_s": 1.0e12,
  "stress_convention": "force_work_tau_eff",
  "vstar_b3": 10.0,
  "x_dagger_rule": "vstar_over_b_squared",
  "event_increment_rule": "mechanism_specific"
}
```

A full strain-hardening model will contain several of these blocks:

```text
peierls_glide
forest_depinning
junction_zip
junction_unzip
cross_slip
activated_annihilation_or_climb
source_activation
```

## Development sequence

### Phase 1: audit-only wrapper

Use the stock ExaDiS strain-hardening example without modifying C++ kernels.  For each time step or saved state, compute the TST quantities that would correspond to the observed forces and rates:

- local candidate stress or force-work coordinate
- EXP-floor barrier
- per-candidate hazard
- aggregate hazard
- expected events per step
- actual plastic strain increment
- implied active mechanism fraction

This phase tests whether stock ExaDiS already samples stress scales consistent with our Arrhenius Taylor branch.

### Phase 2: Python-level candidate driver

Construct a Python driver that uses the ExaDiS network, force calculation, and strain-rate loading, but uses external TST hazards to decide which candidate events are permitted or suppressed.  This will be slower than native ExaDiS but is the correct way to verify the physics before editing C++.

### Phase 3: native ExaDiS mechanism hooks

Promote the Python-level law into native modules:

- `ArrheniusMobilityLaw`
- `ArrheniusCrossSlip`
- `ArrheniusTopology`
- optional `ArrheniusCollision` for activated short-range reactions

The native implementation should preserve the same parameter schema and produce the same audit columns as Phase 1.

### Phase 4: strain-hardening campaign

Run the existing FCC strain-hardening example under:

- stock ExaDiS mobility/topology
- Arrhenius Peierls only
- Arrhenius Peierls + forest depinning
- Arrhenius Peierls + forest depinning + junction topology
- full Arrhenius topology/cross-slip model

The acceptance criterion is not just matching a stress-strain curve.  The model must report the mechanism-resolved hazards and show that hardening emerges from the activated event population and evolving network geometry, not from hidden deterministic thresholds.

## Non-negotiable physics rules

1. No deterministic Taylor stress is allowed as the flow rule.
2. No deterministic depinning threshold is allowed as the release rule.
3. Stress caps may be retained only as diagnostic safety checks, not production physics.
4. Athermal geometric operations must be explicitly labeled numerical geometry maintenance or core-overlap annihilation.
5. Every activated mechanism must report its barrier, local stress convention, prefactor, and event probability.
6. The model must distinguish applied stress, local resolved stress, and force-work equivalent stress.

## Current status of this branch

This branch currently adds the shared TST law utilities and an ExaDiS adapter scaffold.  It is not yet a native ExaDiS replacement module.  The next code step is to run the audit-only wrapper against the stock ExaDiS FCC strain-hardening example and record whether the observed stress scales are compatible with the reduced Arrhenius Taylor laws.
