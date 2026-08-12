# Taylor DDD physics handoff

## Objective

The reduced OpenDiS-style calculations were built to test whether Taylor-like hardening can emerge from explicit forest-contact geometry and thermally activated contact destruction, rather than from an imposed analytic Taylor stress. The transferable lesson is that Taylor physics requires mechanically load-bearing contacts and a consistent local force/work balance for the Arrhenius barrier.

## Current useful formulation

The current useful reduced formulation is the v17-style model:

```text
plastic-strain-source = actual
tau_local_cap_mode = none
tau_local_length_mode = feed
crossing-drive-mode = force_work
cross-force-scale-mode = line_tension
cross-force-scale-factor = 0.25
expfit-cross-entropy-kB = -9.25
expfit-cross-floor-frac = 0.50
```

The exact entropy, floor, and force-scale values are reduced-model working values and should not be treated as universal material constants. The more important transfer is the mechanics: local barrier stress is derived from the pin/contact force over a physical feeding length.

## Lessons learned

### Contacts must be load-bearing

A forest contact that has formed must remain a constraint until it is explicitly destroyed, absorbed, annihilated, or topologically transformed. If obstacles are transparent, the model can count many events without producing Taylor hardening.

### Pinned lines must still feel force

Kinematic pinning should not zero the force on a pinned/contacted line. Projection schemes used to suppress artificial center-of-mass motion must not project away the backstress or line-tension force that loads the contact.

### Local barrier stress must come from force/work balance

The reduced model uses an equivalent local stress of the form

```text
tau_local = F_contact / (b L_eff)
```

where `L_eff` must be a physical length, such as the feeding segment length. A tiny numerical contact length creates singular stresses and forces the use of artificial caps.

### Stress caps are not physics

The useful v17 runs disable the active local stress cap. A cap can remain as a diagnostic reference, but production physics should not depend on clipping the local barrier stress.

### Plastic strain must use actual swept motion

A major correction was to compute plastic strain from actual unwrapped line motion after all capture, depinning, snapping, relaxation, and topology updates. Stress feedback should not be based only on the free-glide predictor.

### Event probabilities must be resolved

For each event class, track `R*dt`, the fraction of events near the barrier floor, and the event-rate distribution. If `R*dt` is large, the burst statistics can become timestep-controlled.

## Avalanche analysis interpretation

The previous v17 uncapped density sweep showed the most obvious burst-like behavior near the point where local pin stress first became substantial. This supported a cautious interpretation of intermittent depinning or pin-onset burstiness, not a definitive scale-free avalanche claim. Stronger avalanche claims require post-transient CCDFs, stress-drop distributions, null comparisons, and robustness to timestep and burst definitions.

## 3-D implementation implication

For future 3-D junction-resolved DDD, junctions should be real topological/mechanical objects. Each junction should store participating segment IDs, Burgers vectors, local force vector, effective reaction-coordinate length, barrier, rate, and destruction event history. The Arrhenius barrier should decide when a junction is destroyed, but the driving force must come from the actual 3-D force balance.
