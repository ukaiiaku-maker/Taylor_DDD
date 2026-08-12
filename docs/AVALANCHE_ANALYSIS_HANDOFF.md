# Avalanche / burst analysis handoff

## Purpose

The avalanche analysis was developed to distinguish ordinary independent Arrhenius depinning from correlated, burst-like depinning. The main diagnostic is not the total number of depinning events, but the distribution of depinning clusters.

A burst is defined as a cluster of active timesteps. A timestep is active if it contains depinning events, a plastic strain increment larger than the imposed strain increment, or a stress drop above threshold. The default settings were:

```text
cluster_gap_steps = 1
active_plastic_ratio = 1.0
stress_drop_threshold_MPa = 0.0
```

For each burst, the analysis records three main avalanche-size measures:

```text
S_N   = total depinning events in the burst
S_eps = total plastic strain released in the burst
S_tau = sum of negative stress drops in the burst
```

The CCDF,

```text
P(S' >= S)
```

is then plotted on log-log axes. A broad, slowly decaying CCDF is consistent with avalanche-like bursting, while rapid curvature/cutoff indicates finite stochastic bursts or ordinary depinning.

## Startup transient warning

The low-density simulations can produce misleading CCDFs if the initial yield/loading transient is included. A single startup burst can dominate the CCDF and appear avalanche-like despite negligible steady-state pinning.

For physical interpretation, use a post-transient cutoff such as:

```text
eps_total > 0.001
```

or

```text
eps_total > 0.0015
```

The CCDF plotting script supports this through:

```bash
--eps-min 0.001
```

## Null model

The analyzer compares clustering against a randomized null model. It keeps the same nonzero depinning counts but randomly redistributes them in time. This asks whether the observed burst clustering is stronger than expected from independent event counts with the same marginal distribution.

Important summary quantities include:

```text
null_largest_fraction_p_ge_obs
null_fano_p_ge_obs
largest_burst_fraction
p90_burst_depin
p99_burst_depin
loglik_power_minus_exp
```

## Recommended usage

Analyze an entire root:

```bash
python3 analyze_depin_burst_statistics.py \
  --root results/MY_ROOT \
  --cluster-gap-steps 1 \
  --active-plastic-ratio 1.0 \
  --stress-drop-threshold-MPa 0.0 \
  --n-boot 200 \
  --show-table
```

Plot all CCDFs in a root with a startup filter:

```bash
python3 plot_depin_burst_ccdfs.py \
  --root results/MY_ROOT \
  --eps-min 0.001 \
  --metrics depin plastic stress duration
```

## Interpretation used in the previous analysis

In the v17 uncapped density sweep at 1100 K and 0.45 s^-1, the response began to look most avalanche-like near rho ~ 1e15 m^-2, which was also near the point where the tail median local pin stress first increased substantially. Higher densities showed more frequent but smaller bursts, and the largest-burst fraction decreased. That suggested an onset/intermittency window rather than a fully established scale-free multi-hit avalanche regime.

The correct interpretation was therefore cautious: `bursty / intermittent depinning` is justified; `scale-free avalanche dynamics` requires stronger evidence.

For future 3-D junction simulations, preserve event-level records with step, time, event type, junction ID, local force, equivalent local stress, effective length, barrier, rate, R*dt, plastic increment, and before/after stress.
