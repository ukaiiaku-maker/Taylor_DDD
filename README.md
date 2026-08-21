# Taylor_DDD

OpenDiS-style reduced DDD drivers for studying Taylor hardening with Arrhenius Peierls and Taylor/forest-contact barriers.

This repository collects the most recent self-contained driver used in the reduced Taylor-hardening studies, analysis scripts for flow-stress and avalanche/burst diagnostics, unattended campaign launchers, and archived intermediate versions that document the development path.

## Native ExaDiS A1 model

The branch `arrhenius-exadis-strain-hardening` also carries an ordered native
patch stack pinned to LLNL ExaDiS SHA
`20ea2e82cdb919581c0611c338a6e46f6ad3f008`. It provides the native audit,
shared EXP-floor kernel, and gate-passed A1 Arrhenius Peierls mobility while
leaving stock ExaDiS as the default. Topology and cross slip remain fail-closed
because their native replacement gates did not pass.

See [the staged native report](docs/EXADIS_NATIVE_ARRHENIUS_STAGE_REPORT.md) for
the fit evidence, exact source blockers, validation metrics, and reproduction
commands.

## Current recommended driver

Use the v17 driver for new reduced-model work:

```bash
python3 clean_arrhenius_taylor_ddd_v17.py \
  --outdir results/test_transfer/bs_on/T1100_rho1e14 \
  --temperature-K 1100 \
  --strain-rate 4.5 \
  --target-strain 1e-4 \
  --dt 1e-8 \
  --forest-rho-m2 1e14 \
  --backstress-mobility on \
  --backstress-com-projection external_drive \
  --capture-mode swept_crossing \
  --snap-swept-capture-to-obstacle \
  --crossing-drive-mode force_work \
  --cross-force-scale-mode line_tension \
  --cross-force-scale-factor 0.25 \
  --max-free-dx-reduced 0.5 \
  --plastic-strain-source actual \
  --expfit-peierls-entropy-kB 0.0 \
  --expfit-cross-entropy-kB -9.25 \
  --expfit-cross-floor-frac 0.50 \
  --expfit-peierls-floor-frac 0.0 \
  --tau-local-cap-mode none \
  --tau-local-length-mode feed \
  --tau-local-L-eff-reduced 1.0
```

The key v17 setting is the uncapped local barrier-stress calculation:

```text
tau_local_cap_mode = none
tau_local_length_mode = feed
```

This avoids making the local stress cap into hidden physics. The former cap is retained only as a diagnostic reference in some outputs.

## Setup

```bash
python3 -m venv .venv-opendis
source .venv-opendis/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

Quick compile check:

```bash
python3 -m py_compile clean_arrhenius_taylor_ddd_v17.py
python3 -m py_compile analyze_v6_results.py
python3 -m py_compile analyze_depin_burst_statistics.py
python3 -m py_compile plot_depin_burst_ccdfs.py
```

## Repository layout

```text
clean_arrhenius_taylor_ddd_v17.py     current top-level driver
analyze_v6_results.py                 current top-level sweep analyzer
analyze_depin_burst_statistics.py      current top-level burst/avalanche analyzer
plot_depin_burst_ccdfs.py              current top-level CCDF plotting script
run_v17_systematic_rate_temp_campaign.sh unattended campaign launcher

drivers/                              organized driver copies
  archive/                            older v12-v16 drivers
analysis/                             analysis utilities
scripts/                              campaign and sanity-run launchers
patches/                              historical patch scripts/diffs
docs/                                 handoffs and method notes
examples/                             small example summaries/logs, not full simulations
```

Large simulation outputs should stay outside the repository under `results/` and are ignored by `.gitignore`.

## Avalanche / burst analysis

Run burst detection over a root:

```bash
python3 analyze_depin_burst_statistics.py \
  --root results/MY_ROOT \
  --cluster-gap-steps 1 \
  --active-plastic-ratio 1.0 \
  --stress-drop-threshold-MPa 0.0 \
  --n-boot 200 \
  --show-table
```

Plot CCDFs while removing the startup transient:

```bash
python3 plot_depin_burst_ccdfs.py \
  --root results/MY_ROOT \
  --eps-min 0.001 \
  --metrics depin plastic stress duration
```

See `docs/AVALANCHE_ANALYSIS_HANDOFF.md` for interpretation details.

## Systematic campaign

The long unattended campaign sweeps temperature, strain rate, and forest density with the fixed v17 uncapped barrier-stress model:

```bash
nohup ./run_v17_systematic_rate_temp_campaign.sh \
  > run_v17_systematic_rate_temp_campaign.nohup.log 2>&1 &
```

It is resumable: completed run directories with `run_summary.txt` are skipped.

## Development notes

The reduced model development established several implementation constraints that should carry into future 3-D junction-resolved DDD work:

1. Forest contacts/junctions must be mechanically load-bearing, not transparent event counters.
2. Pinned/contacted lines must still accumulate local force and backstress.
3. The local barrier stress must be derived from a physical local force/work balance.
4. Artificial local stress caps are diagnostic fuses, not production physics.
5. Plastic strain must be computed from actual swept line motion after contact creation, destruction, relaxation, and topology updates.
6. Burst/avalanche claims require CCDFs and null comparisons, not just event counts.

See `docs/TAYLOR_DDD_HANDOFF.md` for a fuller summary.
