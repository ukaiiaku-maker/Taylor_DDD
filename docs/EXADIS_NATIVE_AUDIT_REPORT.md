# Native ExaDiS audit report

## Provenance

```text
Taylor_DDD base SHA: 68a51f30fa66184d45c26ff4c823f49247e274fa
ExaDiS source SHA:   20ea2e82cdb919581c0611c338a6e46f6ad3f008
Python:              /opt/anaconda3/envs/opendis311/bin/python
Kokkos:              OpenMP+Serial
OpenMP workers:      1
MPI:                 no
```

The source build used:

```bash
cmake -S core/exadis -B core/exadis/build-audit-llvm \
  -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ \
  -DEXADIS_ENABLE_EVENT_AUDIT=ON \
  -DEXADIS_PYTHON_BINDING=ON -DEXADIS_BUILD_EXAMPLES=OFF \
  -DEXADIS_FFT=ON -DFFTW_INC_DIR=/usr/local/include \
  -DFFTW_LIB_DIR=/usr/local/lib -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DPYTHON_EXECUTABLE=/opt/anaconda3/envs/opendis311/bin/python \
  -DPYEXADIS_OUTPUT_DIR="$PWD/core/exadis/build-audit-llvm/python"
cmake --build core/exadis/build-audit-llvm --target pyexadis -j4
```

A second build with `EXADIS_ENABLE_EVENT_AUDIT=OFF` completed and reported
`NATIVE_EVENT_AUDIT_COMPILED=false`; its `Driver` has no audit methods.

## Full stock invariance gate

Command:

```bash
PYTHON_BIN=/opt/anaconda3/envs/opendis311/bin/python \
BUILD_DIR=core/exadis/build-audit-llvm OMP_NUM_THREADS=1 \
OUTDIR=results/exadis_native_audit MAX_STRAIN=1e-5 CROSS_SLIP=0 \
bash scripts/run_exadis_native_audit.sh
```

| Quantity | Audit off | Audit on | Difference |
|---|---:|---:|---:|
| Strain | 1.021951400751762e-5 | 1.021951400751762e-5 | 0 |
| Stress (Pa) | 1,307,390.4867648063 | 1,307,390.4867648063 | 0 |
| Plastic strain | 1.1768853406074455e-6 | 1.1768853406074455e-6 | 0 |
| Density (m^-2) | 1.1485301879790017e12 | 1.1485301879790017e12 | 0 |
| Time (s) | 1.0219514007517618e-8 | 1.0219514007517618e-8 | 0 |
| Final dt (s) | 4.79219999055934e-10 | 4.79219999055934e-10 | 0 |
| Step | 23 | 23 | 0 |
| Nodes | 11,088 | 11,088 | 0 |
| Segments | 11,940 | 11,940 | 0 |

The complete node/segment array SHA-256 is identical:

```text
c8a04cdfa8b9ee27f09ed910592d54facee2e65fae1915b7a5e3224518617d52
```

The 24-row stress-strain-density files, initial configurations, and normalized
initial restarts also hash identically. The raw restart hash is not compared
because ExaDiS embeds wall-clock `date_and_time`.

Full enabled audit rows:

| Mechanism | Rows | Accepted | Rejected |
|---|---:|---:|---:|
| FCC_0 mobility | 529,620 | n/a | n/a |
| Topology split | 24,798 | 964 | 23,834 |
| Collision | 3,764 | 1,081 | 2,683 |
| Driver stage | 414 | n/a | n/a |

All collision rows in this run are explicitly classified deterministic
core-overlap geometry.

## Cross-slip candidate gate

The two-step `run_exadis_native_candidate_smoke.sh` run also passed exact
audit-off/on invariance. Its candidate counts are:

| Mechanism | Rows | Accepted | Rejected |
|---|---:|---:|---:|
| Cross slip | 359 | 59 | 300 |
| Topology split | 2,275 | 122 | 2,153 |
| Collision | 271 | 107 | 164 |
| FCC_0 mobility | 31,740 | n/a | n/a |

There are 27,931 nonzero mobility velocity rows; maximum absolute projected
velocity is 158.155 m/s and maximum absolute diagnostic `F/(bL)` stress is
0.979008503 GPa.

## Calibrated-surrogate fit gate

The native two-step trace was fitted at 900 K with
`exadis_calibration/anisotropic_hazard_fit.py`.

Mobility returned the reduced signed-work values
`eta0=1.1599481917312744e12 s^-1`, `v*=9.6994437434688 b^3`, and
`jump=0.9916127862599651 b`, with 91.5% sign agreement. Its RMS logarithmic
velocity error is 3.4157 decades, so it fails the replacement gate.

The cross-slip EXP-floor surrogate returned:

```text
H                 0.3948787692 eV
S                -9 kB (fixed initialization)
sigma_c           0.01000000049 GPa
floor fraction    0.01879654094
a                 19.99999969
n                 1.003428191
eta0              3.8869570973702516e13 s^-1
a_nn              2.999998927
a_mm              2.999999665
a_np             -2.999999996
median R dt       0.007997291128
95th percentile   0.5229589376
median P          0.007965397872
```

Several parameters are on optimization bounds, the data cover only one
temperature, and the fit does not identify site multiplicity or event strain.
It is therefore a calibrated mechanism surrogate, not a universal barrier.

Collision fitting reports `insufficient_candidate_labels` after excluding all
271 deterministic-geometry rows.

## Decision

```text
native_replacement_authorized: false
arrhenius_replacements_connected: false
delivery label: instrumented stock ExaDiS
```

No native replacement module was created. Independent future pathways must be
combined in hazard space; sequential obstacles require renewal/residence-time
treatment. A multi-temperature/state campaign with held-out trajectory checks
is required before any surrogate can enter the production simulation path.

The exact machine-readable evidence is under `results/exadis_native_audit`,
`results/exadis_native_candidate_smoke`, and
`results/exadis_native_hazard_fit`.
