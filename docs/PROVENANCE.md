# Sources and interpretation

This repository contains a selected archive from **Summer Research 2024** and a small analysis entry point prepared in October 2026. The original source files, two Excel datasets, and poster were copied without changes. File mappings and SHA-256 checksums are recorded in [source_manifest.json](source_manifest.json), using paths relative to the original research folder.

The poster, **Simulating Entropy Evolution**, is credited to **Junyu Hu and Calvin Stubbins**, Department of Physics & Astronomy, Franklin & Marshall College. Its references and attribution are preserved. No blanket license is assigned to the archived research materials.

## Selected experiments

| Example | Original saved dataset | What it demonstrates |
| --- | --- | --- |
| Random walk | `RW1/N=1k/l=50/RWdata300.xlsx` | 50 entropy trajectories, 1,001 recorded time points from 0 to 0.01 s. The folder identifies 1,000 particles. The source models isotropic fixed-length steps in a reflecting cube. |
| Ideal gas | `IG3/N=50k/data3.xlsx` | 50 entropy trajectories, 1,000 recorded time points from 0 to 0.0999 s. The accompanying `data3.txt` specifies 50,000 particles, a 0.1 s run, and a 0.0001 s step. The source models passage through a square opening in a central partition. |

Each workbook has `time`, 50 trial columns (`entropy`, `entropy_1` through `entropy_49`), `entropy_ave`, and `entropy_std`. The analysis recomputes the last two columns from all 50 trials and checks agreement. The shaded band is one **population standard deviation** (`ddof=0`) across realizations; it is not a confidence interval or standard error.

These are selected archived runs, not necessarily the datasets used in the final poster. In particular, the poster's random-walk illustration reports a different particle count. Current archived source defaults also differ from saved run settings, so running the source unchanged is not an exact reproduction procedure for these workbooks.

## Entropy convention

The historical code stores

$$S_{\mathrm{stored}} = 10^{22} k_B H,\qquad H=-\sum_i p_i\ln p_i,$$

with `k_B = 1.38e-23` in the source. The replay divides stored values by **0.138** and plots the dimensionless entropy **H** (natural logarithms, nats). For the gas example, the two probabilities are left/right occupancy, so the maximum is `ln(2)`. For the selected random-walk example, probabilities are spatial-bin occupancies, with the general finite-sample bound `H <= ln(1000)`. These two different coarse grainings should not be compared as equivalent absolute thermodynamic entropies.

Some historical filenames, comments, poster text, and figures use “Gibbs energy” or “Gibbs free energy.” The implemented quantity is coarse-grained entropy; it does not include an extensive factor of particle count. New documentation and plots use this definition explicitly.

## What the new entry point does

`analyze.py` is a **2026 presentation and verification helper**. It reads the unchanged Excel files, verifies the recorded statistics, converts their scale, and generates the figure and numerical summary. It does not rerun the dynamics, refit a model, or alter historical data. The recorded time values are retained: the original routines advance the particles before saving their first sample, although that sample is labeled zero.

## Historical code and limitations

The archive is useful for examining the research method and its implementation:

- `GibbsEnergy*.py` implements particle motion, occupancy counting, and entropy measurement.
- `Static*.py` runs repeated stochastic experiments and summarizes variability.
- `SolveEq.py` fits an autonomous polynomial model, `dS/dt = sum(a_i * (S-S_f)^i)`, by nesting SciPy numerical integration inside nonlinear least squares. Polynomial order is configurable; the poster discusses a second-order model, while some scripts and saved runs use third order.

The originals retain their research-stage limitations. They use global parameters, unseeded randomness, and slow Python particle loops. IG3 also has circular imports and approximate collision handling when several boundaries are crossed within a time step. Its aperture parameter is a half-width. RW1 allocates more histogram cells than its indexing requires, and its temperature parameter does not control the fixed step length. Fit routines do not consistently check solver convergence, and reported coefficient errors do not account for temporal correlation. Some time arrays and initial-value choices differ between simulation and fit code.

Accordingly, the quickstart uses saved data. The original simulators have **not** been validated or rerun as part of this portfolio preparation. Their larger configurations can require billions of particle updates. The source fit and root outputs should be treated as exploratory diagnostics, not as proof of a universal entropy law or validated equilibrium estimates. Running the archived code additionally requires SciPy and tqdm, and may require corrections to imports, parameters, and entry points.

The replay verifies stored data consistency; it cannot establish that the historical physical trajectories were computed correctly. Useful future work includes seeded vectorized dynamics, boundary-event tests, matched time coordinates, and comparisons of fitted models on held-out runs.
