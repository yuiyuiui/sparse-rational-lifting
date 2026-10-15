# Fixed reproduction cases

First install the package with `.venv/bin/python -m pip install .`, as shown
in the root README. All scripts use that installed package.

Run from the repository root with a fresh output directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python numerical_case/benchmark.py --output numerical_case/raw/reproduction
```

The script runs exactly four families at n=4,8,16,32, all with c=2, and three
warm timing repeats. It caps the whole run at 60 seconds. References use
independent exact rational Gaussian elimination; accuracy checks compare
squared norms using exact fractions, without floating-point tolerances. It
tests the stronger error bound n^-c/10 and the residual bit-width bound.

The diagonal family has condition number exactly n^2. The small-scale
diagonal family has condition number one. For the indefinite tridiagonal
family A=D+E, all diagonal magnitudes are two and ||E||_2<=1/2, giving
kappa_2(A)<=5/3. For the nonnormal bidiagonal family A=I+E and
||E||_2<=1/4, the same bound follows. All displayed numerator/denominator
sizes satisfy B=2 ceil(log2(2n)); the right-hand side is (-1)^i/(2i+1).

Reports include source hashes, environment, exact pass/fail inequalities,
lift/product counts, widths, and actual solver timings. Timings exclude exact
reference construction and validation. These finite cases validate software;
they do not prove the all-input theorem or its asymptotic bit complexity.
Raw output is ignored; the small selected report in `published/` records the
delivered version. No large simulation arrays are required by this method.

## Result figures

The gallery is generated from the selected 16-case report and one independently
checked trace of the existing conditioning-boundary case at `n=32`, `c=2`.
The trace checks each callback state/digit against the exact integer recurrence,
then validates the final solution against independent exact Gaussian
elimination and the recorded case. Its 195 points include the initial state
and a final residual derived from the last callback. It is diagnostic work,
outside the timed solver core.

| Figure | What it shows | Interpretation limit |
| --- | --- | --- |
| [Solution accuracy](figures/accuracy.svg) | Relative Euclidean error versus dimension, with the required `n⁻²/10` bound | Four sampled dimensions; connecting lines do not prove a rate |
| [Lifting widths](figures/lifting_width.svg) | Encoded residual width versus completed lifts, alongside accumulated digit capacity | Residual width is not forward error; capacity is not final numerator size |
| [Sparse work](figures/sparse_work.svg) | Actual scalar product counts for diagonal/bidiagonal/tridiagonal storage | Dense-entry curve is an operation-count reference, not a timed solver comparison or total bit cost |
| [Recorded timings](figures/timing.svg) | Median and observed min–max of three warm repeats for each case | Small cases on macOS arm64/Python 3.13.3/python-flint 0.8.0; ranges are not confidence intervals or evidence of general speedup |

The root README selects accuracy, lifting widths and sparse work to explain
the method's main features. These figures visualize validated finite data;
they do not certify the all-input theorem or asymptotic complexity.

### Reproduce the figures

Install the optional plotting dependency in your own environment. It is not
required to use the solver. Plotting reads selected data and does not rerun
the timing experiment:

```sh
.venv/bin/python -m pip install -r numerical_case/requirements-plot.txt
MPLCONFIGDIR=numerical_case/raw/mpl-config .venv/bin/python \
  numerical_case/plot_results.py --output numerical_case/raw/figures-reproduction --preview
```

Use a fresh output directory. SVGs are the primary repository figures;
`--preview` also creates PNGs for visual inspection. The generated
`FIGURE_MANIFEST.json` records data/script hashes, renderer versions,
validation checks and SVG hashes. With the recorded renderer environment,
SVG output is deterministic: no generation date or random SVG identifiers
are embedded. Layout can vary with other renderer versions.

To independently regenerate the selected trace, without rerunning the timing
cases, use another fresh directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python numerical_case/trace_case.py --output numerical_case/raw/trace-reproduction
MPLCONFIGDIR=numerical_case/raw/mpl-config .venv/bin/python \
  numerical_case/plot_results.py \
  --trace numerical_case/raw/trace-reproduction/lifting_trace.json \
  --output numerical_case/raw/figures-new-trace --preview
```

Both scripts reject changed solver/reference source hashes or failed
validation flags. The trace collector has a 60-second cap and runs one
fixed case. Selected inputs are `published/results.json` and
`published/lifting_trace.json`; the original timing data remain unchanged.
The figure manifest belongs with `figures/`. Raw reproduction and preview
outputs stay in the ignored `raw/` directory.

## API/default-policy comparisons (0.2.0)

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  numerical_case/benchmark_policies.py --output numerical_case/raw/policies-reproduction
```

This fixed 60-second-capped run compares `automatic_default`, `explicit_flint`
and `explicit_euclid` on the same four families and dimensions. For the default
path, the justified common input promise is `condition_bound=n**2`; explicit
paths supply c=2. Thus all modes have identical input, rounding and requested
error contracts. It is a deliberately matched comparison, not a search for
an optimal condition estimate or an accuracy-confounded speedup.

The script rotates the three execution orders over three warm repeats. Timings
include conversion/normalization, automatic selection or explicit checks,
backend preflight, rounding, lifting, reconstruction and result creation.
Reference construction, exact error checks and separate tracemalloc calls are
outside production times. Outputs must agree as exact rational pairs and meet
the n^-2/10 error bound against independent Gaussian elimination.

[policy_results.json](published/policy_results.json) stores actual effective
configuration/selection reasons, source hashes, Python/FLINT/package versions,
thread settings, exact acceptance flags, timings, CPU time, digit-storage
counts and memory measurements. [POLICY_SUMMARY.md](published/POLICY_SUMMARY.md)
provides the compact timing table. Tracemalloc excludes native allocations;
whole-process peak RSS includes references/instrumentation and is not per-mode
memory. No seed is used. These small fixed data do not establish general
optimality or asymptotic scaling.

All selected case/trace data and the four figures in the 0.2.0 inventory were
regenerated using the current package. The historical 0.1.0 inventories remain
unchanged at repository root and refer to their pinned earlier revision. For
current file verification use `releases/v0.2.0/SHA256SUMS` from repository root.
