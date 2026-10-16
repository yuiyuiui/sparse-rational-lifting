# Core optimization benchmark

The baseline is commit `9ba888176a85ec725a7c63d6da2a0af21c3373fa` (0.2.0).
The candidate shares FLINT reconstruction bounds across coordinates, fuses
residual-width accounting into updates, and releases consumed digit buffers.
All paired results have identical exact rational output digests and pass an
independent exact forward-error test against backward substitution.

| n | Baseline median ms | Candidate median ms | Baseline/candidate | Baseline Python peak bytes | Candidate Python peak bytes |
|---|---:|---:|---:|---:|---:|
| 8 | 0.228 | 0.172 | 1.32 | 10520 | 9076 |
| 64 | 12.199 | 11.007 | 1.11 | 300136 | 230572 |
| 128 | 64.586 | 55.222 | 1.17 | 1228048 | 940632 |

These are finite observations on the recorded machine, not asymptotic or
general superiority claims. Raw samples, cumulative profiles, effective settings,
source hashes and environment are in [the result](results/core_v0.2.1_final.json).
The lifting loop remains the principal cost at larger n after optimization.

## Reproduce

Install the candidate and its runtime dependency. Create a baseline checkout
from the public history, then run from this repository:

```sh
git worktree add --detach ../sparse-baseline 9ba888176a85ec725a7c63d6da2a0af21c3373fa
python bench/benchmark_core.py --baseline-source ../sparse-baseline/src --output bench/results/repeat.json
```

Use a new result filename. The fixed three cases use seven warm samples each,
with one numerical worker and one BLAS/OpenMP thread. Baseline/candidate order
alternates between cases. Input preparation/import and the first public call
are reported separately; validation and instrumentation are outside timing.
The benchmark supports no unbounded size sweep.

Python traced peak is one instrumented solve and excludes untraced native
allocations. RSS is the absolute fresh-worker high-water mark after timing,
including imports, inputs, outputs and validation, and can vary with allocator
and process startup. It is not a per-call native allocation measurement.
Profiling is a separate pass and does not contaminate timed samples.
