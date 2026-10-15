# Matched API/default-policy benchmark

NUMERICAL EVIDENCE: finite checks, no optimality/scaling or asymptotic speedup claim

All 16 input cases × 3 policies passed exact error and width checks, and returned identical rational outputs.

| Family | n | Default median ms | Explicit FLINT ms | Explicit Euclid ms |
| --- | ---: | ---: | ---: | ---: |
| conditioning_boundary | 4 | 0.0822 | 0.0816 | 0.0571 |
| conditioning_boundary | 8 | 0.2100 | 0.2135 | 0.1781 |
| conditioning_boundary | 16 | 0.6631 | 0.6929 | 0.6211 |
| conditioning_boundary | 32 | 2.3085 | 2.2629 | 3.0060 |
| small_rational_scale | 4 | 0.0812 | 0.0811 | 0.0659 |
| small_rational_scale | 8 | 0.2028 | 0.2028 | 0.2111 |
| small_rational_scale | 16 | 0.6372 | 0.6288 | 0.8539 |
| small_rational_scale | 32 | 2.2556 | 2.2563 | 3.7120 |
| indefinite_tridiagonal | 4 | 0.1084 | 0.1038 | 0.1347 |
| indefinite_tridiagonal | 8 | 0.3306 | 0.3047 | 0.7957 |
| indefinite_tridiagonal | 16 | 1.1116 | 1.0740 | 6.7288 |
| indefinite_tridiagonal | 32 | 4.2562 | 4.2270 | 64.3526 |
| nonnormal_bidiagonal | 4 | 0.0879 | 0.0885 | 0.0898 |
| nonnormal_bidiagonal | 8 | 0.2864 | 0.2597 | 0.5130 |
| nonnormal_bidiagonal | 16 | 0.8941 | 0.8085 | 3.8401 |
| nonnormal_bidiagonal | 32 | 3.2155 | 3.1517 | 36.0813 |

Three rotated warm repeats; full solve including normalization, selection, backend check, rounding, lifting and recovery. Reference and checks excluded.

Separate tracemalloc call per policy excludes native allocations; RSS is entire-process high-water including references/instrumentation, not per-mode memory.

Environment: {"python": "3.13.3", "system": "Darwin", "machine": "arm64", "package": "0.2.0", "python_flint": "0.8.0", "flint": "3.3.1", "logical_cpus": 14, "native_threads": 1, "openblas_threads": "1", "omp_threads": "1"}.
Whole-process peak RSS: 40189952 bytes.
Timing differences on small inputs include noise and do not establish a generally fastest policy.
