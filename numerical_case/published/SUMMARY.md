# Reproduced cases

NUMERICAL EVIDENCE: finite implementation validation; not an asymptotic proof.

| Family | n | m | Lifts | Sparse products | Relative error | Median seconds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| conditioning_boundary | 4 | 4 | 26 | 104 | 2.594e-06 | 0.000097 |
| conditioning_boundary | 8 | 8 | 50 | 400 | 5.053e-08 | 0.000209 |
| conditioning_boundary | 16 | 16 | 98 | 1568 | 8.655e-10 | 0.000680 |
| conditioning_boundary | 32 | 32 | 194 | 6208 | 1.406e-11 | 0.002166 |
| small_rational_scale | 4 | 4 | 26 | 104 | 1.015e-06 | 0.000089 |
| small_rational_scale | 8 | 8 | 50 | 400 | 4.677e-07 | 0.000236 |
| small_rational_scale | 16 | 16 | 98 | 1568 | 7.443e-09 | 0.000594 |
| small_rational_scale | 32 | 32 | 194 | 6208 | 1.165e-10 | 0.002056 |
| indefinite_tridiagonal | 4 | 10 | 26 | 260 | 5.135e-07 | 0.000110 |
| indefinite_tridiagonal | 8 | 22 | 50 | 1100 | 6.063e-09 | 0.000298 |
| indefinite_tridiagonal | 16 | 46 | 98 | 4508 | 5.930e-11 | 0.001132 |
| indefinite_tridiagonal | 32 | 94 | 194 | 18236 | 7.396e-13 | 0.004006 |
| nonnormal_bidiagonal | 4 | 7 | 26 | 182 | 4.554e-07 | 0.000078 |
| nonnormal_bidiagonal | 8 | 15 | 50 | 750 | 6.404e-09 | 0.000238 |
| nonnormal_bidiagonal | 16 | 31 | 98 | 3038 | 5.186e-11 | 0.000738 |
| nonnormal_bidiagonal | 32 | 63 | 194 | 12222 | 7.100e-13 | 0.002914 |

All 16 cases passed the exact stronger forward-error inequality and the residual-width check.
Times describe this fixed Python/FLINT implementation on the recorded machine.
