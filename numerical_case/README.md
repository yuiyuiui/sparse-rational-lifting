# Fixed reproduction cases

Run from the repository root with a fresh output directory:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python numerical_case/benchmark.py --output numerical_case/raw/reproduction
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
