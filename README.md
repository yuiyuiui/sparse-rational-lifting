# AIM 044 sparse rational lifting

A Python implementation of the rational-input specialization of Jonathan A.
Kelner's identity-modulus sparse linear solver. The algorithmic mechanism and
fast bounded rational reconstruction are due to Kelner; this package supplies
explicit rational rounding, a native reconstruction binding and finite checks.

For the mathematical source, see Kelner, [Solving Linear Systems in
~O(mn log(κ/ε)) Bit Operations](https://arxiv.org/abs/2609.38101v1),
Theorem 2.3 and Section 3, and the
[pinned AIM 044 statement](https://github.com/MColbrook/AIM/blob/c9929805cece044705abaaebe1cebddac1ba11a6/problems/044-sparse-linear-solve-bit-complexity.md).
This is an implementation and rational-input specialization, with no claim of
algorithmic priority or upstream acceptance.

## Input and output

For every fixed positive integer `c`, the promised input is a nonsingular
`n × n` rational matrix supplied as unique nonzero `(row, column, value)`
entries, a nonzero rational vector `b`, and `κ₂(A) ≤ nᶜ`. Indices are zero
based. Each numerator and positive denominator has at most
`B = c ceil(log₂(2n))` bits. The solver returns reduced rational pairs for an
approximation with relative Euclidean solution error at most `n⁻ᶜ / 10`.
The mathematical construction uses deterministic `~O_c(nm)` bit operations,
where `m` is the number of nonzeros and hidden factors are logarithmic in `2n`.

The function checks dimensions, sparse indices, duplicate/zero entries, the
nonzero right-hand side and input bit bounds. Nonsingularity and the condition
bound are promises; it does not certify them. The operation bound includes the
solver's rounding, short sparse products, one digit assembly per coordinate,
fast bounded reconstruction and rational output. Exact Gaussian references,
optional tracing and diagnostic Fraction conversions are outside that bound.

## Install and run

Python 3.11 or newer and the pinned binary wheel `python-flint==0.8.0` are
required. The native binding requires one bundled FLINT library and the LP64
ABI on 64-bit macOS or Linux. Validation was performed on macOS arm64, Python
3.13.3, python-flint 0.8.0 / FLINT 3.3.1; Linux was not executed in this run.
The ABI is checked before native calls. Windows and custom non-wheel FLINT
builds are unsupported by this binding.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
PYTHONPATH=src .venv/bin/python -m aim044_lift numerical_case/example.json
```

Output numerators and denominators use hexadecimal strings, so long output
integers do not require decimal conversion. To use the API:

```python
from fractions import Fraction
from aim044_lift import solve

result = solve(2, [(0, 0, Fraction(-1)), (1, 1, Fraction(1, 3))], [1, 1], c=2)
for value in result.x:
    print(hex(value.numerator), hex(value.denominator))
```

## Tests and bounded reproduction

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python -m unittest discover -s test -v
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python numerical_case/benchmark.py --output numerical_case/raw/reproduction
```

Use a fresh reproduction output directory. The script caps execution at
60 seconds and uses four admissible families at `n=4,8,16,32`, `c=2`, with
three warm solver timing repeats. An independent exact Gaussian solve checks
the squared relative-error inequality without floating-point tolerances.
The tests also check signed reconstruction, rejected reconstruction domains,
zero coordinates, residual widths and a 32,768-bit modulus.

All 14 tests and 16 fixed cases passed in the recorded run. The selected small
reports in [numerical_case/published](numerical_case/published) record source
hashes, dependency versions, exact acceptance flags, widths and scalar timings.
See [case definitions and reproduction limits](numerical_case/README.md).
These are finite software checks, not an experimental proof of asymptotic
complexity or an all-input certificate. The diagonal conditioning endpoint,
small-scale diagonal, indefinite tridiagonal and nonnormal bidiagonal families
have analytically checked input promises.

## Reconstruction and representation

The default `flint` backend calls FLINT 3.3.1's
[`fmpq_reconstruct_fmpz_2`](https://github.com/flintlib/flint/blob/v3.3.1/src/fmpq/reconstruct_fmpz_2.c).
Its large-integer branch uses recursive half-GCD in
[`get_cfrac_helpers.c`](https://github.com/flintlib/flint/blob/v3.3.1/src/fmpq/get_cfrac_helpers.c).
Base-16 transfers avoid long decimal conversions. Digits are stored in byte
arrays and concatenated once. The radix exponent is rounded upward to a
multiple of eight, adding at most seven bits while preserving the proved
width and complexity bounds. Returned denominators are odd, so scaling the
numerators by the power-of-two radix does not introduce a long Python gcd.

`solve(..., reconstruction="euclid")` selects a separate classical Euclidean
reference. That mode has a different bit complexity and is not the fast backend
used in the theorem. The package does not vendor FLINT source or binaries;
FLINT and python-flint retain their own licenses.

## Attribution and AI use

Manuscript author: Kaiwen Jin. OpenAI Codex assisted with literature retrieval,
derivation, proof checking, implementation and finite-case validation.
No independent human review or formal verification is asserted. The lifting
method and fast reconstruction theorem are attributed to Kelner.
