# Callable API and parameter policy (0.2.0)

Install the distribution `sparse-rational-lifting` from a fixed repository
version; import `aim044_lift`. The public names are `solve`, `Rational`,
`SolveResult`, and `SolveConfiguration`. The legacy import and c-based calls
remain supported. No condition-number estimator or general singularity
certificate is included.

```python
solve(n, entries, rhs, c=None, *, condition_bound=None,
      relative_tolerance=None, reconstruction=None, trace=None)
```

## Problem inputs and mathematical promises

All quantities are dimensionless; this library assigns no physical units.
Matrix entries and the right-hand side must use mutually consistent units
before being passed to a linear solve.

| Parameter | Type and range | Default, role and recommendations |
| --- | --- | --- |
| `n` | Positive Python `int`, excluding `bool` | Required dimension of a square matrix and length of `rhs`; zero-based indices range from 0 through n-1. |
| `entries` | Iterable of `(row, column, value)` triples; integer indices and nonzero exact rational values | Required sparse nonzero list. No duplicate indices or explicitly stored zeros. Input order may vary; the solver sorts it. At least n nonzeros are necessary. Values retain the previous `Fraction(value)` conversion semantics; prefer `int` or `Fraction`. |
| `rhs` | Iterable of n rational values, not identically zero | Required right-hand side. Exact `Fraction` conversion, as for entries. |
| `c` | Positive Python `int`, excluding `bool`, or `None` | Default `None` selects c automatically using the rule below. A supplied c retains the legacy compound promise: nonsingular A, kappa_2(A)<=n**c, each reduced input numerator/denominator at most c*ceil(log2(2n)) bits, and conditional relative solution error <=n**(-c)/10. It is problem information, not an unchecked safe tuning default. |
| `condition_bound` | Finite `Fraction`-compatible value >=1, excluding `bool`, or `None` | Default `None`. Required for automatic selection when n>1; obtain a justified Euclidean condition upper bound from the application or an independent calculation. It is not inferred from general matrix entries and not verified. For a nonzero 1-by-1 scalar the exact condition number is 1, so the bound can be omitted; a larger supplied upper bound is safely tightened to 1. With explicit c, an optional supplied bound must fit n**c when n>1. |
| `relative_tolerance` | Finite `Fraction`-compatible value in (0,1], excluding `bool`, or `None` | Default `None` means the conditional n**(-c)/10 bound after c is selected. An explicit value is a requested relative Euclidean solution error. Automatic c increases until its guaranteed bound meets the request; explicit c fails if insufficient. For n=1 this policy only guarantees 1/10; smaller requests fail clearly even if a particular output happens to be more accurate. |

The nonsingularity and matrix-conditioning assumptions remain promises.
Representation checks are not their verification. An invalid matrix promise
can yield an output to which the advertised forward-error bound does not apply.
`result.config.matrix_promises_verified` stays false; a successful solve is
not certification.

Finite floats are converted to their **exact binary rational values**, not to
intended decimal fractions. Thus `0.1` generally consumes many more input bits
than `Fraction(1, 10)`. Use exact integers, `Fraction(p, q)` or decimal/rational
strings to represent intended rational data. This follows Python's
[`Fraction` conversion semantics](https://docs.python.org/3.13/library/fractions.html).
No rounding, truncation or denominator limiting is applied during conversion.

## Method settings

| Parameter | Supported settings and actual default | Recommendation and consequences |
| --- | --- | --- |
| `reconstruction` | `None` (default) resolves to `"flint"`; explicit `"flint"` or `"euclid"` | FLINT is the fast bounded-reconstruction implementation used in the bit-cost argument. Its pinned wheel/LP64 native ABI is checked before lifting. Missing or incompatible FLINT raises an error; no automatic fallback. Euclid is a slower exact reference for small comparisons and can run without loading the native library; it has a different bit complexity. Both recover the same exact working-system solution, rather than different accuracy levels. `"auto"` is not a backend name. |
| `trace` | `None` (default) or a callable `(iteration, residual_tuple, digit_tuple)` | Leave disabled for production timing. Called immediately before each of L lifts, starting at iteration zero; the final residual is not a callback state. It receives immutable integer tuples. Return values are ignored and exceptions propagate. Logging can add unbounded user work/storage, outside the solver's bit-cost claim. |

No random seed is used. Native FLINT uses one thread, set when this binding
initializes its shared library. This is process-global native state, not a
per-solve thread reservation; concurrent manipulation by other FLINT clients
has not been validated. There is no automatic reduced-precision mode.

## Exact automatic rule

All selection uses integer/rational comparisons, with no floating logarithms.
Let `u=ceil(log2(2n))`, and let `b` be the largest bit length of any reduced
input numerator or denominator (absolute numerator, positive denominator).

1. Set `c_input=max(1, ceil(b/u))`.
2. Set `c_condition` to the least positive exponent with
   `n**c_condition >= condition_bound`. For n=1 use the exact scalar condition
   number 1 and `c_condition=1`.
3. With an explicit requested tolerance, choose the least positive exponent
   with `1/(10*n**c_accuracy) <= relative_tolerance`. Otherwise use
   `c_accuracy=1`. The n=1 accuracy floor is checked separately.
4. Automatic c is `max(c_input, c_condition, c_accuracy)`. Explicit c is
   retained and validated against all supplied requirements.

No c=2 default is assumed. It is suitable only when the actual input/promise
and requested accuracy fit that exponent. A supplied upper bound is a claim
by the caller, not a verification by the selector.

With the selected c, the existing algorithm fixes `B=c*u`,
`h=4+B+(2c+1)*ceil(log2(n))`, rounds each input times `2**h` individually,
and sets `T=max(1, |rounded entries|, |rounded rhs|)`. The radix exponent is
the smallest multiple of eight at least `ceil(log2(2*n*T))`; `L=6*n+2`
lifts are used. Modulus and reconstruction height are `2**(L*s)` and
`2**(3*n*s)`, where s is the selected radix bit length. These are automatic
algorithm constants, not user knobs.

Larger c (from more input bits, a more conservative condition promise or a
stricter tolerance) generally increases h, s, digit storage and reconstruction
work. Rounding and reconstruction remain exact integer operations. The
forward-error theorem is conditional on the stated input promises, and its
`~O_c(nm)` notation treats c as fixed. No uniform bound independent of an
arbitrarily increasing c is advertised. Conversion/normalization is included
in measured solve time; the short-input theorem does not bound arbitrarily
large unreduced strings or caller-side rational construction for free. The automatic policy is not
claimed to be optimal. Excessively large promises/accuracy requests can use
large memory or time; there is no implicit resource truncation or precision
relaxation. Plan large runs explicitly.

## Outputs and errors

`SolveResult.x` is a tuple of reduced `Rational(numerator, denominator)` pairs
with positive odd denominators; numerators may be signed or zero. Its
`as_fraction()` convenience method may perform long-integer normalization
outside the production solver core.

`SolveResult.config` exposes the selected c, its source, supplied/effective
condition bound and source, requested tolerance, exact conditional error
bound, observed input bits, backend and selection reasons. It does not label
the matrix promises verified. Existing result fields expose n, input and
converted nonzeros, input bit bound, rounding/radix widths, lifts, scalar
sparse-product count, maximum residual width and its bound, and backend.
All sparse-product counts are scalar products, not matrix-vector calls.

`result.as_dict()` and `result.config.as_dict()` return JSON-ready forms.
Rational values use `numerator_hex`/`denominator_hex` strings without long
binary-to-decimal conversion. `json.dumps(dataclasses.asdict(result))` is not
the serialization API because configurations contain exact `Fraction`s.
The dataclasses are frozen containers; manually constructing them does not
validate a solve or a rational pair. Manual legacy `SolveResult` construction
can leave the new optional config field as `None`.

Invalid dimensions, indices, zero/duplicate entries, zero rhs, incompatible
c/bounds/tolerance, unsupported backend names and non-callable tracing produce
`ValueError`. Invalid rational construction can raise its ordinary conversion
error. Native dependency/ABI failures raise `RuntimeError` or `OSError`,
with no retry at weaker accuracy or different backend. Callback exceptions
propagate. Memory/time failures are not turned into successful approximate
outputs. The original c-based positional/keyword calls keep their arithmetic
and output; `reconstruction=None` now explicitly means the prior FLINT default.

## Installed examples and CLI

```python
from fractions import Fraction
from aim044_lift import solve

# This diagonal matrix has kappa_2=3, justified from its diagonal magnitudes.
entries = [(0, 0, -1), (1, 1, Fraction(1, 3))]
result = solve(2, entries, [1, 1], condition_bound=3)
assert result.config.c == 2
print(result.config.as_dict())

strict = solve(2, entries, [1, 1], condition_bound=3,
               relative_tolerance=Fraction(1, 1000))
reference_mode = solve(2, entries, [1, 1], condition_bound=3,
                       reconstruction="euclid")
legacy = solve(2, entries, [1, 1], 2)
assert result.x == reference_mode.x == legacy.x
```

The JSON CLI accepts the existing n/c/entries/b format, or n/entries/b plus
`condition_bound`, optional `relative_tolerance`, and optional
`reconstruction`. Rational problem parameters may be `[numerator, denominator]` integer pairs
(excluding booleans),
finite JSON numbers or rational strings. Entry values retain the four-element
`[row,column,numerator,denominator]` representation; b uses rational pairs.
CLI flags `--c`, `--condition-bound`, `--relative-tolerance` and
`--reconstruction` override corresponding JSON values. Absent fields/flags
are passed to the same API selector, without duplicate defaults. Missing
information and unsupported settings produce a nonzero exit and explanation.

```sh
.venv/bin/python -m aim044_lift numerical_case/example-auto.json
.venv/bin/python -m aim044_lift numerical_case/example-auto.json --relative-tolerance 1/1000
.venv/bin/python -m aim044_lift numerical_case/example-auto.json --reconstruction euclid
```

See [policy tests](../test/test_api.py), [exact reference tests](../test/test_solver.py),
and [matched benchmark/reproduction instructions](../numerical_case/README.md).
