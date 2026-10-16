"""Promised-input AIM 044 solver; exact integer arithmetic after rounding."""

from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import chain

from .parameters import SolveConfiguration, select_contract
from .reconstruct import ensure_backend, recover_euclid, recover_many_flint


@dataclass(frozen=True)
class Rational:
    """A reduced rational pair returned by the reconstruction backend."""
    numerator: int
    denominator: int

    def as_fraction(self):
        """Convenience conversion for checks, outside the timed solver core."""
        return Fraction(self.numerator, self.denominator)


@dataclass(frozen=True)
class SolveResult:
    x: tuple[Rational, ...]
    n: int
    nnz: int
    converted_nnz: int
    input_bit_bound: int
    rounding_bits: int
    radix_bits: int
    lifts: int
    sparse_products: int
    max_residual_bits: int
    residual_bit_bound: int
    reconstruction: str
    config: SolveConfiguration | None = None

    def as_dict(self):
        """JSON-ready output and effective configuration; exact ratios use hex."""
        result = asdict(self)
        result["x"] = [{"numerator_hex": hex(x.numerator),
                        "denominator_hex": hex(x.denominator)} for x in self.x]
        result["config"] = None if self.config is None else self.config.as_dict()
        return result


def _round_scaled(value, bits):
    # Exact round-to-nearest, with half ties toward positive infinity.
    return (2 * (value.numerator << bits) + value.denominator) // (2 * value.denominator)


def solve(n, entries, rhs, c=None, *, condition_bound=None,
          relative_tolerance=None, reconstruction=None, trace=None):
    """Solve exact rational sparse data under an explicit mathematical promise.

    Existing ``solve(n, entries, rhs, c)`` calls retain their contract. Omit c
    and supply a justified condition_bound for deterministic automatic c
    selection from observed input bits, conditioning and optional accuracy.
    For a nonzero 1-by-1 scalar the condition bound is known exactly for free.
    The default backend is FLINT, with no silent reference fallback.

    See docs/API.md for types, all parameters, bounds, errors and costs.
    ``result.config`` records the effective contract and selection reasons.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    if trace is not None and not callable(trace):
        raise ValueError("trace must be callable or None")
    rhs = tuple(Fraction(value) for value in rhs)
    if len(rhs) != n or not any(rhs):
        raise ValueError("rhs must have length n and be nonzero")
    checked = []
    for row, col, value in entries:
        if (not isinstance(row, int) or not isinstance(col, int)
                or isinstance(row, bool) or isinstance(col, bool)
                or not 0 <= row < n or not 0 <= col < n):
            raise ValueError("Sparse index outside the matrix")
        value = Fraction(value)
        if not value:
            raise ValueError("The sparse nonzero list contains a zero")
        checked.append((row, col, value))
    checked.sort(key=lambda item: (item[0], item[1]))
    if any(checked[i][:2] == checked[i-1][:2] for i in range(1, len(checked))):
        raise ValueError("Duplicate sparse index")
    if len(checked) < n:
        raise ValueError("A nonsingular n-by-n matrix needs at least n nonzeros")
    observed_bits = max(max(abs(value.numerator).bit_length(), value.denominator.bit_length())
                        for value in chain(rhs, (value for _, _, value in checked)))
    config = select_contract(n, observed_bits, c, condition_bound,
                             relative_tolerance, reconstruction)
    c = config.c
    reconstruction = config.reconstruction
    # Fail before allocating digit blocks or running callbacks if native recovery is unavailable.
    ensure_backend(reconstruction)
    bit_bound = c * (2 * n - 1).bit_length()
    ell = (n - 1).bit_length()
    h = 4 + bit_bound + (2 * c + 1) * ell
    converted = [(row, col, _round_scaled(value, h)) for row, col, value in checked]
    converted = [(row, col, value) for row, col, value in converted if value]
    d = tuple(_round_scaled(value, h) for value in rhs)
    height_input = max([1] + [abs(value) for _, _, value in converted]
                       + [abs(value) for value in d])
    # A byte-aligned power-of-two radix permits one linear digit concatenation.
    s = 8 * (((2 * n * height_input - 1).bit_length() + 7) // 8)
    radix = 1 << s
    mask = radix - 1
    lifts = 6 * n + 2
    blocks = [bytearray() for _ in range(n)]
    residual = list(d)
    max_bits = max(abs(value).bit_length() for value in residual)
    digit_bytes = s // 8
    row_sums = [0] * n
    for row, _, value in converted:
        row_sums[row] += abs(value)
    residual_bound = max(max(abs(value) for value in d), 1 + radix * max(row_sums))
    for iteration in range(lifts):
        digit = [value & mask for value in residual]
        if trace is not None:
            trace(iteration, tuple(residual), tuple(digit))
        product = [0] * n
        for row, col, value in converted:
            product[row] += value * digit[col]
        for i in range(n):
            blocks[i].extend(digit[i].to_bytes(digit_bytes, "little"))
            residual[i] = ((residual[i] - digit[i]) >> s) - product[i]
            bits = abs(residual[i]).bit_length()
            if bits > max_bits:
                max_bits = bits
    modulus = 1 << (lifts * s)
    height = 1 << (3 * n * s)
    # Consume digit buffers in original coordinate order, releasing each after
    # conversion instead of retaining the complete digit array through recovery.
    blocks.reverse()
    residues = (int.from_bytes(blocks.pop(), "little") for _ in range(n))
    pairs = (recover_many_flint(residues, modulus, height) if reconstruction == "flint"
             else [recover_euclid(value, modulus, height) for value in residues])
    output = []
    for index, (numerator, denominator) in enumerate(pairs):
        output.append(Rational(numerator << s, denominator))
        pairs[index] = None  # Release the unscaled numerator immediately.
    return SolveResult(tuple(output), n, len(checked), len(converted), bit_bound,
                       h, s, lifts, lifts * len(converted), max_bits,
                       residual_bound.bit_length(), reconstruction, config)
