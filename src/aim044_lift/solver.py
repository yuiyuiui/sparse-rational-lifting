"""Promised-input AIM 044 solver; exact integer arithmetic after rounding."""

from dataclasses import dataclass
from fractions import Fraction

from .reconstruct import recover_euclid, recover_flint


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


def _round_scaled(value, bits):
    # Exact round-to-nearest, with half ties toward positive infinity.
    return (2 * (value.numerator << bits) + value.denominator) // (2 * value.denominator)


def solve(n, entries, rhs, c, *, reconstruction="flint", trace=None):
    """Solve a promised AIM 044 input supplied as unique (row, col, rational) entries.

    Indices are zero based. The condition-number and nonsingularity promises
    are assumptions, not tested by this function. Basic representation and
    input-bit checks are performed. Default recovery is the fast FLINT backend;
    ``euclid`` is a slow reference mode with a different bit complexity.
    ``trace`` receives (iteration, residual_tuple, digit_tuple) for diagnostics.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    if not isinstance(c, int) or isinstance(c, bool) or c < 1:
        raise ValueError("c must be a positive integer")
    recover = {"flint": recover_flint, "euclid": recover_euclid}.get(reconstruction)
    if recover is None:
        raise ValueError("reconstruction must be flint or euclid")
    rhs = tuple(Fraction(value) for value in rhs)
    if len(rhs) != n or not any(rhs):
        raise ValueError("rhs must have length n and be nonzero")
    bit_bound = c * (2 * n - 1).bit_length()

    def check_bits(value):
        if max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > bit_bound:
            raise ValueError("Input exceeds the promised numerator/denominator bit bound")

    for value in rhs:
        check_bits(value)
    checked = []
    for row, col, value in entries:
        if (not isinstance(row, int) or not isinstance(col, int)
                or isinstance(row, bool) or isinstance(col, bool)
                or not 0 <= row < n or not 0 <= col < n):
            raise ValueError("Sparse index outside the matrix")
        value = Fraction(value)
        if not value:
            raise ValueError("The sparse nonzero list contains a zero")
        check_bits(value)
        checked.append((row, col, value))
    checked.sort(key=lambda item: (item[0], item[1]))
    if any(checked[i][:2] == checked[i-1][:2] for i in range(1, len(checked))):
        raise ValueError("Duplicate sparse index")
    if len(checked) < n:
        raise ValueError("A nonsingular n-by-n matrix needs at least n nonzeros")
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
    max_bits = 0
    row_sums = [0] * n
    for row, _, value in converted:
        row_sums[row] += abs(value)
    residual_bound = max(max(abs(value) for value in d), 1 + radix * max(row_sums))
    for iteration in range(lifts):
        max_bits = max(max_bits, max(abs(value).bit_length() for value in residual))
        digit = [value & mask for value in residual]
        if trace is not None:
            trace(iteration, tuple(residual), tuple(digit))
        product = [0] * n
        for row, col, value in converted:
            product[row] += value * digit[col]
        for i in range(n):
            blocks[i].extend(digit[i].to_bytes(s // 8, "little"))
            residual[i] = ((residual[i] - digit[i]) >> s) - product[i]
    max_bits = max(max_bits, max(abs(value).bit_length() for value in residual))
    modulus = 1 << (lifts * s)
    height = 1 << (3 * n * s)
    output = []
    for block in blocks:
        numerator, denominator = recover(int.from_bytes(block, "little"), modulus, height)
        output.append(Rational(numerator << s, denominator))
    return SolveResult(tuple(output), n, len(checked), len(converted), bit_bound,
                       h, s, lifts, lifts * len(converted), max_bits,
                       residual_bound.bit_length(), reconstruction)
