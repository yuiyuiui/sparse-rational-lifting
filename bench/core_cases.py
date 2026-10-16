"""Fixed nonnormal sparse systems with independent exact substitution checks."""
from fractions import Fraction
import hashlib
from aim044_lift import solve


def make(case):
    n = {'small': 8, 'medium': 64, 'large': 128}[case]
    entries = [(i, i, 1) for i in range(n)]
    entries += [(i, i + 1, Fraction(1, 4)) for i in range(n - 1)]
    reference = [Fraction(1)] * n
    for i in range(n - 2, -1, -1):
        reference[i] -= reference[i + 1] / 4

    def call():
        return solve(n, entries, [1] * n, 2)

    def validate(result):
        values = [v.as_fraction() for v in result.x]
        error2 = sum((x - y)**2 for x, y in zip(values, reference))
        norm2 = sum(x*x for x in reference)
        if error2 * 100 * n**4 > norm2 or result.max_residual_bits > result.residual_bit_bound:
            raise ArithmeticError('Independent exact error or residual-width check failed')
        digest = hashlib.sha256()
        for x in result.x:
            digest.update((hex(x.numerator) + '/' + hex(x.denominator)).encode())
        return {'passed': True, 'exact_squared_error_test': True,
                'output_sha256': digest.hexdigest()}

    return call, validate, {'n': n, 'nnz': len(entries), 'c': 2, 'family': 'upper bidiagonal'}
