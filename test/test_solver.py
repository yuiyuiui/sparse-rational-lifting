"""Independent exact references and adversarial representation tests."""

import unittest
from fractions import Fraction
from math import gcd

from aim044_lift import solve
from aim044_lift.reconstruct import recover_euclid, recover_flint


def gaussian_reference(matrix, rhs):
    """Small exact Gaussian elimination, independent of modular lifting."""
    n = len(rhs)
    rows = [[Fraction(value) for value in row] + [Fraction(rhs[i])]
            for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [value / scale for value in rows[col]]
        for row in range(n):
            if row != col:
                factor = rows[row][col]
                rows[row] = [x - factor * y for x, y in zip(rows[row], rows[col])]
    return [row[-1] for row in rows]


def sparse(matrix):
    return [(i, j, Fraction(value)) for i, row in enumerate(matrix)
            for j, value in enumerate(row) if value]


class ReconstructionTests(unittest.TestCase):
    def test_signed_reduced_fractions_and_zero(self):
        modulus, height = 4096, 16
        for numerator in range(-9, 10):
            for denominator in range(1, 16, 2):
                if gcd(numerator, denominator) == 1:
                    residue = numerator * pow(denominator, -1, modulus) % modulus
                    with self.subTest(numerator=numerator, denominator=denominator):
                        pair = numerator, denominator
                        self.assertEqual(recover_flint(residue, modulus, height), pair)
                        self.assertEqual(recover_euclid(residue, modulus, height), pair)

    def test_large_half_gcd_case(self):
        modulus, height = 1 << 32768, 1 << 14000
        fraction = Fraction((1 << 13500) + 1234567, (1 << 13000) + 17)
        for sign in (-1, 1):
            residue = sign * fraction.numerator * pow(fraction.denominator, -1, modulus) % modulus
            self.assertEqual(recover_flint(residue, modulus, height),
                             (sign * fraction.numerator, fraction.denominator))

    def test_rejects_nonexistent_fraction(self):
        for recover in (recover_flint, recover_euclid):
            with self.assertRaises(ValueError):
                recover(17, 65536, 8)

    def test_parameter_domain(self):
        for parameters in ((0, 512, 16), (-1, 4096, 16), (4096, 4096, 16), (1, 4000, 16)):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                recover_flint(*parameters)


class SolverTests(unittest.TestCase):
    def assert_case(self, matrix, rhs, c):
        n = len(rhs)
        result = solve(n, sparse(matrix), rhs, c)
        reference = gaussian_reference(matrix, rhs)
        output = [value.as_fraction() for value in result.x]
        squared_error = sum((x-y)**2 for x, y in zip(output, reference))
        squared_norm = sum(x*x for x in reference)
        self.assertLessEqual(squared_error * 100 * n**(2*c), squared_norm)
        self.assertLessEqual(result.max_residual_bits, result.residual_bit_bound)
        self.assertEqual(result.lifts, 6*n+2)
        self.assertEqual(result.sparse_products, result.lifts*result.converted_nnz)
        self.assertEqual(result.reconstruction, "flint")
        return result

    def test_dimension_one(self):
        self.assert_case([[Fraction(1, 3)]], [Fraction(2, 3)], 2)

    def test_indefinite_and_nondyadic(self):
        self.assert_case([[-1, 0], [0, Fraction(1, 3)]], [1, 1], 2)

    def test_nonnormal(self):
        self.assert_case([[1, Fraction(1, 4), 0], [0, 1, Fraction(1, 4)], [0, 0, 1]],
                         [Fraction(1, 3), -1, Fraction(2, 5)], 2)

    def test_maximum_condition_promise(self):
        self.assert_case([[1, 0, 0, 0], [0, Fraction(1, 16), 0, 0],
                          [0, 0, Fraction(1, 4), 0], [0, 0, 0, 1]], [1, -1, 2, 1], 2)

    def test_small_matrix_and_rhs_scale(self):
        matrix = [[Fraction(int(i==j), 511) for j in range(4)] for i in range(4)]
        self.assert_case(matrix, [Fraction(1, denominator) for denominator in (509, 507, 505, 503)], 3)

    def test_non_power_of_two_dimension_and_large_c(self):
        matrix = [[Fraction(int(i==j), 32749) for j in range(3)] for i in range(3)]
        self.assert_case(matrix, [Fraction(1, 32719), Fraction(-1, 32717), Fraction(2, 32713)], 5)

    def test_zero_output_coordinate(self):
        result = self.assert_case([[1, Fraction(1, 3)], [0, 1]], [1, 0], 2)
        self.assertEqual(result.x[1].numerator, 0)
        self.assertEqual(result.x[1].denominator, 1)

    def test_working_system_exact_solution(self):
        matrix = [[2, Fraction(-1, 3)], [Fraction(1, 5), -2]]
        rhs = [Fraction(1, 7), Fraction(-2, 11)]
        result = self.assert_case(matrix, rhs, 3)
        scale, radix = 1 << result.rounding_bits, 1 << result.radix_bits
        def rounding(x):
            x = Fraction(x)
            return (2 * scale * x.numerator + x.denominator) // (2 * x.denominator)
        working = [[radix*rounding(matrix[i][j])+int(i==j) for j in range(2)] for i in range(2)]
        exact_z = gaussian_reference(working, [rounding(x) for x in rhs])
        self.assertEqual([x.as_fraction() for x in result.x], [radix*x for x in exact_z])

    def test_classical_reference_agrees(self):
        entries = sparse([[1, Fraction(-1, 3)], [0, 1]])
        fast = solve(2, entries, [1, -1], 2)
        slow = solve(2, entries, [1, -1], 2, reconstruction="euclid")
        self.assertEqual(fast.x, slow.x)

    def test_invalid_input_representation(self):
        cases = [(0, [], [1], 1), (1, [(0,0,1)], [1], 0),
                 (1, [(0,0,1)], [0], 2), (1, [(0,0,1)], [1,1], 2),
                 (1, [(0,0,0)], [1], 2), (1, [(1,0,1)], [1], 2),
                 (1, [(0,0,1),(0,0,2)], [1], 2), (1, [(0,0,8)], [1], 2)]
        for case in cases:
            with self.subTest(case=case), self.assertRaises(ValueError):
                solve(*case)


if __name__ == "__main__":
    unittest.main()
