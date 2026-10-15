"""Installed-package policy regressions checked against exact Gaussian solves."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

from aim044_lift import solve
from test_solver import gaussian_reference, sparse

ROOT = Path(__file__).resolve().parents[1]


class ParameterPolicyTests(unittest.TestCase):
    def assert_accuracy(self, result, matrix, rhs, tolerance=None):
        reference = gaussian_reference(matrix, rhs)
        error = sum((value.as_fraction()-target)**2 for value, target in zip(result.x, reference))
        norm = sum(value**2 for value in reference)
        tolerance = result.config.relative_error_bound if tolerance is None else tolerance
        self.assertLessEqual(error, tolerance**2 * norm)
        self.assertLessEqual(result.max_residual_bits, result.residual_bit_bound)
        return result

    def test_minimal_automatic_call_and_legacy_agree(self):
        matrix = [[-1, 0], [0, Fraction(1, 3)]]
        result = self.assert_accuracy(solve(2, sparse(matrix), [1, 1], condition_bound=3), matrix, [1, 1])
        legacy = solve(2, sparse(matrix), [1, 1], 2)
        explicit = solve(2, sparse(matrix), [1, 1], c=2, reconstruction="flint")
        self.assertEqual(result.x, legacy.x)
        self.assertEqual(result.x, explicit.x)
        self.assertEqual(result.config.c, 2)
        self.assertEqual(result.config.relative_error_bound, Fraction(1, 40))
        self.assertEqual(result.config.reconstruction, "flint")
        self.assertFalse(result.config.matrix_promises_verified)

    def test_input_precision_controls_automatic_selection(self):
        matrix, rhs = [[1, 0], [0, 1]], [Fraction(1, 2**40+1), 1]
        result = self.assert_accuracy(solve(2, sparse(matrix), rhs, condition_bound=1), matrix, rhs)
        self.assertGreater(result.config.c, 1)
        with self.assertRaisesRegex(ValueError, "Input exceeds"):
            solve(2, sparse(matrix), rhs, c=result.config.c-1)

    def test_condition_endpoint_and_nonpower_upper_bound(self):
        matrix = [[1, 0, 0], [0, Fraction(1, 9), 0], [0, 0, 1]]
        outputs = []
        for bound in (9, Fraction(28, 3)):
            result = self.assert_accuracy(solve(3, sparse(matrix), [1, -1, 1], condition_bound=bound),
                                          matrix, [1, -1, 1])
            outputs.append(result.config.c)
        self.assertEqual(outputs, [2, 3])

    def test_requested_accuracy_is_met_and_changes_work(self):
        matrix, rhs = [[1, Fraction(1, 4)], [0, 1]], [1, -1]
        ordinary = solve(2, sparse(matrix), rhs, condition_bound=2)
        target = Fraction(1, 10000)
        strict = self.assert_accuracy(solve(2, sparse(matrix), rhs, condition_bound=2,
                                           relative_tolerance=target), matrix, rhs, target)
        self.assertLessEqual(strict.config.relative_error_bound, target)
        self.assertGreater(strict.rounding_bits, ordinary.rounding_bits)
        self.assertGreater(strict.radix_bits, ordinary.radix_bits)

    def test_missing_condition_information_is_actionable(self):
        with self.assertRaisesRegex(ValueError, "Provide condition_bound"):
            solve(2, sparse([[1, 0], [0, 1]]), [1, 1])

    def test_scalar_automatic_boundary_and_accuracy_floor(self):
        matrix, rhs = [[Fraction(1, 3)]], [Fraction(2, 3)]
        for options in ({}, {"condition_bound": 17}, {"relative_tolerance": Fraction(1, 10)}):
            with self.subTest(options=options):
                result = self.assert_accuracy(solve(1, sparse(matrix), rhs, **options), matrix, rhs)
                self.assertEqual(result.config.condition_bound, 1)
        with self.assertRaisesRegex(ValueError, "For n=1"):
            solve(1, sparse(matrix), rhs, relative_tolerance=Fraction(1, 100))

    def test_invalid_and_incompatible_options_fail(self):
        matrix, rhs = [[1, 0], [0, 1]], [1, 1]
        cases = [{"c": value} for value in (False, 0, -1, 1.5)]
        cases += [{"c": 2, "condition_bound": value}
                  for value in (True, 0, Fraction(9, 10), float("inf"), float("nan"), [])]
        cases += [{"c": 2, "relative_tolerance": value}
                  for value in (True, 0, -1, 2, "bad", float("nan"))]
        cases += [{"c": 2, "reconstruction": value} for value in (False, [], "auto", "missing")]
        cases += [{"c": 2, "trace": False}, {"c": 1, "condition_bound": 3},
                  {"c": 2, "relative_tolerance": Fraction(1, 10000)}]
        for options in cases:
            with self.subTest(options=options), self.assertRaises(ValueError):
                solve(2, sparse(matrix), rhs, **options)

    def test_unavailable_backend_fails_before_trace_and_never_falls_back(self):
        trace = []
        with patch("aim044_lift.reconstruct._library", side_effect=RuntimeError("Unavailable FLINT")):
            with self.assertRaisesRegex(RuntimeError, "Unavailable FLINT"):
                solve(2, sparse([[1, 0], [0, 1]]), [1, 1], condition_bound=1,
                      trace=lambda *state: trace.append(state))
        self.assertEqual(trace, [])

    def test_explicit_reference_needs_no_native_backend(self):
        matrix, rhs = [[1, Fraction(-1, 3)], [0, 1]], [1, -1]
        fast = solve(2, sparse(matrix), rhs, condition_bound=2)
        with patch("aim044_lift.reconstruct._library", side_effect=RuntimeError("Unavailable FLINT")):
            slow = self.assert_accuracy(solve(2, sparse(matrix), rhs, condition_bound=2,
                                             reconstruction="euclid"), matrix, rhs)
        self.assertEqual(fast.x, slow.x)
        self.assertEqual(slow.config.reconstruction, "euclid")

    def test_deterministic_selection_and_generator_inputs(self):
        matrix, rhs = [[1, Fraction(1, 3)], [0, 1]], [1, -1]
        first = solve(2, sparse(matrix), rhs, condition_bound=2)
        second = solve(2, iter(reversed(sparse(matrix))), iter(rhs), condition_bound=2)
        self.assertEqual(first, second)
        self.assertEqual(first.config.as_dict(), second.config.as_dict())

    def test_trace_preserves_output_and_callback_failure_propagates(self):
        entries, rhs = sparse([[1, 0], [0, 1]]), [1, 1]
        plain = solve(2, entries, rhs, condition_bound=1)
        states = []
        traced = solve(2, entries, rhs, condition_bound=1,
                       trace=lambda *state: states.append(state))
        self.assertEqual(plain, traced)
        self.assertEqual(len(states), traced.lifts)
        for iteration, residual, digit in states:
            self.assertTrue(all(0 <= value < 1 << traced.radix_bits for value in digit))
            self.assertTrue(all(abs(value).bit_length() <= traced.residual_bit_bound for value in residual))
        def fail(*state):
            raise LookupError("callback stopped")
        with self.assertRaisesRegex(LookupError, "callback stopped"):
            solve(2, entries, rhs, condition_bound=1, trace=fail)

    def test_long_output_and_configuration_remain_exact_json(self):
        result = solve(1, [(0, 0, 1)], [1], c=15000)
        encoded = json.loads(json.dumps(result.as_dict()))
        self.assertGreater(result.x[0].numerator.bit_length(), 14000)
        self.assertEqual(int(encoded["x"][0]["numerator_hex"], 16), result.x[0].numerator)
        bound = encoded["config"]["relative_error_bound"]
        self.assertEqual(Fraction(int(bound["numerator_hex"], 16), int(bound["denominator_hex"], 16)), Fraction(1, 10))

    def test_finite_float_inputs_keep_their_exact_binary_value(self):
        matrix, rhs = [[1, 0], [0, 1]], [0.1, 1]
        result = self.assert_accuracy(solve(2, sparse(matrix), rhs, condition_bound=1), matrix, rhs)
        self.assertGreater(result.config.observed_input_bits, 4)


class InstalledCLITests(unittest.TestCase):
    def run_cli(self, *args):
        env = dict(os.environ)
        env.pop("PYTHONPATH", None)
        return subprocess.run([sys.executable, "-m", "aim044_lift", *map(str, args)],
                              cwd=ROOT, env=env, text=True, capture_output=True, timeout=15)

    def test_default_json_and_cli_overrides_share_api(self):
        matrix, rhs = [[-1, 0], [0, Fraction(1, 3)]], [1, 1]
        for flags, options in (([], {}), (["--reconstruction", "euclid", "--relative-tolerance", "1/1000"],
                                         {"reconstruction": "euclid", "relative_tolerance": Fraction(1, 1000)})):
            with self.subTest(flags=flags):
                process = self.run_cli(ROOT/"numerical_case/example-auto.json", *flags)
                self.assertEqual(process.returncode, 0, process.stderr)
                expected = solve(2, sparse(matrix), rhs, condition_bound=3, **options)
                self.assertEqual(json.loads(process.stdout), expected.as_dict())

    def test_legacy_json_and_explicit_c_override(self):
        process = self.run_cli(ROOT/"numerical_case/example.json", "--c", "3")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(json.loads(process.stdout)["config"]["c"], 3)

    def test_invalid_cli_reports_actionable_error(self):
        process = self.run_cli(ROOT/"numerical_case/example-auto.json", "--reconstruction", "auto")
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("invalid choice", process.stderr)
        raw = ROOT/"numerical_case/raw"
        raw.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=raw) as directory:
            path = Path(directory)/"missing.json"
            data = json.loads((ROOT/"numerical_case/example-auto.json").read_text())
            data.pop("condition_bound")
            path.write_text(json.dumps(data))
            process = self.run_cli(path)
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("Provide condition_bound", process.stderr)
            for bad in (True, float("inf"), float("nan"), [3], [3, 1, 2], [True, 1], [3.0, 1]):
                data["condition_bound"] = bad
                path.write_text(json.dumps(data))
                process = self.run_cli(path)
                self.assertNotEqual(process.returncode, 0)
                self.assertIn("condition_bound", process.stderr)
                self.assertNotIn("Traceback", process.stderr)


if __name__ == "__main__":
    unittest.main()
