"""Collect one fixed, independently checked lifting trace; never benchmark it."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import signal
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from benchmark import ROOT, gaussian_reference, matrix_family, sparse
from aim044_lift import solve


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    signal.alarm(60)
    published_path = ROOT / "numerical_case/published/results.json"
    published = json.loads(published_path.read_text())
    if not all(case["exact_stronger_error_check"] and case["width_check"]
               for case in published["cases"]):
        raise AssertionError("Published benchmark contains a failed validation flag")
    for name, expected_hash in published["source_sha256"].items():
        if digest(ROOT / name) != expected_hash:
            raise AssertionError(f"Published benchmark source has changed: {name}")

    family, n, c = "conditioning_boundary", 32, 2
    matrix, rhs = matrix_family(family, n)
    reference = gaussian_reference(matrix, rhs)
    h = 4 + c * (2 * n - 1).bit_length() + (2 * c + 1) * (n - 1).bit_length()

    def rounded(value):
        shifted = value * (1 << h) + Fraction(1, 2)
        return shifted.numerator // shifted.denominator

    converted = [(row, col, rounded(value)) for row, col, value in sparse(matrix)]
    converted = [(row, col, value) for row, col, value in converted if value]
    expected_residual = [rounded(value) for value in rhs]
    height = max([1] + [abs(v) for _, _, v in converted] + list(map(abs, expected_residual)))
    radix_bits = 8 * (((2 * n * height - 1).bit_length() + 7) // 8)
    radix = 1 << radix_bits
    rows = []

    def collect(iteration, residual, digit):
        nonlocal expected_residual
        if iteration != len(rows) or tuple(expected_residual) != residual:
            raise AssertionError("Trace disagrees with the independently computed recurrence")
        if tuple(value % radix for value in residual) != digit:
            raise AssertionError("Trace digit is not the least nonnegative radix residue")
        rows.append({
            "completed_lifts": iteration,
            "integer_residual_bits": max(abs(value).bit_length() for value in residual),
            "stored_digit_capacity_bits_per_coordinate": iteration * radix_bits,
        })
        product = [0] * n
        for row, col, value in converted:
            product[row] += value * digit[col]
        expected_residual = [(value - digit[i]) // radix - product[i]
                             for i, value in enumerate(residual)]

    result = solve(n, sparse(matrix), rhs, c, trace=collect)
    rows.append({
        "completed_lifts": result.lifts,
        "integer_residual_bits": max(abs(value).bit_length() for value in expected_residual),
        "stored_digit_capacity_bits_per_coordinate": result.lifts * radix_bits,
    })
    output = [value.as_fraction() for value in result.x]
    error_squared = sum((x - y) ** 2 for x, y in zip(output, reference))
    norm_squared = sum(x * x for x in reference)
    relative_error = math.sqrt(float(error_squared / norm_squared))
    old_case = next(case for case in published["cases"]
                    if case["family"] == family and case["n"] == n)
    checks = {
        "all_callback_states_and_digits_match_independent_recurrence": len(rows) == result.lifts + 1,
        "exact_stronger_error_inequality": error_squared * 100 * n ** (2 * c) <= norm_squared,
        "all_integer_residual_widths_within_bound": all(
            row["integer_residual_bits"] <= result.residual_bit_bound for row in rows),
        "summary_matches_solver": max(row["integer_residual_bits"] for row in rows) == result.max_residual_bits
            and radix_bits == result.radix_bits and h == result.rounding_bits,
        "matches_published_case": all(getattr(result, key) == old_case[key] for key in (
            "lifts", "sparse_products", "radix_bits", "rounding_bits", "max_residual_bits", "residual_bit_bound"))
            and relative_error == old_case["relative_l2_error"],
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    report = {
        "schema_version": 1,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "case": {"family": family, "n": n, "c": c, "nnz": result.nnz},
        "radix_bits": radix_bits, "lifts": result.lifts,
        "residual_bit_bound": result.residual_bit_bound,
        "relative_l2_error": relative_error,
        "environment": {"python": platform.python_version(), "system": platform.system(),
                        "machine": platform.machine(),
                        "python_flint": importlib.metadata.version("python-flint")},
        "source_sha256": {**published["source_sha256"],
                          "numerical_case/trace_case.py": digest(Path(__file__))},
        "benchmark_sha256": digest(published_path),
        "validation": checks,
        "interpretation": "Encoded integer lifting residual, not physical forward error. Stored digit capacity is k*s bits per coordinate, not the final rational numerator size. The last residual is derived from the checked last callback by the exact recurrence. Tracing and Gaussian validation are outside solver timings.",
        "trace": rows,
    }
    (args.output / "lifting_trace.json").write_text(json.dumps(report, indent=2) + "\n")
    signal.alarm(0)
    print(json.dumps({"case": report["case"], "trace_points": len(rows), "validation": checks}))


if __name__ == "__main__":
    main()
