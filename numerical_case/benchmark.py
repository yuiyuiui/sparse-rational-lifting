"""Fixed, bounded accuracy and timing reproduction on admissible sparse families."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import signal
import statistics
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "test"))
from test_solver import gaussian_reference, sparse
from aim044_lift import solve


def matrix_family(name, n):
    matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        if name == "conditioning_boundary":
            matrix[i][i] = Fraction(1, n*n) if i == n-1 else Fraction(1)
        elif name == "small_rational_scale":
            matrix[i][i] = Fraction(1, 2*n*n-1)
        elif name == "indefinite_tridiagonal":
            matrix[i][i] = Fraction(2 if i % 2 else -2)
            if i:
                matrix[i][i-1] = Fraction(1, 4)
            if i+1 < n:
                matrix[i][i+1] = Fraction(1, 4)
        elif name == "nonnormal_bidiagonal":
            matrix[i][i] = Fraction(1)
            if i+1 < n:
                matrix[i][i+1] = Fraction(1, 4)
        else:
            raise ValueError(name)
    rhs = [Fraction((-1)**i, 2*i+1) for i in range(n)]
    return matrix, rhs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    signal.alarm(60)
    started = datetime.now(timezone.utc).isoformat()
    cases = []
    solve(1, [(0, 0, Fraction(1))], [1], 2)  # warm native loader
    for family in ("conditioning_boundary", "small_rational_scale",
                   "indefinite_tridiagonal", "nonnormal_bidiagonal"):
        for n in (4, 8, 16, 32):
            matrix, rhs = matrix_family(family, n)
            entries = sparse(matrix)
            reference = gaussian_reference(matrix, rhs)
            times = []
            result = None
            for _ in range(3):
                start = time.perf_counter()
                result = solve(n, entries, rhs, 2)
                times.append(time.perf_counter()-start)
            output = [value.as_fraction() for value in result.x]
            error_squared = sum((x-y)**2 for x, y in zip(output, reference))
            norm_squared = sum(x*x for x in reference)
            exact_pass = error_squared * 100 * n**4 <= norm_squared
            width_pass = result.max_residual_bits <= result.residual_bit_bound
            if not exact_pass or not width_pass:
                raise AssertionError((family, n, exact_pass, width_pass))
            cases.append({
                "family": family, "n": n, "c": 2, "nnz": result.nnz,
                "lifts": result.lifts, "sparse_products": result.sparse_products,
                "rounding_bits": result.rounding_bits, "radix_bits": result.radix_bits,
                "max_residual_bits": result.max_residual_bits,
                "residual_bit_bound": result.residual_bit_bound,
                "relative_l2_error": math.sqrt(float(error_squared/norm_squared)),
                "exact_stronger_error_check": exact_pass, "width_check": width_pass,
                "seconds": times, "median_seconds": statistics.median(times),
                "reconstruction": result.reconstruction,
                "effective_configuration": result.config.as_dict(),
            })
    signal.alarm(0)
    relevant = list((ROOT/"src").rglob("*.py")) + [ROOT/"test/test_solver.py", Path(__file__)]
    report = {
        "schema_version": 1, "started_at": started,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": platform.python_version(), "platform": platform.system(),
                        "machine": platform.machine(),
                        "python_flint": importlib.metadata.version("python-flint"),
                        "library_threads": 1},
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in relevant},
        "timing_contract": "Three warm repeats; solver call only, excludes exact reference and error checks.",
        "evidence_status": "NUMERICAL EVIDENCE: finite implementation validation; not an asymptotic proof.",
        "cases": cases,
    }
    (args.output/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    rows = ["# Reproduced cases", "", report["evidence_status"], "",
            "| Family | n | m | Lifts | Sparse products | Relative error | Median seconds |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for case in cases:
        rows.append(f'| {case["family"]} | {case["n"]} | {case["nnz"]} | {case["lifts"]} | '
                    f'{case["sparse_products"]} | {case["relative_l2_error"]:.3e} | '
                    f'{case["median_seconds"]:.6f} |')
    rows += ["", "All 16 cases passed the exact stronger forward-error inequality and the residual-width check.",
             "Times describe this fixed Python/FLINT implementation on the recorded machine."]
    (args.output/"SUMMARY.md").write_text("\n".join(rows)+"\n")
    print(json.dumps({"cases":len(cases),"all_exact_checks_passed":True,
                      "max_relative_error":max(c["relative_l2_error"] for c in cases),
                      "max_median_seconds":max(c["median_seconds"] for c in cases)}))


if __name__ == "__main__":
    main()
