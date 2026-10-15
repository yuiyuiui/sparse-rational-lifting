"""Fixed matched default/explicit-backend comparisons with exact validation."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import resource
import signal
import statistics
import subprocess
import time
import tracemalloc
from datetime import datetime, timezone
from pathlib import Path

import flint
from aim044_lift import solve
from benchmark import ROOT, gaussian_reference, matrix_family, sparse

POLICIES = ("automatic_default", "explicit_flint", "explicit_euclid")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    signal.alarm(60)
    started = datetime.now(timezone.utc).isoformat()
    solve(1, [(0, 0, 1)], [1])  # native-loader warm-up, outside comparisons
    cases = []
    for family in ("conditioning_boundary", "small_rational_scale",
                   "indefinite_tridiagonal", "nonnormal_bidiagonal"):
        for n in (4, 8, 16, 32):
            matrix, rhs = matrix_family(family, n)
            entries = sparse(matrix)
            reference = gaussian_reference(matrix, rhs)
            norm_squared = sum(x*x for x in reference)
            # All four families satisfy this common conservative problem promise.
            options = {"automatic_default": {"condition_bound": n*n},
                       "explicit_flint": {"c": 2, "reconstruction": "flint"},
                       "explicit_euclid": {"c": 2, "reconstruction": "euclid"}}
            outputs, seconds, cpu_seconds = {}, {key: [] for key in POLICIES}, {key: [] for key in POLICIES}
            for repeat in range(3):
                order = POLICIES[repeat:] + POLICIES[:repeat]
                for policy in order:
                    wall, cpu = time.perf_counter(), time.process_time()
                    result = solve(n, entries, rhs, **options[policy])
                    cpu_seconds[policy].append(time.process_time()-cpu)
                    seconds[policy].append(time.perf_counter()-wall)
                    if policy in outputs and outputs[policy].x != result.x:
                        raise AssertionError("Repeated output changed")
                    outputs[policy] = result
            if not all(outputs[key].x == outputs[POLICIES[0]].x for key in POLICIES):
                raise AssertionError("Matched policies returned different rational outputs")
            rows = []
            for policy in POLICIES:
                result = outputs[policy]
                if result.config.c != 2:
                    raise AssertionError("Comparison accuracy/rounding contract changed")
                error_squared = sum((value.as_fraction()-target)**2
                                    for value, target in zip(result.x, reference))
                exact_pass = error_squared * 100 * n**4 <= norm_squared
                width_pass = result.max_residual_bits <= result.residual_bit_bound
                if not exact_pass or not width_pass:
                    raise AssertionError((family, n, policy))
                # Separate instrumented call: excluded from production timings.
                tracemalloc.start()
                measured = solve(n, entries, rhs, **options[policy])
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                if measured.x != result.x:
                    raise AssertionError("Instrumented output changed")
                rows.append({"policy": policy, "effective_configuration": result.config.as_dict(),
                             "rounding_bits": result.rounding_bits, "radix_bits": result.radix_bits,
                             "lifts": result.lifts, "sparse_products": result.sparse_products,
                             "digit_storage_bytes": n*result.lifts*(result.radix_bits//8),
                             "python_tracemalloc_peak_bytes": peak,
                             "seconds": seconds[policy], "median_seconds": statistics.median(seconds[policy]),
                             "cpu_seconds": cpu_seconds[policy],
                             "relative_l2_error": math.sqrt(float(error_squared/norm_squared)),
                             "exact_error_check": exact_pass, "width_check": width_pass,
                             "rational_output_identical_across_policies": True})
            cases.append({"family": family, "n": n, "nnz": len(entries), "policies": rows})
    signal.alarm(0)
    relevant = sorted((ROOT/"src").rglob("*.py")) + [ROOT/"test/test_solver.py", ROOT/"test/test_api.py",
                                                           ROOT/"numerical_case/benchmark.py", Path(__file__)]
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() != "Darwin":
        peak_rss *= 1024
    report = {
        "schema_version": 1, "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": platform.python_version(), "system": platform.system(),
                        "machine": platform.machine(), "package": importlib.metadata.version("sparse-rational-lifting"),
                        "python_flint": flint.__version__, "flint": flint.__FLINT_VERSION__,
                        "logical_cpus": os.cpu_count(), "native_threads": 1,
                        "openblas_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
                        "omp_threads": os.environ.get("OMP_NUM_THREADS")},
        "source_base_commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in relevant},
        "contract": {"input_promise": "kappa_2(A)<=n^2, c=2 input precision; identical inputs for all modes",
                     "accuracy": "Relative Euclidean error <=n^-2/10; exact squared-fraction acceptance",
                     "timing": "Three rotated warm repeats; full solve including normalization, selection, backend check, rounding, lifting and recovery. Reference and checks excluded.",
                     "memory": "Separate tracemalloc call per policy excludes native allocations; RSS is entire-process high-water including references/instrumentation, not per-mode memory.",
                     "seeds": "None: deterministic fixed cases", "timeout_seconds": 60,
                     "evidence": "NUMERICAL EVIDENCE: finite checks, no optimality/scaling or asymptotic speedup claim"},
        "whole_process_peak_rss_bytes": peak_rss, "cases": cases,
    }
    (args.output/"policy_results.json").write_text(json.dumps(report, indent=2)+"\n")
    lines = ["# Matched API/default-policy benchmark", "", report["contract"]["evidence"], "",
             "All 16 input cases × 3 policies passed exact error and width checks, and returned identical rational outputs.", "",
             "| Family | n | Default median ms | Explicit FLINT ms | Explicit Euclid ms |", "| --- | ---: | ---: | ---: | ---: |"]
    for case in cases:
        values = [row["median_seconds"]*1000 for row in case["policies"]]
        lines.append(f'| {case["family"]} | {case["n"]} | '+" | ".join(f"{value:.4f}" for value in values)+" |")
    lines += ["", report["contract"]["timing"], "", report["contract"]["memory"], "",
              "Environment: "+json.dumps(report["environment"])+".",
              f"Whole-process peak RSS: {peak_rss} bytes.",
              "Timing differences on small inputs include noise and do not establish a generally fastest policy."]
    (args.output/"POLICY_SUMMARY.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"input_cases": len(cases), "policy_cases": 3*len(cases),
                      "all_exact_checks_passed": True, "identical_rational_outputs": True,
                      "whole_process_peak_rss_bytes": peak_rss}))


if __name__ == "__main__":
    main()
