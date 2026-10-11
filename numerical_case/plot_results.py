"""Render four vector figures from selected validated data, without running solvers."""

import argparse
import hashlib
import json
import math
import platform
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = {
    "conditioning_boundary": ("Diagonal: conditioning boundary", "#2066a8", "o"),
    "small_rational_scale": ("Diagonal: small rational scale", "#a14d12", "s"),
    "indefinite_tridiagonal": ("Indefinite tridiagonal", "#277d60", "^"),
    "nonnormal_bidiagonal": ("Nonnormal bidiagonal", "#8751a0", "D"),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "numerical_case/published/results.json")
    parser.add_argument("--trace", type=Path, default=ROOT / "numerical_case/published/lifting_trace.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview", action="store_true", help="Also render PNGs for local visual inspection")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    data, trace = json.loads(args.results.read_text()), json.loads(args.trace.read_text())
    if trace["benchmark_sha256"] != digest(args.results) or not all(trace["validation"].values()):
        raise AssertionError("Trace validation or benchmark identity failed")
    if (trace["case"] != {"family": "conditioning_boundary", "n": 32, "c": 2, "nnz": 32}
            or trace["lifts"] != 194 or trace["radix_bits"] != 48 or trace["residual_bit_bound"] != 90):
        raise AssertionError("Expected the fixed conditioning-boundary trace")
    for report in (data, trace):
        for name, expected in report["source_sha256"].items():
            if digest(ROOT / name) != expected:
                raise AssertionError(f"Source identity failed: {name}")
    cases = data["cases"]
    if len(cases) != 16 or {(row["family"], row["n"]) for row in cases} != {
            (name, n) for name in FAMILIES for n in (4, 8, 16, 32)}:
        raise AssertionError("Expected the fixed four-by-four benchmark")
    for row in cases:
        if not row["exact_stronger_error_check"] or not row["width_check"] or row["c"] != 2:
            raise AssertionError("A published validation flag failed")
        if not 0 < row["relative_l2_error"] <= row["n"] ** -2 / 10:
            raise AssertionError("Invalid relative error for log axes")
        if (row["lifts"] != 6 * row["n"] + 2
                or row["sparse_products"] != row["lifts"] * row["nnz"]
                or row["max_residual_bits"] > row["residual_bit_bound"]):
            raise AssertionError("Unexpected fixed-case work or width counts")
        if len(row["seconds"]) != 3 or not all(math.isfinite(t) and t > 0 for t in row["seconds"]):
            raise AssertionError("Invalid timing sample")
        if statistics.median(row["seconds"]) != row["median_seconds"]:
            raise AssertionError("Stored timing median disagrees with repeats")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 15, "axes.titleweight": "bold",
                         "svg.fonttype": "none", "svg.hashsalt": "sparse-rational-lifting-v1",
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    files = []

    def save(fig, name):
        svg_path = args.output / f"{name}.svg"
        fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
        # Keep generated path lines clean in Git; newlines retain token separators.
        svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
        files.append(f"{name}.svg")
        if args.preview:
            fig.savefig(args.output / f"{name}.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

    def configure(ax):
        ax.set_xscale("log", base=2)
        ax.set_xticks([4, 8, 16, 32])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xlabel("Matrix dimension n")
        ax.grid(axis="y", alpha=0.18)
        ax.set_axisbelow(True)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    for name, (label, color, marker) in FAMILIES.items():
        subset = sorted((row for row in cases if row["family"] == name), key=lambda row: row["n"])
        ax.plot([row["n"] for row in subset], [row["relative_l2_error"] for row in subset],
                color=color, marker=marker, label=label, linewidth=2, markersize=6)
    ns = [4, 8, 16, 32]
    ax.plot(ns, [n ** -2 / 10 for n in ns], "--", color="#4c4c4c", linewidth=1.8,
            label=r"Required bound: $n^{-2}/10$")
    configure(ax)
    ax.set_yscale("log")
    ax.set_ylabel(r"Relative solution error $\|x_{out}-x_*\|_2/\|x_*\|_2$")
    ax.set_title("Accuracy on four sparse rational families", loc="left", pad=16)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    fig.subplots_adjust(bottom=0.17, top=0.88)
    fig.text(0.125, 0.025, "c = 2; 16 cases checked against exact rational Gaussian elimination.\n"
             "Lines connect sampled dimensions; finite cases do not establish an asymptotic rate.", fontsize=9, color="#555555")
    save(fig, "accuracy")

    points = trace["trace"]
    if ([row["completed_lifts"] for row in points] != list(range(trace["lifts"] + 1))
            or any(row["stored_digit_capacity_bits_per_coordinate"] != row["completed_lifts"] * trace["radix_bits"]
                   or row["integer_residual_bits"] > trace["residual_bit_bound"] for row in points)):
        raise AssertionError("Invalid trace shape or bound")
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.2), sharex=True)
    steps = [row["completed_lifts"] for row in points]
    axes[0].plot(steps, [row["integer_residual_bits"] for row in points], color="#2066a8", linewidth=1.8,
                 label="Maximum integer residual width")
    axes[0].axhline(trace["residual_bit_bound"], linestyle="--", color="#4c4c4c", label="Validated width bound: 90 bits")
    axes[0].set_ylim(0, 105)
    axes[0].set_ylabel("Encoded residual width (bits)")
    axes[0].legend(frameon=False, fontsize=9, loc="lower right")
    axes[0].set_title("Short residuals while digit blocks grow", loc="left", pad=16)
    axes[1].plot(steps, [row["stored_digit_capacity_bits_per_coordinate"] / 1024 for row in points],
                 color="#277d60", linewidth=2)
    axes[1].set_ylabel("Digit capacity / coordinate (Kibit)")
    axes[1].set_xlabel("Completed lifts k")
    axes[1].set_xlim(0, trace["lifts"])
    axes[1].set_xticks([0, 50, 100, 150, trace["lifts"]])
    axes[1].annotate("194 × 48 = 9,312 bits", xy=(194, 9312 / 1024), xytext=(76, 7.9),
                     arrowprops={"arrowstyle": "->", "color": "#555555"}, fontsize=10)
    for ax in axes:
        ax.grid(axis="y", alpha=0.18)
        ax.set_axisbelow(True)
    fig.subplots_adjust(bottom=0.19, hspace=0.25, top=0.89)
    fig.text(0.125, 0.03, "Conditioning-boundary diagonal, n = 32, c = 2; radix = 2⁴⁸.\n"
             "Encoded residual ≠ forward error; stored digit capacity ≠ final numerator size.\n"
             "Final residual point is derived from the checked last callback by the exact recurrence.", fontsize=9, color="#555555")
    save(fig, "lifting_width")

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    for name, label in (("conditioning_boundary", "Both diagonal families: m = n"),
                        ("nonnormal_bidiagonal", "Bidiagonal: m = 2n − 1"),
                        ("indefinite_tridiagonal", "Tridiagonal: m = 3n − 2")):
        _, color, marker = FAMILIES[name]
        subset = sorted((row for row in cases if row["family"] == name), key=lambda row: row["n"])
        ax.plot(ns, [row["sparse_products"] for row in subset], color=color, marker=marker,
                label=label, linewidth=2, markersize=6)
    ax.plot(ns, [(6 * n + 2) * n ** 2 for n in ns], "--", color="#4c4c4c", linewidth=1.8,
            label="Dense-entry loop reference: (6n + 2)n²")
    configure(ax)
    ax.set_yscale("log")
    ax.set_ylabel("Scalar products in the lifting loop")
    ax.set_title("Lifting work follows the stored nonzeros", loc="left", pad=16)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.subplots_adjust(bottom=0.17, top=0.88)
    fig.text(0.125, 0.025, "L = 6n + 2 lifts; exactly Lm scalar products in these cases.\n"
             "Dashed curve is an operation-count reference, not measured dense-solver timing or total bit cost.", fontsize=9, color="#555555")
    save(fig, "sparse_work")

    fig, ax = plt.subplots(figsize=(8.4, 5.5))
    for name, (label, color, marker) in FAMILIES.items():
        subset = sorted((row for row in cases if row["family"] == name), key=lambda row: row["n"])
        medians = [row["median_seconds"] * 1000 for row in subset]
        errors = [[(row["median_seconds"] - min(row["seconds"])) * 1000 for row in subset],
                  [(max(row["seconds"]) - row["median_seconds"]) * 1000 for row in subset]]
        ax.errorbar(ns, medians, yerr=errors, color=color, marker=marker, label=label,
                    linewidth=1.8, markersize=6, capsize=3)
    configure(ax)
    ax.set_yscale("log")
    ax.set_ylabel("Warm solver-call time (milliseconds)")
    ax.set_title("Recorded timings on small fixed cases", loc="left", pad=16)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.subplots_adjust(bottom=0.2, top=0.88)
    env = data["environment"]
    fig.text(0.125, 0.025, f"{env['platform']} {env['machine']}; Python {env['python']}; python-flint {env['python_flint']}; 1 library thread.\n"
             "Markers: median of 3 warm repeats. Whiskers: observed min–max, not confidence intervals.\n"
             "Excludes exact reference and tracing; finite timings do not establish scaling or speedup.", fontsize=9, color="#555555")
    save(fig, "timing")

    manifest = {
        "schema_version": 1, "evidence_status": "NUMERICAL EVIDENCE: fixed finite cases only",
        "data_sha256": {"published/results.json": digest(args.results),
                        "published/lifting_trace.json": digest(args.trace)},
        "plot_script_sha256": digest(Path(__file__)),
        "environment": {"python": platform.python_version(), "matplotlib": matplotlib.__version__},
        "validation": {"benchmark_flags_and_source_hashes": True, "trace_flags_source_hashes_and_recurrence_summary": True,
                       "work_count_identities": True, "timing_samples_and_medians": True},
        "figures_sha256": {name: digest(args.output / name) for name in files},
    }
    (args.output / "FIGURE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"figures": files, "validation": manifest["validation"]}))


if __name__ == "__main__":
    main()
