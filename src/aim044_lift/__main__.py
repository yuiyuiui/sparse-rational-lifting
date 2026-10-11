"""JSON input and hexadecimal rational output for reproducible sparse solves."""

import argparse
import json
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path

from .solver import solve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--reconstruction", choices=["flint", "euclid"], default="flint")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    result = solve(data["n"], [(i, j, Fraction(a, b)) for i, j, a, b in data["entries"]],
                   [Fraction(a, b) for a, b in data["b"]], data["c"],
                   reconstruction=args.reconstruction)
    info = asdict(result)
    info["x"] = [{"numerator_hex": hex(x.numerator), "denominator_hex": hex(x.denominator)}
                 for x in result.x]
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
