"""Installed JSON CLI sharing all parameter-selection rules with solve()."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

from .parameters import RECONSTRUCTIONS
from .solver import solve


def _rational(value, name):
    if isinstance(value, list):
        if len(value) != 2 or any(not isinstance(part, int) or isinstance(part, bool) for part in value):
            raise ValueError(f"{name} JSON pair must contain two integers")
        return Fraction(*value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--c", type=int, help="Explicit legacy problem exponent; overrides JSON c")
    parser.add_argument("--condition-bound", help="Justified rational upper bound on kappa_2(A)")
    parser.add_argument("--relative-tolerance", help="Requested relative solution error, e.g. 1/1000")
    parser.add_argument("--reconstruction", choices=RECONSTRUCTIONS,
                        help="Recovery backend (API default: flint, never a silent fallback)")
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text())
        options = {"c": data.get("c"), "condition_bound": data.get("condition_bound"),
                   "relative_tolerance": data.get("relative_tolerance"),
                   "reconstruction": data.get("reconstruction")}
        for key in options:
            value = getattr(args, key)
            if value is not None:
                options[key] = value
        for key in ("condition_bound", "relative_tolerance"):
            if options[key] is not None:
                options[key] = _rational(options[key], key)
        result = solve(data["n"], [(i, j, Fraction(a, b)) for i, j, a, b in data["entries"]],
                       [Fraction(a, b) for a, b in data["b"]], **options)
    except (KeyError, TypeError, ValueError, ZeroDivisionError, RuntimeError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(result.as_dict(), indent=2))


if __name__ == "__main__":
    main()
