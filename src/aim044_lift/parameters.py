"""Exact, shared problem-contract selection for the API and JSON CLI."""

from dataclasses import asdict, dataclass
from fractions import Fraction

DEFAULT_RECONSTRUCTION = "flint"
RECONSTRUCTIONS = ("flint", "euclid")


def _fraction(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite rational value, not bool")
    try:
        return Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite rational value") from error


def _exponent(n, target):
    """Smallest positive c with n**c >= target, using exact comparisons."""
    lo, hi = 1, max(1, target.numerator.bit_length())
    while lo < hi:
        mid = (lo + hi) // 2
        if n**mid * target.denominator >= target.numerator:
            hi = mid
        else:
            lo = mid + 1
    return lo


def _ratio_dict(value):
    if value is None:
        return None
    return {"numerator_hex": hex(value.numerator),
            "denominator_hex": hex(value.denominator)}


@dataclass(frozen=True)
class SolveConfiguration:
    """Effective problem contract; matrix promises are not certified."""

    c: int
    c_source: str
    condition_bound: Fraction
    condition_bound_source: str
    requested_condition_bound: Fraction | None
    requested_relative_tolerance: Fraction | None
    relative_error_bound: Fraction
    observed_input_bits: int
    reconstruction: str
    selection_reasons: tuple[str, ...]
    matrix_promises_verified: bool = False

    def as_dict(self):
        """JSON-ready configuration, with exact hexadecimal rational pairs."""
        result = asdict(self)
        for name in ("condition_bound", "requested_condition_bound",
                     "requested_relative_tolerance", "relative_error_bound"):
            result[name] = _ratio_dict(getattr(self, name))
        result["selection_reasons"] = list(self.selection_reasons)
        return result


def select_contract(n, observed_bits, c, condition_bound, relative_tolerance,
                    reconstruction):
    """Resolve one deterministic policy; called by solve, never duplicated in CLI."""
    if c is not None and (not isinstance(c, int) or isinstance(c, bool) or c < 1):
        raise ValueError("c must be a positive integer or None")
    bound = None if condition_bound is None else _fraction(condition_bound, "condition_bound")
    tolerance = (None if relative_tolerance is None
                 else _fraction(relative_tolerance, "relative_tolerance"))
    if bound is not None and bound < 1:
        raise ValueError("condition_bound must be at least 1")
    if tolerance is not None and not 0 < tolerance <= 1:
        raise ValueError("relative_tolerance must lie in (0, 1]")
    if (reconstruction is not None
            and (not isinstance(reconstruction, str) or reconstruction not in RECONSTRUCTIONS)):
        raise ValueError("reconstruction must be None, 'flint', or 'euclid'; no fallback is used")
    backend = DEFAULT_RECONSTRUCTION if reconstruction is None else reconstruction
    reasons = ["Default fast FLINT reconstruction; unavailable means an error."
               if reconstruction is None else f"Explicit {backend} reconstruction."]
    if c is None and bound is None and n != 1:
        raise ValueError("Provide condition_bound (a justified upper bound on kappa_2(A)) "
                         "for automatic selection, or supply the legacy c promise; "
                         "conditioning is not inferred from the entries")
    if n == 1 and tolerance is not None and tolerance < Fraction(1, 10):
        raise ValueError("For n=1 this policy only guarantees relative error 1/10; "
                         "a smaller relative_tolerance cannot be promised by increasing c")
    unit_bits = (2 * n - 1).bit_length()
    input_c = max(1, (observed_bits + unit_bits - 1) // unit_bits)
    if n == 1:
        effective_bound = Fraction(1)
        bound_source = "exact nonzero-scalar condition number"
        condition_c = 1
        reasons.append("For a validated nonzero 1-by-1 scalar, kappa_2 is exactly 1.")
    elif bound is not None:
        effective_bound = bound
        bound_source = "caller-supplied, unchecked matrix promise"
        condition_c = _exponent(n, bound)
    else:
        effective_bound = Fraction(n**c)
        bound_source = "implicit legacy n**c promise, unchecked"
        condition_c = 1
    accuracy_c = (1 if tolerance is None or n == 1
                  else _exponent(n, 1 / (10 * tolerance)))
    if c is None:
        c = max(input_c, condition_c, accuracy_c)
        c_source = "automatic"
        reasons.append(f"Smallest admissible positive c: max(input={input_c}, "
                       f"condition={condition_c}, accuracy={accuracy_c}).")
    else:
        c_source = "explicit legacy contract"
        if input_c > c:
            raise ValueError("Input exceeds the promised numerator/denominator bit bound; "
                             "increase c or use automatic selection with condition_bound")
        if condition_c > c:
            raise ValueError("condition_bound exceeds n**c; increase c or omit c for automatic selection")
        if accuracy_c > c:
            raise ValueError("relative_tolerance is stricter than n**(-c)/10; "
                             "increase c or omit c with a justified condition_bound")
        reasons.append("Explicit c retained: no change to the legacy rounding or error contract.")
    error_bound = Fraction(1, 10 * n**c)
    reasons.append("Default accuracy is exactly the conditional n**(-c)/10 guarantee."
                   if tolerance is None else "Selected c meets the requested relative_tolerance exactly by rational comparison.")
    reasons.append("h=4+B+(2c+1)*ceil(log2(n)); radix bits round upward to a byte; L=6n+2.")
    return SolveConfiguration(c, c_source, effective_bound, bound_source, bound,
                              tolerance, error_bound, observed_bits, backend, tuple(reasons))
