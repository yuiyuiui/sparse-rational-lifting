# 0.2.1 — measured core and delivery improvements

- Reuse native modulus/height conversions within each solve.
- Release digit buffers during recovery and remove repeated residual scans.
- Preserve the caller's FLINT thread setting.
- Add independent batch/cleanup/state tests and paired time/memory profiles.
- Add installed-wheel CI, citation and contribution metadata.
- Public signatures, promises and exact rational outputs are preserved.

# Changes

## 0.2.0 — 15 October 2026

The reusable solve API, exact lifting arithmetic, rational output, fixed
rounding/height rules, real FLINT/Euclid modes and independent references were
already present. This maintenance adds a safe parameter-selection and
installed-use contract without changing the underlying solver or a scientific
problem's status.

- Existing `solve(n, entries, rhs, c, reconstruction=...)` calls retain their
  numerical output and mathematical promises. c may now be omitted in favor of
  an explicitly justified `condition_bound`. A nonzero 1-by-1 scalar supplies
  its exact condition number one automatically.
- `relative_tolerance` is an optional problem target. Exact selection covers
  input bits, the supplied condition bound and this target. Missing general
  conditioning information or incompatible explicit c fails clearly; n=1's
  supported guarantee remains 1/10. No c=2 default or costly condition estimate
  is invented.
- `reconstruction=None` explicitly denotes the previous fast FLINT default.
  Both explicit backend names remain supported. Native availability is checked
  before lifting; no reference fallback or precision downgrade occurs.
- `SolveResult.config`, `SolveConfiguration` and `as_dict()` expose effective
  parameters, exact rational bounds and selection reasons. JSON adds a config
  object; strict schema consumers must accept this additive field. Existing
  result data fields remain. Manual positional result construction can omit
  the new optional config. Dataclass equality/repr includes it; compare x when
  checking numerical output. Use `as_dict()` for JSON serialization.
- Callback types are checked before work. CLI overrides and optional JSON
  fields use the same resolver as the callable API. CLI invalid/non-finite
  parameter errors are actionable.
- The package is tested from a wheel installation: 30 tests. A matched
  16-input × 3-policy benchmark records exact output agreement, conditional
  error checks, complete production timings, effective configuration, source
  identity and qualified memory metrics. Existing cases, trace and four
  required figures are regenerated for the current source.
- README/API documentation now give every parameter's meaning, default,
  supported choices, recommendations and cost/guarantee limits, and explicitly
  disclose creation and maintenance by AI agents under the author's direction.

No independent human review, formal certification or general default-optimality
claim is added. Prior immutable commits and historical inventories are preserved.
Dependencies, native reconstruction algorithm and precision contract are unchanged.
