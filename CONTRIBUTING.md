# Contributing

Start with the API contract in `docs/API.md`. Include a small reproducer,
Python/python-flint versions, exact rational inputs, conditioning promise and
expected/actual behavior in a numerical issue. Do not post private inputs.

Install `.[dev]`, run `python -m unittest discover -s test -v`, and build with
`python -m build`. CLI tests deliberately require an installed package.
Run tests from outside the source directory against the built wheel as well.
Use `python -m ruff check src test bench` for static checks.

Performance contributions should include a profile, independent exact output
checks and matched time/memory evidence using `bench/benchmark_core.py`.
Preserve the input promises, default accuracy and backend failure behavior.
Large inputs and raw benchmark campaigns are optional, explicit work.

Public API changes follow semantic versioning; behavior changes need migration
notes. Bug reports and focused pull requests are welcome. Support is best effort.
Credit original algorithms and reused code. AI assistance is recorded in the
README; it is not independent mathematical certification.
