# 0.2.1 readiness evidence

| Gate | Status | Evidence |
|---|---|---|
| Value and scope | READY | Existing short-rational sparse solver; original Kelner attribution and caller conditioning promise preserved |
| Reliability | READY | 33 tests, exact Gaussian/back-substitution oracles, batch cleanup and process-thread-state controls |
| Usability | READY | Fresh installed wheel from outside the checkout; all 33 tests including CLI passed |
| Performance | READY | `bench/results/core_v0.2.1_final.json`; paired source hashes, raw samples and profiles; identical exact outputs on all three cases |
| Distribution | READY | Wheel and sdist built; wheel metadata checked; wheel rebuilt from sdist |
| Community and credit | READY | MIT license, original-method attribution, AI assistance disclosure, CITATION.cff and contribution guide |
| Hosted checks | NEEDS_WORK | GitHub PR/CI results are reported separately; workflow configuration alone is not passing CI |
| DOI / paper | NOT_APPLICABLE | No DOI registration or new paper requested |

Local validation used CPython 3.13.3 and python-flint 0.8.0 on macOS arm64.
Linux validation belongs to the configured CI job; other platforms are not
claimed as locally tested. Source tests and the source-tree CLI alone were
not treated as installed-artifact validation.

The first source-only test invocation found four CLI failures because the
old development environment did not have the package installed. After explicit
installation, and then in a separate wheel consumer, all tests passed. No
accuracy tolerance was relaxed.
