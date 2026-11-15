# Root refinement and endpoint reuse

An altitude crossing search solves a scalar residual such as
`sin(altitude(t)) − sin(threshold)`. With a continuous residual, opposite
endpoint signs bracket at least one root. An even number of roots or a
tangency can leave both signs equal. Search completeness therefore depends on
how brackets are discovered and validated.

The [typed Brent solver](../../src/event/search/root_finding.rs) combines
bisection, secant and inverse quadratic interpolation. Its `brent` wrapper
evaluates both endpoints. `brent_with_values` and `brent_tol` accept values
already sampled by the caller. When both endpoints are available, that avoids
two function calls per bracket relative to a wrapper that recomputes them.
This follows from the call structure. Elapsed savings depend on the residual
cost, number of brackets and fallback behavior.

The [generic interval helpers](../../src/event/search/intervals.rs) carry values
from one scan step to the next and pass them into refinement. The directed
helper derives rising or setting from ordinary sign changes. The general
`label_crossings` helper samples either side of a candidate when it needs
direction. Handling zeros at interval endpoints can require additional probes.

The [altitude scan fallback](../../src/event/search/scan_fallback.rs) uses its
own scalar `brent_f64` with supplied endpoint values. The
[Chebyshev search](../../src/event/search/crossings.rs) validates polynomial
candidates against the workload's residual and can fall back to this scan.
Solar prediction can also use a Newton/secant polish before using Brent.
See the [solar engine report](solar_altitude_perf_report.md).

The typed solver's default of `1e-9` days is 86.4 µs when the domain is time
in days. The altitude search configuration also defaults to `1e-9` days.
Residual criteria and coordinate or model errors affect the final event time.
A small root tolerance alone supplies no physical timing bound, especially
for a shallow crossing.

The solver, scan fallback and interval modules contain tests for convergence,
endpoint zeros, direction and period assembly. The altitude benchmark suites
time complete workloads. To attribute timing differences to refinement, collect
residual evaluation counts for matching brackets and verify crossing counts,
directions and times before comparing elapsed cost. No standalone root solver
timing is recorded here.
