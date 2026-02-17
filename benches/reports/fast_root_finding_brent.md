# Fast Root Refinement in Siderust (Brent + Endpoint Reuse)

> Note: This report lives under `benches/reports/` because it documents
> benchmark-driven performance work.

Many of Siderust’s “find events over time” problems reduce to the same shape:

1. scan a time window to **bracket** crossings of a scalar function, then  
2. **refine** each bracket to a precise timestamp, then  
3. assemble the resulting events into periods.

This report documents why the crate’s numerical core centers around Brent’s
method and why “endpoint reuse” is one of the highest-leverage optimizations in
these pipelines.

## Where this lives in the codebase

The numerical building blocks are under:

- `../../src/calculus/math_core/root_finding.rs` (Brent + bisection)
- `../../src/calculus/math_core/intervals.rs` (scan→refine→label→assemble)

Solar/Lunar altitude routines (and the higher-level altitude API) delegate into
these generic pieces rather than implementing their own solvers.

## Why Brent is a good fit here

In this codebase, the expensive part is usually not arithmetic inside the root
solver. It is the function evaluation itself: “Sun altitude at `t`”, “Moon
altitude at `t`”, “hour angle at `t`”, and so on. Those evaluations typically
trigger ephemeris work and trigonometric transforms.

Brent’s method is a pragmatic choice because it is:

it is derivative-free (no finite-difference probes), it maintains a bracket
like bisection (so it is hard to “explode” numerically), and it tends to
converge with relatively few *new* function evaluations.

That makes it a strong default when the caller already has a sign-change bracket
from scanning.

## The key optimization: reuse scan endpoint evaluations

When you scan a window, you must evaluate the function at both endpoints of
each step to detect sign changes:

```text
g(t) = f(t) - threshold
```

If `g(t0)` and `g(t1)` have opposite signs, there is a root in `[t0, t1]`. At
that moment you already have the values `g(t0)` and `g(t1)` in hand.

A “naive” refiner will recompute them inside the solver, paying for two
redundant expensive evaluations per crossing. Over long horizons, those add up.

To avoid this, the crate provides:

- `brent_with_values(period, g_lo, g_hi, g)` in `root_finding.rs`

and the interval engine uses it directly when a sign-change is detected. The
hot path in `intervals::find_crossings` is intentionally structured as:

- evaluate `g(t)` and `g(next_t)` during the scan
- if signs differ, call `brent_with_values(Period::new(t, next_t), prev, next_v, g)`

So refinement starts without re-paying the endpoint cost.

## Crossing direction classification (and avoiding extra probes)

After finding a crossing time, many algorithms want to classify it as “entering”
or “exiting” an above-threshold region. The generic interval engine can do this
by probing `f(t−ε)` and `f(t+ε)` (see `PROBE_DT` in `intervals.rs`). That is
simple and robust, but it costs two additional expensive evaluations per
crossing.

In cases where the scan already knows the sign direction across the bracket, we
can label the crossing without extra probes. The Moon altitude pipeline uses
this approach via `find_and_label_crossings` in
`../../src/calculus/lunar/moon_cache.rs`.

## Tolerance: a deliberate speed/precision trade-off

Root refinement terminates when the bracket is sufficiently small in the time
domain. The default tolerance in `math_core::root_finding` is very strict (on
the order of microseconds when the domain unit is days), which is appropriate
for correctness tests but may be tighter than planning-grade outputs need.

If a workload is dominated by refinement iterations, using `brent_tol` to relax
the time tolerance is a straightforward lever: it reduces iterations without
changing the overall scan/bracketing logic.

## Why this matters for benchmarks

This is most visible in long-horizon benchmarks (30/365 days) where the scan
produces many brackets and each bracket needs refinement.

For example:

- Sun night/day windows in `../solar_altitude.rs`
- Moon above/below-horizon windows in `../moon_altitude.rs`
