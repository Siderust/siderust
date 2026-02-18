# Moon Altitude Bench — Performance Notes (`find_moon_above_horizon_365day`)

> Note: This report lives under `benches/reports/` because it documents
> benchmark-driven performance work.

This report describes what the Moon altitude benchmarks in `benches/moon_altitude.rs`
are actually measuring, why Moon window-finding is inherently heavier than the
Sun case, and where optimizations tend to pay off.

## What is being measured

The benchmark suite covers both “single altitude evaluation” and “window
finding over time horizons”. In particular, the year-long cases exercise the
complete scan→refine→assemble pipeline:

- `Moon.above_threshold(site, period, 0°)` (above-horizon windows)
- `Moon.below_threshold(site, period, -0.5°)` (below-horizon windows)
- `Moon.altitude_periods(query)` (within an altitude band)

The suite also includes an explicit single-point cost benchmark
(`moon_altitude_single/compute_altitude`) to anchor the cost model.

## Current call path (today)

At a high level, long-horizon Moon window finding goes through:

1. the altitude API (`AltitudePeriodsProvider` on `Moon`),  
2. Moon-specific closures in `../../src/calculus/lunar/altitude_periods.rs`,  
3. generic scan/refine logic in `../../src/calculus/math_core/intervals.rs`,  
4. and per-sample altitude evaluation based on `Moon::get_horizontal(...)`.

Unlike the Sun, topocentric correctness for the Moon is not optional: lunar
parallax is on the order of **degrees** near the horizon. That means every
altitude evaluation must incorporate the observer site and Earth-rotation state
(UT1/EOP), not just an inertial ephemeris.

## Why Moon altitude is expensive

One “Moon altitude at time `t`” evaluation includes:

1) a Moon ephemeris query (geocentric position),  
2) frame transforms into an equatorial-of-date basis,  
3) topocentric translation for parallax (observer-dependent), and  
4) equatorial→horizontal conversion to extract altitude.

If you compute the Moon ephemeris from scratch each time, the ELP2000 series
summation dominates the runtime in year-long searches because the interval
finder needs thousands of evaluations (scan points + refinement iterations).

## What makes it fast in this repo

The production path used by the benchmarks relies on two ideas:

### 1) Coarse scan with refinement

The Moon period finders use a 2-hour scan step (see `SCAN_STEP` in
`../../src/calculus/lunar/altitude_periods.rs`). The scan is only responsible
for bracketing sign changes of:

```text
g(t) = altitude(t) - threshold
```

Once a bracket is found, Brent refinement produces precise boundaries without
needing a dense scan.

### 2) Make per-evaluation ephemeris cheap (Chebyshev cache)

Moon window finding is still evaluation-count heavy over 365 days. To keep it
fast, the Moon path uses a Chebyshev interpolation cache implemented in:

- `../../src/calculus/lunar/moon_cache.rs`

That cache precomputes Moon geocentric ecliptic coordinates at Chebyshev nodes
inside short segments, and evaluates intermediate times with Clenshaw
recurrence. It also caches nutation values (IAU 2000B) on a coarse grid and
interpolates them.

The net effect is that the “ephemeris part” of altitude evaluation becomes
microsecond-scale, so the remaining cost is dominated by transforms and a much
smaller number of expensive operations.

## What to optimize first

If a regression shows up, the highest-return checks are:

1. Did we accidentally increase the **scan evaluation count** (step size, extra probes)?
2. Are we still using `brent_with_values` (endpoint reuse) for refinement?
3. Is the Moon ephemeris cache being exercised (and not bypassed)?

If further speed is needed, a next step is to explicitly relax root tolerance
for planning-grade outputs using `brent_tol` (currently the Moon path uses the
default tolerance inherited from `math_core::root_finding`).

## How to run and profile

Run:

```bash
cargo bench --bench moon_altitude
```

For profiling guidance (Linux `perf` / flamegraphs), use:

- `benches/reports/profiling.md`

