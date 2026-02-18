# Solar Altitude Bench — Performance Investigation (`find_night_periods_365day`)

> Note: This report lives under `benches/reports/` because it documents
> benchmark-driven performance work.

This report explains where time goes in the solar-night benchmarks in
`benches/solar_altitude.rs`, and which levers matter most when you want the
365‑day horizon to be fast without sacrificing physically meaningful altitude
results.

## What is being measured

The benchmark group `solar_altitude_periods` measures the end-to-end cost of
finding periods where the Sun is **below** a twilight threshold (astronomical
night uses −18°), for horizons of 1/7/30/365 days.

At the API level, the benchmark calls:

```rust
Sun.below_threshold(site, period, twilight::ASTRONOMICAL)
```

That call includes both:

1) scanning/bracketing potential crossings over the full time window, and  
2) root refinement + interval assembly to produce precise boundaries.

## Current call path (today)

At a high level, the Sun-night pipeline is:

1. `Sun.below_threshold(...)` (the altitude API)  
2. Sun-specific altitude closure in `src/calculus/solar/altitude_periods.rs`  
3. Generic scan→refine→assemble logic in `src/calculus/math_core/intervals.rs`  
4. Per-sample Sun altitude computed by `Sun::get_horizontal(...)` in
   `src/calculus/solar/sun_equations.rs`

The *important* detail is how `Sun::get_horizontal` obtains a geocentric/topocentric
Sun direction. It starts from a heliocentric “Sun at origin” position and shifts
it to the geocentric frame via a center transform:

- `cartesian::position::EclipticMeanJ2000::<U, Heliocentric>::CENTER`
- `.transform(jd)` into a geocentric equatorial J2000 position

That center transform requires Earth’s heliocentric position, which is provided
by the analytical VSOP87 backend (Earth VSOP87A). So even “Sun altitude” is
largely an **Earth ephemeris** workload.

## The cost model: evaluation count × per-evaluation cost

For this benchmark, the time spent is dominated by a single question:

How many times do we evaluate “Sun altitude at time `t`” and how expensive is
one evaluation?

### Evaluation count

The current Sun altitude period finder uses a fixed **2-hour scan step**
(`SCAN_STEP` in `src/calculus/solar/altitude_periods.rs`). That produces about
12 altitude evaluations per day for the coarse scan, plus a small number of
refinement calls near each sunrise/sunset and near each twilight crossing.

This matters more than micro-optimizations: if you double the evaluation count,
you will usually double the runtime.

### Per-evaluation cost

A single altitude evaluation (`Sun::get_horizontal`) includes:

- Earth VSOP87A evaluation (thousands of trigonometric terms)
- a chain of rotations/transforms (precession/nutation, equatorial→horizontal)
- conversion to spherical altitude

So the performance floor is set by “how expensive is Earth VSOP87A + transforms”.

## What tends to show up in profilers

On typical x86_64 CPUs, profiles of the 365‑day benchmark are usually
compute-bound and dominated by:

- trigonometric range reduction and polynomial approximations (SIMD and scalar)
- series accumulation (mul-add heavy loops)
- nutation/precession argument generation and trig

If you see `wide` SIMD trig frames dominating, that generally means “we are
paying for a lot of VSOP87 terms”. If you see `mul_add` lowering through helper
symbols, that usually means the build is not taking advantage of `+fma` on the
local CPU.

## What to optimize first (in practice)

If this benchmark regresses or becomes a bottleneck, the most reliable ordering
is:

1. **Reduce evaluation count** (scan step, reuse precomputed values, avoid extra probes).  
2. **Reduce per-evaluation work** (cache slow-changing transforms, cheaper solar model for bracketing).  
3. **Only then** consider micro-optimizations inside VSOP87 loops.

The generic interval engine already does one high-leverage thing: it reuses
scan endpoint values when refining roots (via `brent_with_values`), so a bracket
does not pay for redundant endpoint evaluations.

## How to reproduce and profile

Run the benchmark:

```bash
cargo bench --bench solar_altitude
```

Profile it (Linux, `perf`) using the workflow in `benches/reports/profiling.md`.

If you are doing local, machine-specific performance work, compiling with native
CPU features is often meaningful:

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench solar_altitude
```

That build is not portable, but it is a good way to determine whether your
bottleneck is “math kernel throughput” vs “evaluation count / algorithm”.

