# Benchmark Reports (Performance / Profiling / Investigations)

This folder contains engineering reports and performance investigations for the
Criterion benchmarks under `benches/`.

These documents are **not** treated as stable architecture docs. Architectural
decisions and guarantees live under `doc/`.

## Benchmark taxonomy

Benchmarks are organized around two concepts:

- **`public_api/{solar,moon,star}/{30d,184d,365d}`** — stable Option A API calls
  (`below_threshold`, `above_threshold`, `altitude_ranges`, `crossings`).
  These reflect what users measure.
- **`engines/{solar,moon}/{30d,184d,365d}`** — internal engine comparisons
  (daily predictor vs Chebyshev baseline vs scan+Brent). Require `bench-internals`.
- **`altitude/single_eval`** and **`public_api/comparison/*`** — cross-body
  comparisons in `altitude_comparison`.

## How to run

```bash
# All benches (public_api + engines, requires bench-internals)
cargo bench --features bench-internals

# Gated per-body benches
cargo bench --features bench-internals --bench solar_altitude
cargo bench --features bench-internals --bench moon_altitude

# Non-gated benches (star, altitude_comparison, vsop87, …)
cargo bench --bench star_altitude
cargo bench --bench altitude_comparison
```

## Index

If you are chasing a regression, start with `profiling.md` to set up a reliable
profiling loop, then pick the subsystem report that matches the benchmark.

- `benches/reports/profiling.md`
  How to profile Criterion benchmarks (symbols, `cargo flamegraph`, noise
  reduction, and common pitfalls).

- `benches/reports/solar_altitude_perf_report.md`
  Groups `public_api/solar` and `engines/solar`: solar daily predictor vs
  Chebyshev-first and scan+Brent baselines; batch altitude-range path.
  Requires `bench-internals`.

- `benches/reports/moon_altitude_perf_report.md`
  Groups `public_api/moon` and `engines/moon`: `MoonAltitudeContext` +
  Chebyshev-first vs scan+Brent baseline. Requires `bench-internals`.

- `benches/reports/moon_altitude_chebyshev_cache.md`
  Design and rationale for the Moon Chebyshev ephemeris cache and nutation
  interpolation used to accelerate long-horizon searches.

- `benches/reports/fast_root_finding_brent.md`
  How Brent root refinement and endpoint reuse reduce expensive function
  evaluations in scan→refine pipelines.
