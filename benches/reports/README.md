# Benchmark Reports (Performance / Profiling / Investigations)

This folder contains engineering reports and performance investigations for the
Criterion benchmarks under `benches/`.

These documents are **not** treated as stable architecture docs. Architectural
decisions and guarantees live under `doc/`.

## Index

If you are chasing a regression, start with `profiling.md` to set up a reliable
profiling loop, then pick the subsystem report that matches the benchmark.

- `benches/reports/profiling.md`
  How to profile Criterion benchmarks (symbols, `cargo flamegraph`, noise
  reduction, and common pitfalls).

- `benches/reports/solar_altitude_perf_report.md`
  Where time goes in `benches/solar_altitude.rs` (Sun altitude windows), and
  the levers that matter most (evaluation count vs per-evaluation cost).

- `benches/reports/moon_altitude_perf_report.md`
  Performance characteristics of Moon altitude window finding
  (`benches/moon_altitude.rs`), including why topocentric parallax makes this
  path heavier than the Sun case.

- `benches/reports/moon_altitude_chebyshev_cache.md`
  Design and rationale for the Moon Chebyshev ephemeris cache and nutation
  interpolation used to accelerate long-horizon searches.

- `benches/reports/fast_root_finding_brent.md`
  How Brent root refinement and endpoint reuse reduce expensive function
  evaluations in scan→refine pipelines.
