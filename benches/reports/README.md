# Benchmark reports

Start with the [benchmark catalog](../README.md) for Cargo targets, features and
case selection.

## Plan, reproduce or diagnose an experiment

- [Measurement and reproducibility](methodology.md) covers timing boundaries,
  baseline capture, uncertainty, scientific validation and reporting.
- [Profiling](profiling.md) shows how to isolate a Criterion case and inspect
  CPU stacks, hardware counters and heap allocation on Linux.

## Investigate a workload

- [Coverage and validation](coverage.md) maps coordinates, time, ephemerides,
  observation events, SPK and propagation workloads to their reference tests
  and measurement gaps.
- [Solar altitude](solar_altitude_perf_report.md) describes the daily predictor,
  search baselines and the five-threshold crossing comparison.
- [Lunar altitude](moon_altitude_perf_report.md) explains what the public and
  internal cases measure, including cache construction in both paths.
- [Lunar interpolation](moon_altitude_chebyshev_cache.md) describes the position
  and nutation caches and the tests that bound sampled interpolation error.
- [Root refinement](fast_root_finding_brent.md) explains bracket, endpoint and
  crossing-direction reuse, and where function counts can be measured.
