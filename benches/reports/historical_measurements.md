# Legacy altitude and ELP2000 measurements

These observations document the workloads earlier developers investigated.
Missing execution metadata prevents exact reproduction or direct comparison
with the current cases.

## Solar daily predictor

The [report at 5d0a600](https://github.com/Siderust/siderust/blob/5d0a600/benches/reports/solar_altitude_perf_report.md)
reported a 365-day 2026 window at Roque de los Muchachos, Linux 6.17, an AMD
Ryzen CPU and `cargo bench --features bench-internals`.

| Public operation | Reported time | Reported predecessor |
| --- | ---: | ---: |
| Below −18° | 471 ms | 0.85 s |
| Above 0° | 467 ms | Unknown |
| −18° to −12° band | 955 ms | 1.58 s |
| Horizon crossings | 552 ms | Unknown |

For a separate 30-day −18° case, the report described 60 crossings, about
261 precise evaluations, about 51 accepted Newton solves and zero bracket
failures or scan fallback days. The record lacks the measured baseline
revision, exact CPU, compiler, dependency state, run date, sampling overrides
and confidence intervals. These counts give a hypothesis for profiling the
daily predictor. They supply no current evaluation budget or timing baseline.

## Lunar interpolation

The [former report](https://github.com/Siderust/siderust/blob/c88f169/benches/reports/moon_altitude_chebyshev_cache.md)
recorded approximately 17.33 s for a 10-minute scan, 2.81 s
for a 2-hour uncached scan and 0.257 s for a 2-hour scan with interpolation.
All were described as a 365-day above-horizon workload at Roque starting
2026-01-01, with the historical filter `365day`. The date describes the
workload epoch. The timing date, exact revision pair, host, compiler, feature
resolution and uncertainty are absent. Scan cadence changed across rows.
Current default and forced-scan cases both build a lunar cache.

## ELP2000 before SIMD

The former [source note](https://github.com/Siderust/siderust/blob/c88f169/benches/elp2000.rs)
listed about 937 µs for J2000 and 1.66 ms for the offset case before SIMD. It
gave no machine, build or uncertainty and did not identify the measured
revision. The current harness passes the returned position to Criterion's
compiler barrier. Collect new baselines.

## Harness changes that require fresh baselines

- `vsop87` and `ephemeris_comparison` now time a fixed epoch and retain output.
- `elp2000` retains output and prepares its epoch outside timing.
- Runtime BSP cases in `de441` and `spk_perf` now use a validated fixed epoch
  instead of advancing indefinitely. Synthetic SPK case IDs are unchanged.
- `sidereal` and the sidereal cases in `de441` now obscure more inputs to the
  optimizer.
- Solar twilight batch timings now use an output-equivalent independent
  crossing workload. Use `independent_crossings_x5` in place of the old
  `independent_below_threshold_x5` ID.
