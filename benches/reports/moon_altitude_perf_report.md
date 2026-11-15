# Lunar altitude workloads

The [lunar harness](../06_moon_altitude.rs) uses Roque de los Muchachos, windows
of 30, 184 and 365 days starting 2026-01-01 UTC for long searches. Searches use `SearchOpts::default()`.

| Criterion group | Cases |
| --- | --- |
| `evaluation/moon` | Cached position and altitude, Earth orientation, EOP lookup, nutation, and context construction |
| `event_windows/moon/{1d,7d,30d}` | Rise/set, horizon threshold, 0°–30° band, and a high-latitude horizon threshold in January 2024 TT |
| `event_windows/moon/shallow_horizon` | Two shallow horizon crossings near lunar lower culmination at the high-latitude site |
| `public_api/moon/{30d,184d,365d}` | `above_threshold/horizon`, `below_threshold/horizon`, `altitude_ranges/observation_window` (0°–30°) |
| `engines/moon/{30d,184d,365d}/above_threshold` | `chebyshev_context` and `scan_brent_baseline` at 0° |
| `engines/moon/{30d,184d,365d}/altitude_ranges` | Same engine labels for the 0°–30° band |

`chebyshev_context` measures the public default. Both engines build a
`MoonAltitudeContext` within each timed call and use its cached topocentric
signal. The default searches the residual with Chebyshev candidate discovery,
validation and local scan fallback. The internal baseline forces a uniform
10-minute scan and Brent refinement. Both include interval assembly, allocation
and destruction. The shared altitude evaluator interpolates the ephemeris and
evaluates Earth orientation per query.

```sh
cargo bench --locked -p siderust --bench moon_altitude --features bench-internals -- \
  'engines/moon/30d/above_threshold'
cargo bench --locked -p siderust --bench moon_altitude --features bench-internals -- \
  'evaluation/moon'
cargo bench --locked -p siderust --bench moon_altitude --no-default-features \
  --features std,embedded-compact,bench-internals -- 'event_windows/moon'
cargo test --locked -p siderust --features bench-internals --test test_event_search
```

The [event integration test](../../tests/test_event_search.rs) compares horizon
interval counts and endpoints against the cached scan for a 30-day fixture
with a `1e-5` day endpoint tolerance. The [lunar interpolation report](moon_altitude_chebyshev_cache.md)
covers direct comparisons.
