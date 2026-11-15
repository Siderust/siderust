# Lunar altitude workloads

How much of a lunar event search is spent finding crossings once the position
and nutation interpolation strategy is held fixed? The
[lunar harness](../06_moon_altitude.rs) uses Roque de los Muchachos, windows
of 30, 184 and 365 days starting 2026-01-01 UTC, and `SearchOpts::default()`.

| Criterion group | Cases |
| --- | --- |
| `public_api/moon/{30d,184d,365d}` | `above_threshold/horizon`, `below_threshold/horizon`, `altitude_ranges/observation_window` (0°–30°) |
| `engines/moon/{30d,184d,365d}/above_threshold` | `chebyshev_context` and `scan_brent_baseline` at 0° |
| `engines/moon/{30d,184d,365d}/altitude_ranges` | Same engine labels for the 0°–30° band |

`chebyshev_context` measures the public default. Both engines build a
`MoonAltitudeContext` within each timed call and use its cached topocentric
signal. The default searches the residual with Chebyshev candidate discovery,
validation and local scan fallback. The internal baseline forces a uniform
10-minute scan and Brent refinement. Both include interval assembly, allocation
and destruction. Their comparison measures search strategy with shared cache
construction. It cannot quantify the benefit of the cache versus direct ELP2000.

```sh
cargo bench --locked -p siderust --bench moon_altitude --features bench-internals -- \
  'engines/moon/30d/above_threshold'
cargo test --locked -p siderust --features bench-internals --test test_event_search
```

The [event integration test](../../tests/test_event_search.rs) compares horizon
interval counts and endpoints against the cached scan for a 30-day fixture
with a `1e-5` day endpoint tolerance. The [lunar interpolation report](moon_altitude_chebyshev_cache.md)
covers direct comparisons. The current cases have no controlled timings in
this report. Earlier cached-versus-uncached values used different workloads.
