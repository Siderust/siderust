# Solar altitude workloads

The [solar harness](../09_solar_altitude.rs) measures complete twilight and horizon
searches at Roque de los Muchachos. All groups use 30-, 184- and 365-day windows
starting 2026-01-01 UTC with `SearchOpts::default()`. In the group names below
`{N}` means `30d`, `184d` or `365d`.

| Criterion group | Operation and cases |
| --- | --- |
| `public_api/solar/{N}` | `below_threshold/astro_night` at −18°, `above_threshold/horizon` at 0°, `altitude_ranges/astro_to_nautical` at −18° to −12°, and `crossings/horizon` at 0° |
| `engines/solar/{N}/below_threshold` | −18° intervals with `daily_predictor`, `chebyshev_baseline` or `scan_brent_baseline` |
| `engines/solar/{N}/altitude_ranges` | −18° to −12° intervals with `daily_batch`, `chebyshev_baseline` or `scan_brent_baseline` |
| `engines/solar/{N}/twilight_profile` | Labeled crossings grouped by threshold (0°, −0.833°, −6°, −12°, −18°) with `daily_batch_all_thresholds` or `independent_crossings_x5` |

The [daily predictor](../../src/event/solar/daily_events.rs) forms a solar state
per day, proposes crossing times from hour-angle geometry and uses a cheaper
altitude model to guide refinement. It checks roots against the precise
altitude residual using Newton/secant and bracketed refinement as applicable.
Difficult segments can fall back to the generic crossing search. The range
implementation shares the daily work across both thresholds.
[Period assembly](../../src/event/solar/periods.rs), context creation and result
allocation/destruction are included in the interval cases.

The internal Chebyshev baseline disables the daily predictor. The scan baseline
uses a forced 10-minute step. The twilight batch case uses one pass over the
days. Its independent baseline performs five separate single-threshold crossing
searches with the same output structure. Before timing, the harness checks
counts, directions and agreement within the default time tolerance at each
threshold. The former `independent_below_threshold_x5` case produced assembled
intervals, so its measurements require a fresh baseline.

```sh
cargo bench --locked -p siderust --bench solar_altitude --features bench-internals -- \
  'engines/solar/30d/twilight_profile'
cargo test --locked -p siderust --features bench-internals --test test_event_search
cargo test --locked -p siderust --test test_astro_nights_roque_2026
```

[Event integration tests](../../tests/test_event_search.rs) compare selected
interval counts and endpoints with the internal scan and check window bounds.
The daily engine has additional direction, partial-window and polar cases.
The [Roque astronomical night test](../../tests/test_astro_nights_roque_2026.rs)
compares interval boundaries against a checked-in fixture with a three-second
tolerance.
The mid-latitude benchmark primarily samples the normal predictor path.
Fallback frequency and precise evaluation cost vary by site, threshold and
geometry. Current timing claims require a new experiment. The previous timings
and diagnostic counts appear in the [historical record](historical_measurements.md).
