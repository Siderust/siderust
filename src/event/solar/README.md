# `event::solar`, day/night & twilight periods

Sun-specific helpers for computing **time windows** where the Sun's altitude satisfies
some condition (night, day, twilight bands, etc.).

Internal code lives in `src/event/solar/`; period finding delegates to the generic
engine in `src/event/search/`.

---

## Quick start

```rust
use siderust::bodies::Sun;
use siderust::event::altitude::{below_threshold, SearchOpts};
use siderust::event::solar::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Interval};
use qtty::{Degrees, Meter, Quantity};

let site = ObserverSite::new(
    Degrees::new(0.0),       // lon
    Degrees::new(51.4769),   // lat
    Quantity::<Meter>::new(0.0),
);

let period = Interval::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60007.0));
let nights = below_threshold(&Sun, &site, period, twilight::ASTRONOMICAL, SearchOpts::default());
```

See `examples/06_night_events.rs`.

---

## Public API

Night, day, and twilight windows use the unified altitude period functions in
`siderust::event::altitude`. Pass `&Sun` as the subject; twilight constants live in
[`twilight`] (`night_types`).

| Function | Purpose |
|----------|---------|
| [`below_threshold`] | Periods where altitude is **below** a threshold (e.g. astronomical night at −18°). |
| [`above_threshold`] | Periods where altitude is **above** a threshold (e.g. daylight above the horizon). |
| [`altitude_ranges`] | Periods where altitude lies **between** `(min, max)` (e.g. nautical twilight −18° to −12°). |

[`SearchOpts`] exposes a single public field:

- `time_tolerance` — root-refinement tolerance in days (default ~1 µs).

---

## Implementation

Solar period searches use an internal **daily predictor** that brackets rise/set
candidates per day and refines each crossing against `sun_altitude_rad`. When the
analytic model is unreliable, the engine falls back to Chebyshev polynomial fitting
or scan + Brent refinement.

For `altitude_ranges`, the band is computed as `above(min) ∩ complement(above(max))`.

---

## Twilight thresholds (`night_types`)

- `Twilight` enum (`Civil`, `Nautical`, `Astronomical`, `Horizon`, `ApparentHorizon`)
  — `Degrees::from(Twilight)` converts to a threshold angle.
- `twilight::*` constants (Sun **center** altitude):
  - `CIVIL = -6°`, `NAUTICAL = -12°`, `ASTRONOMICAL = -18°`
  - `HORIZON = 0°` (geometric)
  - `APPARENT_HORIZON = -0.833°` (rule-of-thumb refraction + Sun semi-diameter)

Altitude is **geometric** (no refraction). Treat `APPARENT_HORIZON` as a convenient
approximation, not a site-specific refraction model.

---

## Accuracy & edge cases

- **Numerical precision**: Brent refinement with `time_tolerance` typically yields
  sub-millisecond event times.
- **Physical limits**: VSOP87 ephemeris, threshold choice (center vs limb), and
  ignored atmospheric refraction dominate real-world error.
- **Polar day/night**: empty results are normal when the Sun never crosses the
  requested threshold.
- **Grazing events**: tangent crossings near polar boundaries may be missed by
  sign-change bracketing.

Benchmarks: `cargo bench --bench solar_altitude`.

[`below_threshold`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.below_threshold.html
[`above_threshold`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.above_threshold.html
[`altitude_ranges`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.altitude_ranges.html
[`twilight`]: https://docs.rs/siderust/latest/siderust/event/solar/night_types/index.html
[`SearchOpts`]: https://docs.rs/siderust/latest/siderust/event/altitude/struct.SearchOpts.html
