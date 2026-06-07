# `event::solar`, day/night & twilight periods

This folder contains Sun-specific helpers for computing **time windows** where the Sun's
altitude satisfies some condition (night, day, twilight bands, etc.).

The internal implementation lives in:

- `src/event/solar/` (solar altitude closures, `solar_*_impl`)
- `src/event/solar/night_types.rs`

Under the hood they delegate to the generic numerical engine in:

- `src/event/search/intervals.rs` (scan + Brent refinement + crossing classification)
- `src/event/search/root_finding.rs` (Brent, bisection)

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

See runnable examples:

- `examples/06_night_events.rs`

---

## Public API (what each function does)

Night, day, and twilight windows are expressed through the unified altitude period
functions in `siderust::event::altitude`:

- [`below_threshold`] — periods where altitude is **below** a threshold (e.g. astronomical night at −18°).
- [`above_threshold`] — periods where altitude is **above** a threshold (e.g. daylight above the horizon).
- [`altitude_ranges`] — periods where altitude lies **between** `(min, max)` (e.g. nautical twilight band −18° to −12°).

Pass `&Sun` as the subject. Twilight constants live in [`twilight`] (`night_types`).

### Altitude evaluator (internal)

- `sun_altitude_rad(jd, &site) -> Quantity<Radian>`
  - Computes the Sun's **topocentric geometric altitude** (no refraction) at a given `JulianDate`.
  - Used by the search engine as the model function to be thresholded.

### Search options

[`SearchOpts`] controls scan step and refinement behaviour. The default engine fits
Chebyshev polynomials per crossing segment and falls back to scan+Brent when needed.
A 2-hour scan step (`scan_step_days`) yields ~12 VSOP87 evaluations per day and is
fast enough for multi-year sweeps.

---

## Twilight thresholds (`night_types`)

`src/event/solar/night_types.rs` provides:

- `Twilight` enum (`Civil`, `Nautical`, `Astronomical`, `Horizon`, `ApparentHorizon`)
  - `Degrees::from(Twilight)` converts to a threshold angle.
- `twilight::*` constants (all are **Sun center altitude** thresholds):
  - `CIVIL = -6°`, `NAUTICAL = -12°`, `ASTRONOMICAL = -18°`
  - `HORIZON = 0°` (geometric)
  - `APPARENT_HORIZON = -0.833°` (rule-of-thumb refraction + Sun semi-diameter)

Important: the altitude computed here is **geometric**. If you need site-specific refraction
models, treat `APPARENT_HORIZON` as a convenient approximation rather than a guarantee.

---

## Algorithm

All period searches share the same pipeline from `event::search::intervals`:

1. **Coarse scan** at fixed time steps to detect sign changes in
   `altitude(t) − threshold`.
2. **Brent refinement** of each bracket to ~1 µs precision.
3. **Crossing classification** (rising vs setting).
4. **Interval assembly** from consecutive crossing pairs.

For `altitude_ranges` queries, the band is computed as
`above(min) ∩ complement(above(max))`, two interval-algebra passes at near-zero cost.

### Scan step choice

| Mode              | Step     | Evals/day | Use case                    |
|-------------------|----------|-----------|-----------------------------|
| Default           | adaptive | varies    | Production (Chebyshev-first)|
| `scan_step_days`  | 2 hours  | ~12       | Scan+Brent baseline         |
| finer scan step   | 10 min   | ~144      | Validation / comparison     |

The 2-hour step is safe because the shortest daylight arc at 65° latitude is ~5 h,
so every sunrise/sunset is guaranteed to fall within at least one bracket.

---

## Performance

The dominant cost is Sun position evaluation (VSOP87 + transforms + trig).

### Measuring on your machine

Benchmarks live in `benches/solar_altitude.rs`:

```bash
cargo bench --bench solar_altitude
```

---

## Accuracy

### Event-time precision (numerics)

Both scan steps refine crossings using Brent's method with tight tolerances, so
**numerical** root timing precision is typically far below 1 millisecond.

### Physical/model accuracy (dominant in practice)

Your limiting factors are usually:

- the underlying solar/earth ephemeris model used by `Sun::get_horizontal(...)`,
- choice of twilight threshold (Sun center vs limb),
- and atmospheric refraction (ignored by the geometric altitude; approximated only via
  `twilight::APPARENT_HORIZON`).

### Edge cases to be aware of

- **Polar day / polar night**: it's normal to get an empty result for "night" during summer
  at high latitudes (Sun never goes below the requested threshold).
- **Grazing events**: when the Sun *just touches* a threshold (double root / tangent),
  sign-change bracketing can miss the event. This happens near the boundary of polar day/night.

[`below_threshold`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.below_threshold.html
[`above_threshold`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.above_threshold.html
[`altitude_ranges`]: https://docs.rs/siderust/latest/siderust/event/altitude/fn.altitude_ranges.html
[`twilight`]: https://docs.rs/siderust/latest/siderust/event/solar/night_types/index.html
[`SearchOpts`]: https://docs.rs/siderust/latest/siderust/event/altitude/struct.SearchOpts.html
