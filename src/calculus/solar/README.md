# `calculus::solar` — day/night & twilight periods

This folder contains Sun-specific helpers for computing **time windows** where the Sun's
altitude satisfies some condition (night, day, twilight bands, etc.).

The public entrypoints live in:

- `src/calculus/solar/altitude_periods.rs`
- `src/calculus/solar/night_types.rs`

Under the hood they delegate to the generic numerical engine in:

- `src/calculus/math_core/intervals.rs` (scan + Brent refinement + crossing classification)
- `src/calculus/math_core/root_finding.rs` (Brent, bisection)

---

## Quick start

```rust
use siderust::calculus::solar::{find_night_periods, twilight};
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Period};
use qtty::{Degrees, Meter, Quantity};

let site = ObserverSite::new(
    Degrees::new(0.0),       // lon
    Degrees::new(51.4769),   // lat
    Quantity::<Meter>::new(0.0),
);

let period = Period::new(ModifiedJulianDate::new(60000.0), ModifiedJulianDate::new(60007.0));
let nights = find_night_periods(site, period, twilight::ASTRONOMICAL);
```

See runnable examples:

- `examples/astronomical_night.rs`
- `examples/find_night_periods_365day.rs`
- `examples/solar_altitude_culminations.rs`

---

## Public API (what each function does)

### Altitude evaluator

- `sun_altitude_rad(jd, &site) -> Quantity<Radian>`
  - Computes the Sun's **topocentric geometric altitude** (no refraction) at a given `JulianDate`.
  - Implementation: calls `Sun::get_horizontal::<AstronomicalUnit>(jd, site)` and reads `.alt()`.
  - Used by all the period finders as the "truth" function to be thresholded.

### Period finders (recommended — 2-hour scan step)

Functions:

- `find_day_periods(site, period, twilight)`
  - Day = Sun altitude **above** `twilight` (often `twilight::HORIZON` or `twilight::APPARENT_HORIZON`).
- `find_night_periods(site, period, twilight)`
  - Night = Sun altitude **below** `twilight` (e.g. `twilight::ASTRONOMICAL` = -18°).
- `find_sun_range_periods(site, period, (min, max))`
  - Finds windows where Sun altitude is **between** `(min, max)` (inclusive).
  - Useful for twilight "bands", e.g. nautical twilight could be `(-18°, -12°)`.

Implementation: all three delegate to `math_core::intervals` with a **2-hour** scan step,
yielding ~12 VSOP87 evaluations per day — fast enough for multi-year sweeps.

### Period finders (scan-based — 10-minute step, for comparison)

- `find_day_periods_scan(...)`
- `find_night_periods_scan(...)`
- `find_sun_range_periods_scan(...)`

These use a finer **10-minute** scan step via `math_core::intervals`.
Prefer the 2-hour variants unless you're comparing algorithms or validating edge cases.

---

## Twilight thresholds (`night_types`)

`src/calculus/solar/night_types.rs` provides:

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

All period finders share the same pipeline from `math_core::intervals`:

1. **Coarse scan** at fixed time steps to detect sign changes in
   `altitude(t) − threshold`.
2. **Brent refinement** of each bracket to ~1 µs precision.
3. **Crossing classification** (rising vs setting).
4. **Interval assembly** from consecutive crossing pairs.

For `Between {min, max}` queries, the range is computed as
`above(min) ∩ complement(above(max))` — two interval-algebra passes at near-zero cost.

### Scan step choice

| Variant      | Step     | Evals/day | Use case                    |
|--------------|----------|-----------|-----------------------------|
| `find_*`     | 2 hours  | ~12       | Production (fast)           |
| `find_*_scan`| 10 min   | ~144      | Validation / comparison     |

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
