# `calculus::solar` — day/night & twilight periods

This folder contains Sun-specific helpers for computing **time windows** where the Sun’s
altitude satisfies some condition (night, day, twilight bands, etc.).

The public entrypoints live in:

- `src/calculus/solar/altitude_periods.rs`
- `src/calculus/solar/night_types.rs`

Under the hood they reuse generic infrastructure from:

- `src/calculus/events/altitude_periods.rs` (scan + refine)
- `src/calculus/events/find_dynamic_extremas.rs` (culminations / extrema)
- `src/calculus/root_finding/` (Brent, Newton, bisection)

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
  - Computes the Sun’s **topocentric geometric altitude** (no refraction) at a given `JulianDate`.
  - Implementation: calls `Sun::get_horizontal::<AstronomicalUnit>(jd, site)` and reads `.alt()`.
  - Used by all the period finders as the “truth” function to be thresholded.

### Period finders (recommended)

All return `Option<Vec<Period<ModifiedJulianDate>>>`:

- `Some(periods)` when at least one interval satisfies the condition.
- `Some(vec![period])` when the whole input interval satisfies the condition.
- `None` when the condition is never satisfied.

Functions:

- `find_night_periods(site, period, twilight)`
  - Night = Sun altitude **below** `twilight` (e.g. `twilight::ASTRONOMICAL` = -18°).
- `find_day_periods(site, period, twilight)`
  - Day = Sun altitude **above** `twilight` (often `twilight::HORIZON` or `twilight::APPARENT_HORIZON`).
- `find_sun_range_periods(site, period, (min, max))`
  - Finds windows where Sun altitude is **between** `(min, max)` (inclusive).
  - Useful for twilight “bands”, e.g. nautical twilight could be `(-18°, -12°)`.

Implementation note: these three are thin wrappers around
`find_sun_altitude_periods_via_culminations(...)` (see below).

### Period finders (scan-based, mostly for comparison)

- `find_night_periods_scan(...)`
- `find_day_periods_scan(...)`
- `find_sun_range_periods_scan(...)`

These delegate to the **generic** routine:
`calculus::events::altitude_periods::find_altitude_periods(...)`.

Prefer the non-`_scan` variants unless you’re:

- comparing algorithms/behavior,
- using a custom altitude function (non-solar) via the generic API, or
- validating assumptions near edge cases.

### Low-level building block

- `find_sun_altitude_periods_via_culminations(site, period, condition)`
  - “Core” Sun-specific routine that supports `Below`, `Above`, and `Between` via
    `calculus::events::altitude_periods::AltitudeCondition`.

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

## Algorithms (what’s different, and why)

Both approaches ultimately:

1. **bracket** threshold crossings, then
2. **refine** each crossing time with Brent’s method, then
3. **classify** crossings (enter vs exit), and
4. **pair** crossings into contiguous `Period`s.

The difference is how they obtain good brackets with *few* expensive Sun position evaluations.

### 1) Generic scan + refine (`find_*_scan`, and `calculus::events::find_altitude_periods`)

File: `src/calculus/events/altitude_periods.rs`

How it works:

- Sample the altitude function in fixed steps (`SCAN_STEP = 10 minutes`).
- For each boundary (1 for `Below/Above`, 2 for `Between`), watch for `f(t)` sign changes.
- When a sign change is found, refine the root with
  `calculus::root_finding::find_crossing_brent_with_values(...)`.

Pros:

- generic: works for any altitude function you can provide;
- simple mental model;
- refinement is high precision once a crossing is bracketed.

Cons:

- performance scales with `(interval length) / SCAN_STEP`;
- for `Between {min, max}` it scans *twice* (two boundaries);
- any scan-based method can miss events that do **not** produce a sign change
  (a grazing/tangent contact with the threshold).

### 2) Culmination-partition + refine (`find_*` and `find_sun_altitude_periods_via_culminations`)

File: `src/calculus/solar/altitude_periods.rs`

How it works:

1. Find the Sun’s **upper and lower culminations** (meridian crossings) over the interval
   using `calculus::events::find_dynamic_extremas(...)`.
2. Use `[start, culm..., end]` as “key times” and iterate over adjacent pairs.
   Between an upper and a lower culmination the Sun’s altitude is (to a very good
   approximation) **monotonic**, so there can be at most one crossing per boundary.
3. For each boundary and each segment:
   - evaluate altitude at the segment endpoints;
   - refine with Brent only when there is a sign change.

Why this helps:

- For the Sun you expect roughly **two** altitude extrema per day (upper/lower transit),
  so the number of segments grows like `~2 * days`, not `~144 * days`.

Trade-offs / caveats:

- There is extra work up front to find culminations (itself implemented as a scan+refine
  on hour angle in `find_dynamic_extremas`).
- Like the scan method, it primarily detects **sign-changing** crossings; purely tangential
  contacts can be missed unless they occur exactly at segment boundaries.

---

## Performance comparison (rule of thumb)

The dominant cost is usually “Sun position evaluation” (VSOP87 + transforms + trig), not
allocation or bookkeeping. See `doc/solar_altitude_perf_report.md` for a deeper profile.

### Scaling

Let:

- `D` = days in the search period
- `B` = number of boundaries (`1` for `Below/Above`, `2` for `Between`)

Then the *coarse bracketing* evaluation counts scale approximately as:

- **Scan-based**: `~ 144 * D * B` altitude evaluations (`10 min` step)
- **Culmination-based**:
  - `~ 72 * D` evaluations to locate culminations (`20 min` step in `find_dynamic_extremas`)
  - plus `~ O(2 * D * B)` endpoint checks for threshold crossings
  - plus refinement iterations near each actual crossing (typically a small constant factor)

Practical implication:

- For short horizons (≈ 1 day) both can be similar (fixed overheads matter).
- For longer horizons and/or `Between {min, max}` conditions, the culmination-based
  method typically wins because it avoids doing a fixed 10-minute scan for *each* boundary.

### Measuring on your machine

Benchmarks live in `benches/solar_altitude.rs`:

```bash
cargo bench --bench solar_altitude
```

To focus on a subset, pass a filter string:

```bash
cargo bench --bench solar_altitude find_night_periods_scan_7day
```

Example (AMD Ryzen 7 5700X, 2026-02-06; Criterion shortened run):

- `find_night_periods_7day` (culminations): ~14.7 ms
- `find_night_periods_scan_7day` (scan): ~23.5 ms
- speedup: ~1.6×

---

## Accuracy comparison (what “accuracy” means here)

### Event-time precision (numerics)

Both methods refine crossings using Brent’s method with tight tolerances, so **numerical**
root timing precision is typically far below 1 millisecond.

If both algorithms detect the same crossings, they should agree closely because they refine
against the same `sun_altitude_rad` function.

### Physical/model accuracy (dominant in practice)

Your limiting factors are usually:

- the underlying solar/earth ephemeris model used by `Sun::get_horizontal(...)`,
- choice of twilight threshold (Sun center vs limb),
- and atmospheric refraction (ignored by the geometric altitude; approximated only via
  `twilight::APPARENT_HORIZON`).

### Edge cases to be aware of

- **Polar day / polar night**: it’s normal to get `None` for “night” during summer at high
  latitudes (Sun never goes below the requested threshold) or `Some(vec![period])` when
  the Sun never goes above it.
- **Grazing events**: when the Sun *just touches* a threshold (double root / tangent),
  sign-change bracketing can miss the event. This happens near the boundary of polar day/night.
