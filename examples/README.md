# Siderust Examples

Runnable examples organized by theme. Each demonstrates a specific aspect of
the library — from basic coordinate algebra to full observing-session planners.

> **Tip**: Run any example with `cargo run --example <name>`.
> Feature-gated examples require `--features <feat>` (noted below).

---

## Getting Started

These introduce the core types and concepts you will use everywhere else.

### 1. Basic Coordinates (`basic_coordinates.rs`)
Cartesian and spherical coordinate types, reference frames (EclipticMeanJ2000,
EquatorialMeanJ2000, ICRS) and reference centers (Helio-, Geo-, Barycentric).

```bash
cargo run --example basic_coordinates
```

### 2. Coordinate Transformations (`coordinate_transformations.rs`)
Frame transforms (Ecliptic ↔ Equatorial ↔ ICRS), center transforms
(Helio ↔ Geo ↔ Bary), combined transforms and round-trip verification.

```bash
cargo run --example coordinate_transformations
```

### 3. Frame Conversion Matrix (`all_frame_conversions.rs`)
Prints every currently-supported **frame rotation** and a round-trip error
metric (includes time-dependent specialty conversions like Horizontal).

```bash
cargo run --example all_frame_conversions
```

### 4. Center Conversion Matrix (`all_center_conversions.rs`)
Prints every currently-supported **center shift** and a round-trip error metric,
including parameterized centers (`Topocentric`, `Bodycentric`).

```bash
cargo run --example all_center_conversions
```

### 5. Time Periods (`time_periods.rs`)
`Period<T>` for JulianDate, MJD, UTC DateTime — conversions, arithmetic,
and duration queries.

```bash
cargo run --example time_periods
```

---

## Observational Astronomy

These focus on altitude calculations, night-period finding, and observing
session planning at real observatory sites.

### 6. Astronomical Night (`astronomical_night.rs`)
CLI tool: finds astronomical night periods (Sun < −18°) for 7 days from
a given date/location. Accepts CLI arguments with defaults to Greenwich.

```bash
cargo run --example astronomical_night
cargo run --example astronomical_night -- 2026-06-21 51.4769 -0.0005 80
```

### 7. Night Periods — Full Year (`find_night_periods_365day.rs`)
Computes all astronomical night windows for a 365-day horizon at Roque de
los Muchachos. Accepts an optional start-date CLI argument.

```bash
cargo run --example find_night_periods_365day
cargo run --example find_night_periods_365day -- 2026-01-01
```

### 8. Night Quality Scoring (`night_quality_scoring.rs`)
Practical observing planner: scores 30 consecutive nights at Mauna Kea
based on darkness duration and Moon interference.

```bash
cargo run --example night_quality_scoring
```

### 9. Star Observability (`star_observability.rs`)
Observing planner: 6 stars at Greenwich — visibility windows during
astronomical night, peak altitudes, and observing strategy.

```bash
cargo run --example star_observability
```

### 10. Altitude Periods API (`altitude_periods_trait.rs`)
Comprehensive tour of the unified `AltitudePeriodsProvider` trait:
astronomical nights, star visibility, custom ICRS targets, Moon range
queries, circumpolar detection, twilight bands, single-point altitude.

```bash
cargo run --example altitude_periods_trait
```

### 11. Generic Body Comparison (`compare_sun_moon_star.rs`)
Analyzes Sun, Moon, and a star using the **same generic function**,
showcasing the polymorphic trait-based API.

```bash
cargo run --example compare_sun_moon_star
```

---

## Solar System & Coordinates

Working with planets, body-centric views, observer locations, and
the full solar-system coordinate machinery.

### 12. Solar System Bodies (`solar_system_example.rs`)
All 8 planets at J2000, inter-planetary distances, geocentric positions,
time evolution, barycentric vs heliocentric comparison.

```bash
cargo run --example solar_system_example
```

### 13. Observer Coordinates (`observer_coordinates.rs`)
Defining observer locations, horizontal frames, topocentric transforms,
multi-observer parallax, and altitude effects.

```bash
cargo run --example observer_coordinates
```

### 14. Body-Centric Coordinates (`bodycentric_coordinates.rs`)
ISS-centric, Mars-centric, Venus-centric views, round-trip transforms,
and directions as free vectors.

```bash
cargo run --example bodycentric_coordinates
```

---

## Targets & Proper Motion

### 15. Targets + Proper Motion (`targets_proper_motion.rs`)
`Target<T>` as a (position, epoch, proper-motion) bundle, plus applying a
Gaia/Hipparcos-style proper motion model from J2000.0 to a future epoch.

```bash
cargo run --example targets_proper_motion
```

---

## Ephemeris Backends

### 16. JPL Precise Ephemeris (`jpl_precise_ephemeris.rs`)
Compares VSOP87/ELP2000 with JPL DE440 and/or DE441 — Earth and Moon
positions, velocities, and precision differences side-by-side.

Requires at least one JPL feature:
```bash
cargo run --example jpl_precise_ephemeris --features de440
cargo run --example jpl_precise_ephemeris --features de441
cargo run --example jpl_precise_ephemeris --features de440,de441
```

Fast/offline loop (skips runtime DE calls and avoids large downloads):
```bash
SIDERUST_JPL_STUB=all cargo run --example jpl_precise_ephemeris --features de440,de441
```

---

## Serialization

### 17. Serde Serialization (`serde_serialization.rs`)
JSON round-trips for JulianDate, cartesian positions, spherical
directions (with frame-specific field names), complex structs,
file I/O, and collections.

Requires the `serde` feature:
```bash
cargo run --example serde_serialization --features serde
```

---

## Numerical Methods

### 18. Kepler Orbit (`kepler_orbit.rs`)
Solves Kepler’s equation and propagates a simple orbit with `Orbit::kepler_position`.

```bash
cargo run --example kepler_orbit
```

### 19. Brent Root Finding (`brent_root_finding.rs`)
Unit-safe root finding with the shared numerical engine (`math_core`).

```bash
cargo run --example brent_root_finding
```

---

## Key Concepts

### Reference Centers
| Center | Origin |
|--------|--------|
| `Barycentric` | Solar system barycenter |
| `Heliocentric` | Center of the Sun |
| `Geocentric` | Center of the Earth |
| `Topocentric` | Observer's location on Earth's surface |
| `Bodycentric` | Any orbiting body (parameterized) |

### Reference Frames
| Frame | Orientation |
|-------|------------|
| `ICRS` | Fixed to distant quasars (IAU standard) |
| `EclipticMeanJ2000` | Mean ecliptic of J2000.0 |
| `EquatorialMeanJ2000` | Mean equator/equinox of J2000.0 |
| `EquatorialMeanOfDate` | Precession applied |
| `EquatorialTrueOfDate` | Precession + nutation applied |
| `Horizontal` | Local horizon (Alt/Az, IAU convention) |
| `Galactic` | Galactic coordinate system |
| `ECEF` | Earth-Centered Earth-Fixed |

### Type Safety

All coordinate types are parameterized by **Center** (`C`), **Frame** (`F`),
and **Unit** (`U`). The compiler enforces that you never mix incompatible
coordinate systems:

```rust
let pos1: Position<Geocentric, EquatorialMeanJ2000, Km> = /* ... */;
let pos2: Position<Geocentric, EquatorialMeanJ2000, Km> = /* ... */;
let d = pos1.distance_to(&pos2);  // ✓ compiles

let pos3: Position<Heliocentric, EclipticMeanJ2000, Au> = /* ... */;
// pos1.distance_to(&pos3);  // ✗ compile error — different types!
```

---

## Running All Examples

```bash
# Getting started
cargo run --example basic_coordinates
cargo run --example coordinate_transformations
cargo run --example all_frame_conversions
cargo run --example all_center_conversions
cargo run --example time_periods

# Observational astronomy
cargo run --example astronomical_night
cargo run --example find_night_periods_365day
cargo run --example night_quality_scoring
cargo run --example star_observability
cargo run --example altitude_periods_trait
cargo run --example compare_sun_moon_star

# Solar system & coordinates
cargo run --example solar_system_example
cargo run --example observer_coordinates
cargo run --example bodycentric_coordinates

# Targets & proper motion
cargo run --example targets_proper_motion

# Ephemeris backends (feature-gated)
cargo run --example jpl_precise_ephemeris --features de440
cargo run --example jpl_precise_ephemeris --features de441

# Serialization (feature-gated)
cargo run --example serde_serialization --features serde

# Numerical methods
cargo run --example kepler_orbit
cargo run --example brent_root_finding
```
