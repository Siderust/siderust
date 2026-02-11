# Siderust Coordinate System Examples

This directory contains examples demonstrating how to use the siderust coordinate system and transformations.

## Available Examples

### 1. Basic Coordinates (`basic_coordinates.rs`)
Introduction to creating and using different coordinate systems:
- Cartesian coordinates (Position and Direction)
- Spherical coordinates
- Different reference frames (Ecliptic, EquatorialMeanJ2000, ICRS)
- Different reference centers (Heliocentric, Geocentric, Barycentric)

**Run with:**
```bash
cargo run --example basic_coordinates
```

### 2. Coordinate Transformations (`coordinate_transformations.rs`)
Demonstrates transformations between different coordinate systems:
- Frame transformations (Ecliptic ↔ EquatorialMeanJ2000 ↔ ICRS)
- Center transformations (Heliocentric ↔ Geocentric ↔ Barycentric)
- Combined transformations
- Time-dependent transformations

**Run with:**
```bash
cargo run --example coordinate_transformations
```

### 3. Body-Centric Coordinates (`bodycentric_coordinates.rs`)
Shows how to use the body-centric coordinate system for viewing from arbitrary orbiting bodies:
- Satellite-centric coordinates (ISS example)
- Planet-centric coordinates (Mars view of Earth)
- Comet-centric coordinates
- Round-trip transformations

**Run with:**
```bash
cargo run --example bodycentric_coordinates
```

### 4. Observer-Based Coordinates (`observer_coordinates.rs`)
Demonstrates topocentric (observer-based) coordinates:
- Defining observer locations
- Horizontal coordinate system (altitude/azimuth)
- Converting between geocentric and topocentric
- Real observatory examples

**Run with:**
```bash
cargo run --example observer_coordinates
```

### 5. Solar System Bodies (`solar_system_example.rs`)
Working with solar system bodies and their positions:
- Computing planetary positions
- Planet-to-planet views
- Using VSOP87 ephemerides
- Orbital mechanics

**Run with:**
```bash
cargo run --example solar_system_example
```

### 6. DE440 Precise Ephemeris (`de440_precise_ephemeris.rs`)
Demonstrates high-precision ephemeris calculations using JPL DE440:
- Computing Sun, Earth, and Moon positions with meter-level precision
- Working with different coordinate centers (barycentric, heliocentric, geocentric)
- Computing velocities for aberration corrections
- Comparing DE440 accuracy vs VSOP87/ELP2000
- Time series analysis of Earth-Sun distance
- Understanding when to use DE440 vs analytical theories

**Requires the `de440` feature:**
```bash
cargo run --example de440_precise_ephemeris --features de440
```

### 7. Serialization and Deserialization (`serde_serialization.rs`)
Learn how to serialize and deserialize siderust types:
- Julian dates and time types
- Cartesian coordinates (positions and directions)
- Spherical coordinates
- Complex observation data structures
- Saving/loading astronomical data to/from JSON files
- Working with collections of astronomical data

**Run with:**
```bash
cargo run --example serde_serialization --features serde
```

### 7. Solar Altitude Periods (`solar_altitude_culminations.rs`)
Computes a full year (2026) of astronomical night intervals at Roque de los
Muchachos using `find_night_periods`.

**Run with:**
```bash
cargo run --example solar_altitude_culminations
```

### 8. Altitude Periods API (`altitude_periods_trait.rs`)
Comprehensive demonstration of the unified `AltitudePeriodsProvider` trait:
- Finding astronomical night periods (Sun below -18°)
- Computing star visibility windows above specific altitudes
- Using `direction::ICRS` for custom targets
- Moon altitude range queries
- Detecting circumpolar stars
- Twilight band calculations
- Single-point altitude evaluation

**Run with:**
```bash
cargo run --example altitude_periods_trait
```

### 9. Star Observability Planner (`star_observability.rs`)
Practical observing planner that analyzes star visibility:
- Finds best observing windows for multiple stars
- Considers astronomical darkness and altitude constraints
- Computes peak altitudes during observation window
- Provides observing strategy recommendations
- Identifies circumpolar vs. rising/setting stars

**Run with:**
```bash
cargo run --example star_observability
```

### 10. Generic Body Comparison (`compare_sun_moon_star.rs`)
Shows the power of the trait-based API by analyzing Sun, Moon, and stars
using the **same generic function**:
- Demonstrates unified interface across different body types
- Compares Sun (VSOP87), Moon (ELP2000), and Stars (analytical)
- Shows `direction::ICRS` as the lightweight path
- Illustrates how the trait abstracts implementation details

**Run with:**
```bash
cargo run --example compare_sun_moon_star
```

### 12. Night Quality Scoring (`night_quality_scoring.rs`)
Practical tool for planning observing sessions over a month:
- Scores each night based on darkness duration and Moon interference
- Finds optimal nights for deep-sky observation
- Combines Sun and Moon altitude periods
- Demonstrates building real-world planning tools with minimal code

**Run with:**
```bash
cargo run --example night_quality_scoring
```

## Key Concepts

### Reference Centers

- **Barycentric**: Center of mass of the solar system
- **Heliocentric**: Center of the Sun
- **Geocentric**: Center of the Earth
- **Topocentric**: Observer's location on Earth's surface
- **Bodycentric**: Generic center for any orbiting body

### Reference Frames

- **Ecliptic**: The plane of Earth's orbit around the Sun
- **EquatorialMeanJ2000**: mean equator/equinox of J2000.0 (FK5/J2000 mean)
- **EquatorialMeanOfDate**: mean equator/equinox of date (precession applied)
- **EquatorialTrueOfDate**: true equator/equinox of date (precession + nutation)
- **Horizontal**: Local horizon plane (altitude/azimuth)
- **ICRS**: International Celestial Reference System (fixed to distant quasars)
- **ECEF**: Earth-Centered Earth-Fixed (rotates with Earth)

### Type Safety

All coordinate types are parameterized by:
- **Center** (`C`): Where the origin is located
- **Frame** (`F`): How the axes are oriented
- **Unit** (`U`): What kind of vector (Position with length, Direction unitless, Velocity)

This ensures compile-time safety - you can't accidentally mix incompatible coordinate systems!

```rust
// Type-safe - won't compile if centers/frames don't match
let pos1: Position<Geocentric, EquatorialMeanJ2000, Km> = ...;
let pos2: Position<Geocentric, EquatorialMeanJ2000, Km> = ...;
let distance = pos1.distance_to(&pos2);  // ✓ OK

let pos3: Position<Heliocentric, Ecliptic, Au> = ...;
// pos1.distance_to(&pos3);  // ✗ Compile error - different types!
```

## Running All Examples

To run all examples:

```bash
cargo run --example basic_coordinates
cargo run --example coordinate_transformations
cargo run --example bodycentric_coordinates
cargo run --example observer_coordinates
cargo run --example solar_system_example
cargo run --example serde_serialization --features serde
```

## Tips

1. **Units**: The `qtty` crate provides unit-safe quantities. Use `AU`, `KM`, `M`, `DEG`, etc.
2. **Time**: Most transformations require a `JulianDate` for time-dependent calculations
3. **Conversions**: Use `.to_frame()`, `.to_center(jd)`, or `.transform(jd)` for conversions
4. **Precision**: Astronomical calculations use f64 precision throughout

## Further Reading

- See the main library documentation: `cargo doc --open`
- Check the test files in `tests/` for more usage patterns
- Read the module documentation in `src/coordinates/`
