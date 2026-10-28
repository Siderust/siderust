# Siderust Examples

Runnable examples organized by theme. Run any example with:

```bash
cargo run --example <name>
```

Feature-gated examples list the required feature in the command.

## Body-surface Alt/Az (start here for Moon / Mars)

Beginner path for geometric altitude/azimuth on celestial bodies:

1. `22_moon_surface_quickstart` — **start here**. Minimal lunar site + ICRS
   direction → geometric Alt/Az. No downloads. Explains site, ICRS, directions
   (no parallax), epoch, Alt/Az, and `AstroContext`.
2. `21_body_surface_horizontal` — progressive Moon (finite targets with
   parallax) → Mars (direction-only) → Earth (WGS84) workflow.
3. `23_lunar_landing_site_navigator` — typed Sun/Earth lunar Alt/Az at six-hour
   intervals over two days, with a small command-line site/date interface.
4. `19_horizontal_coordinates` — Earth-only equatorial ↔ horizontal using
   `Topocentric` with **GMST** (positions) / **GAST** (directions). Scientifically
   distinct from the generic body-surface ENU pipeline.

Lunar orientation is the NAIF `IAU_MOON` analytical Mean-Earth model
(`MoonIauMeanEarth`), not a JPL DE `MOON_PA` / binary-PCK `MOON_ME` realization.
Longitude is east-positive; Mars SPICE planetographic longitudes are
west-positive and must be converted explicitly. Results are geometric (no
refraction, aberration, light-time, or terrain).

```bash
cargo run --example 22_moon_surface_quickstart
cargo run --example 21_body_surface_horizontal
cargo run --example 23_lunar_landing_site_navigator
cargo run --example 19_horizontal_coordinates
```

## Core Coordinates

- `01_basic_coordinates`: coordinate centers, frames, and typed units.
- `02_coordinate_transformations`: frame and center transforms.
- `03_all_frames_conversions`: supported frame rotations and round-trip checks.
- `04_all_center_conversions`: supported center shifts and round-trip checks.
- `13_coordinate_operations`: angular separation, distances, and displacement algebra.
- `14_nutation_models`: default and custom nutation-model transform paths.
- `19_horizontal_coordinates`: Earth-only equatorial ↔ horizontal (GMST/GAST).
- `21_body_surface_horizontal`: complete Moon/Mars/Earth geometric Alt/Az workflow.
- `22_moon_surface_quickstart`: minimal lunar geometric Alt/Az (no downloads).
- `23_lunar_landing_site_navigator`: lunar Sun/Earth geometric horizon sampling CLI.

```bash
cargo run --example 01_basic_coordinates
cargo run --example 02_coordinate_transformations
cargo run --example 03_all_frames_conversions
cargo run --example 04_all_center_conversions
cargo run --example 13_coordinate_operations
cargo run --example 14_nutation_models
cargo run --example 19_horizontal_coordinates
cargo run --example 22_moon_surface_quickstart
cargo run --example 21_body_surface_horizontal
cargo run --example 23_lunar_landing_site_navigator
```

## Lunar Landing Site Navigator (NASA Space Apps 2026)

`23_lunar_landing_site_navigator` demonstrates the [Lunar Landing Site Navigator
challenge](https://www.spaceappschallenge.org/2026/challenges/lunar-landing-site-navigator/)
using the same typed `SurfaceSite<Moon>`, `position::EclipticMeanJ2000`,
`TransformCenter` and `qtty` quantities as the rest of Siderust.

```bash
# South-polar default site, two days sampled every six hours
cargo run --example 23_lunar_landing_site_navigator

# LAT_DEG LON_EAST_DEG START_JD_TT
cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5
```

Latitude is planetocentric and longitude east-positive, both in degrees.
Time is a Julian Date in **TT** (not UTC). The output directly prints typed
altitude and azimuth quantities with degree units and indicates whether each
object's centre lies above the reference-sphere geometric horizon. To explore
other sites or dates, change the arguments; the sample interval is deliberately
fixed to keep the example focused on the Siderust API.

This is a geometric demonstration, **not** a terrain-aware or mission-grade
landing-site assessment. It does not model LOLA topography, partial solar-disc
visibility, aberration, light-time, antenna masks, power availability, or
communications link quality. The lunar orientation uses analytical `IAU_MOON`
rather than a JPL DE binary-PCK realization.

## Observing Workflows

- `05_target_tracking`: timestamped targets and proper-motion propagation.
- `06_night_events`: night and twilight event calculations.
- `07_moon_properties`: Moon position and illumination helpers.
- `09_star_observability`: star altitude windows at an observatory.

```bash
cargo run --example 05_target_tracking
cargo run --example 06_night_events
cargo run --example 07_moon_properties
cargo run --example 09_star_observability
```

## Solar System And Time

- `08_solar_system`: solar-system body positions and comparisons.
- `10_time_periods`: `tempoch` time scales and periods re-exported by `siderust`.

```bash
cargo run --example 08_solar_system
cargo run --example 10_time_periods
```

## Optional Features

- `11_serde_serialization`: JSON round-trips for time and coordinate values.
- `12_runtime_ephemeris`: runtime JPL ephemeris loading and DE feature checks.
- `18_lisa_pod`: load LISA CCSDS OEM samples and compute a light-time-corrected inter-spacecraft range.
- `20_mars_saturn_orientation`: geocentric Mars pole/bright-limb and Saturn ring orientation via the generic planetary orientation API.

```bash
cargo run --example 11_serde_serialization --features serde
cargo run --example 12_runtime_ephemeris -- /path/to/de440.bsp
cargo run --features runtime-data --example 12_runtime_ephemeris
cargo run --example 18_lisa_pod --features pod
cargo run --example 20_mars_saturn_orientation
```
