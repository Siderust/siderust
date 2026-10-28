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
3. `23_lunar_landing_site_navigator` — solar and terrestrial Alt/Az over a lunar
   site and a time interval, with an optional CSV export for comparing sites.
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

The [Lunar Landing Site Navigator challenge](https://www.spaceappschallenge.org/2026/challenges/lunar-landing-site-navigator/)
asks participants to compare Sun and Earth positions relative to the horizon at
candidate lunar landing sites and dates. Example `23_lunar_landing_site_navigator`
is a no-download Rust CLI starting point: enter a lunar site and a TT Julian
start date, sample a time span, and optionally export CSV for plots or a UI.

```bash
# Defaults: latitude -89 degrees, east longitude 0 degrees; 48 h, 60 min steps
cargo run --example 23_lunar_landing_site_navigator

# LAT_DEG LON_EAST_DEG START_JD_TT [HOURS STEP_MINUTES]
cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5 48 60

# Write rows of JD(TT), Sun/Earth azimuth, altitude, and above-horizon flags
cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5 48 60 --csv > lunar.csv

# Repeat with a different site to compare
cargo run --example 23_lunar_landing_site_navigator -- -88.5 30 2461358.5 48 60 --csv > second_site.csv
```

All coordinates are lunar planetocentric latitude / **east-positive longitude**;
zero height is on the IAU 1737.4 km reference sphere. Input Julian dates are in
**TT, not UTC**. Samples include both endpoints (possibly a shorter last step).
The CSV flags mean only that the Sun's or Earth's *centre* is above the
reference-sphere geometric horizon. This does **not** model terrain/LOLA
obstructions, the solar disc, antenna masks, light-time, aberration, power
availability, or radio-link quality; even uninterrupted positive samples do not
establish real illumination or a direct-to-Earth communications window.

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
