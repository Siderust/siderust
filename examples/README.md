# Siderust Examples

Runnable examples organized by theme. Run any example with:

```bash
cargo run --example <name>
```

Feature-gated examples list the required feature in the command.

## Core Coordinates

- `01_basic_coordinates`: coordinate centers, frames, and typed units.
- `02_coordinate_transformations`: frame and center transforms.
- `03_all_frames_conversions`: supported frame rotations and round-trip checks.
- `04_all_center_conversions`: supported center shifts and round-trip checks.
- `13_coordinate_operations`: angular separation, distances, and displacement algebra.
- `14_nutation_models`: default and custom nutation-model transform paths.

```bash
cargo run --example 01_basic_coordinates
cargo run --example 02_coordinate_transformations
cargo run --example 03_all_frames_conversions
cargo run --example 04_all_center_conversions
cargo run --example 13_coordinate_operations
cargo run --example 14_nutation_models
```

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

```bash
cargo run --example 11_serde_serialization --features serde
cargo run --example 12_runtime_ephemeris --features de440
cargo run --example 12_runtime_ephemeris --features de441

# Fast/offline loop: compile JPL features while stubbing large DE datasets.
SIDERUST_JPL_STUB=all cargo run --example 12_runtime_ephemeris --features de440,de441
```
