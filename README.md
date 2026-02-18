# Siderust

[![Crates.io](https://img.shields.io/crates/v/siderust.svg)](https://crates.io/crates/siderust)
[![Docs.rs](https://docs.rs/siderust/badge.svg)](https://docs.rs/siderust)

> **Precision astronomy & satellite mechanics in safe, fast Rust.**

Siderust aims to be a reference ephemeris and orbit‑analysis library for research‐grade pipelines and ground‑segment tooling. Every algorithm ships with validation tests against authoritative data (JPL Horizons, IMCCE, SOFA). No unsafe blocks, no hidden allocations.

---

## Table of Contents

1. [Supported Feature Flags](#supported-feature-flags)
2. [Features](#features)
3. [Installation](#installation)
4. [Coordinate Systems](#coordinate-systems)
5. [Units & Physical Quantities](#units--physical-quantities)
6. [Quick Start](#quick-start)
7. [Examples](#examples)
8. [Crate Layout](#crate-layout)
9. [Roadmap](#roadmap)
10. [Contributing](#contributing)
11. [License](#license)
12. [Acknowledgments](#acknowledgments)

---

## Supported Feature Flags

| Feature  | Default | What it enables |
|----------|---------|-----------------|
| *(none)* | ✔       | VSOP87 + ELP2000-82B analytical ephemerides, full coordinate/altitude API |
| `de440`  |         | JPL DE440 Chebyshev ephemeris backend (1550–2650 CE) |
| `de441`  |         | JPL DE441 Chebyshev ephemeris backend (extended coverage) |
| `serde`  |         | `Serialize` / `Deserialize` on public types |

> **Note:** `no_std` and `f128` quad‑precision are **not** supported today.
> The crate depends on `std`‑only libraries (`chrono`, `nalgebra`).
> Sub‑crates `qtty` and `qtty-core` do offer `no_std` support independently.

---

## Features

| Category                | What you get                                                                                                                                                                  |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Coordinate Systems**  | `Position` and spherical `Direction` types parameterised by `ReferenceCenter`, `ReferenceFrame`, and `Unit`. Compile‑time guarantees prevent mixing frames by accident.       |
| **Target Tracking**     | `Target<T>` couples any coordinate with an observation epoch and optional `ProperMotion`, enabling extrapolation & filtering.                                                  |
| **Physical Units**      | Strongly typed `Mass`, `Length`, `Angle`, `Velocity`, `Duration` & more via the [`qtty`](https://crates.io/crates/qtty) crate — dimensional correctness at compile time.       |
| **Celestial Mechanics** | Kepler solvers, VSOP87 & ELP2000 theories, Pluto (Meeus/Williams), light‑time & aberration, nutation & precession, apparent Sun & Moon, culmination searches.                |
| **Ephemeris Backends**  | Pluggable `Ephemeris` trait with three backends — `Vsop87Ephemeris` (always available), `De440Ephemeris`, and `De441Ephemeris` (feature-gated JPL DE4xx).                     |
| **Altitude API**        | Unified `AltitudePeriodsProvider` trait for Sun, Moon, stars, and arbitrary ICRS directions — find crossings, culminations, altitude ranges, and above/below‑threshold periods.|
| **Catalogs & Bodies**   | Built‑in Sun→Neptune, asteroids (Ceres, Bennu, Apophis), comets (Halley, Encke, Hale-Bopp), a starter star catalog, + helpers for custom datasets.                           |
| **Observatories**       | Predefined sites (Roque de los Muchachos, El Paranal, Mauna Kea, La Silla) with `ObserverSite` for topocentric transforms.                                                  |

### Astrometry Compliance Note

- Stellar aberration uses the full special-relativistic (Lorentz) formula per IERS Conventions (2020, §7.2); annual uses VSOP87E barycentric Earth velocity and topocentric adds a diurnal `ω×r` term (GMST-based Earth rotation).
- This is not yet a full IAU 2000/2006 "apparent place" pipeline (missing CIO/CIP, polar motion, and gravitational light deflection).

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
siderust = "0.5"
```

Build-time datasets (VSOP87/ELP2000/IERS and optional JPL kernels) are downloaded on demand; see `doc/datasets.md`.

### Ephemeris Backends (Enable / Disable / Combine)

Siderust always includes `Vsop87Ephemeris` (VSOP87 + ELP2000-82B).
Optional features add JPL backends:
- `de440` → `De440Ephemeris` (1550–2650 CE)
- `de441` → `De441Ephemeris` (extended coverage from NAIF `de441_part-2.bsp`)

`DefaultEphemeris` selects the best available:
- `De441Ephemeris` when `de441` is enabled
- otherwise `De440Ephemeris` when `de440` is enabled
- otherwise `Vsop87Ephemeris`

1. VSOP87-only (explicit)

```toml
[dependencies]
siderust = { version = "0.5", default-features = false }
```

2. Enable DE440

```toml
[dependencies]
siderust = { version = "0.5", features = ["de440"] }
```

3. Enable DE441

```toml
[dependencies]
siderust = { version = "0.5", features = ["de441"] }
```

4. Combine backends in one binary

```rust
use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::time::JulianDate;

let jd = JulianDate::J2000;

// Analytical series (always available)
let earth_vsop = Vsop87Ephemeris::earth_heliocentric(jd);
```

```rust
// With `de441` feature enabled:
use siderust::calculus::ephemeris::De441Ephemeris;

let earth_jpl = De441Ephemeris::earth_heliocentric(jd);
```

You can combine ephemeris features with others (for example `serde`):

```toml
[dependencies]
siderust = { version = "0.5", features = ["de441", "serde"] }
```

### JPL Build Modes: Real vs Stubbed

When JPL features are enabled, build scripts may download:
- `de440.bsp` (~120 MB)
- `de441_part-2.bsp` (~1.65 GB)

Default checked-in behavior is **real JPL builds** (no global stubbing).

Real JPL mode (recommended for representative DE440/DE441 behavior):

```bash
unset SIDERUST_JPL_STUB
cargo test --features de440
cargo test --features de441
```

Offline/stub mode (explicit opt-in for fast local loops):

```bash
SIDERUST_JPL_STUB=all cargo check --all-features
```

Supported stub values:
- `de441`: stubs DE441 only (`De441Ephemeris` is mocked to `Vsop87Ephemeris`).
- `de440`: stubs DE440 only.
- `de440,de441` or `all` (also `1`, `true`, `yes`, `on`): stubs both.

Optional local override file (keep untracked):

```toml
# .cargo/config.local.toml
[env]
SIDERUST_JPL_STUB = "all"
```

Use it explicitly:

```bash
cargo --config .cargo/config.local.toml check --all-features
```

Caveat:
- Stubbed DE datasets compile successfully, but low-level DE calls are unavailable at runtime.
- `de441` has a high-level mock backend (`De441Ephemeris -> Vsop87Ephemeris`).
- `de440` is compile-only when stubbed; direct runtime calls to `De440Ephemeris` will panic.

---

## Coordinate Systems

Siderust encodes both the **origin** and the **orientation** of every coordinate at the type level:

```rust
use siderust::coordinates::{cartesian, centers::*, frames::*};
use qtty::Au;

// Position of Mars in the Heliocentric Ecliptic frame
let mars_helio = cartesian::Position::<Heliocentric, EclipticMeanJ2000, Au>::new(1.5, 0.0, 0.0);
```

Impossible states (e.g. adding heliocentric and geocentric positions) simply do not compile.

### Compile-Time vs Runtime Safety

For most centers (`Barycentric`, `Heliocentric`, `Geocentric`) all invariants are
enforced at **compile time** with zero runtime cost (`Params = ()`).

**Parameterized centers** (`Topocentric`, `Bodycentric`) carry runtime data
(e.g., `ObserverSite`). The *center type* is still checked at compile time, but
*parameter equality* (e.g., "are these two positions at the same site?") is
checked at **runtime**:

| API | Behaviour on mismatch |
|-----|-----------------------|
| `pos_a - pos_b` / `distance_to` | `assert!` (panics in all builds) |
| `checked_sub` / `try_distance_to` | Returns `Err(CenterParamsMismatchError)` |
| `ObserverSite::try_new` | Validates lat/lon ranges, returns `Result` |

---

## Units & Physical Quantities

Siderust uses the [`qtty`](https://crates.io/crates/qtty) crate for dimensionally
typed quantities. The compiler prevents mixing incompatible units:

```rust
use qtty::*;

let distance = AstronomicalUnits::new(1.523); // Mars semi-major axis
let period   = Days::new(686.97);

// distance + period → compile error (length + time)
```

Common unit types: `AstronomicalUnit` (`Au`), `Kilometer` (`Km`), `Meter`,
`Degree`, `Radian`, `Day`, `Second`, `AuPerDay`, and many more.

---

## Quick Start

```rust
use siderust::{
    bodies::Mars,
    time::JulianDate,
};
use chrono::prelude::*;

// 1. Select an epoch (UTC now → JD)
let jd = JulianDate::from_utc(Utc::now());

// 2. Compute barycentric ecliptic coordinates via VSOP87
let mars = Mars::vsop87e(jd);

// 3. Print Mars's barycentric ecliptic position
println!("{}", mars.position);
```

---

## Examples

The `examples/` directory is a curated tour of the crate’s major building blocks
(coordinates, transforms, altitude periods, ephemeris backends, serialization).

- Browse: `examples/README.md`
- Run one: `cargo run --example basic_coordinates`

Feature-gated examples:

```bash
# JPL DE4xx (may download large BSP datasets unless you explicitly stub)
cargo run --example jpl_precise_ephemeris --features de440
cargo run --example jpl_precise_ephemeris --features de441

# Fast/offline loop: compile JPL features but skip runtime DE calls
SIDERUST_JPL_STUB=all cargo run --example jpl_precise_ephemeris --features de440,de441

# Serde
cargo run --example serde_serialization --features serde
```

## Crate Layout

```
├─ astro/         # Aberration, nutation, precession, sidereal time
├─ bodies/        # Planet, Star, Satellite, Asteroid, Comet + built-in catalogs
├─ calculus/
│   ├─ altitude/     # Unified altitude API (AltitudePeriodsProvider trait)
│   ├─ ephemeris/    # Ephemeris trait + VSOP87/DE440/DE441 backends
│   ├─ jpl/          # Shared JPL DE4xx infrastructure (Chebyshev evaluation)
│   ├─ math_core/    # Root-finding (Brent/bisection), extrema, interval assembly
│   ├─ solar/        # Sun altitude, night/day/twilight periods
│   ├─ lunar/        # Moon altitude with topocentric parallax
│   ├─ stellar/      # Analytical star altitude engine
│   ├─ vsop87/       # VSOP87 planetary theory
│   ├─ elp2000/      # ELP2000-82B lunar theory
│   ├─ kepler_equations/  # Kepler equation solvers
│   └─ pluto         # Meeus/Williams Pluto ephemeris
├─ coordinates/   # Cartesian/Spherical types, frames, centers, transforms
├─ observatories/ # Predefined observatory locations (Roque, Paranal, Mauna Kea, La Silla)
├─ targets/       # Target<T> with time & ProperMotion
└─ time           # Re-export of tempoch: JulianDate, MJD, Period<S>, time scales
```

---

## Roadmap

* [x] Custom dynamic reference centers (topocentric, bodycentric)
* [x] DE440/DE441 JPL ephemerides
* [x] Unified altitude API (`AltitudePeriodsProvider` trait)
* [x] Serde serialization support
* [ ] Gaia DR3 star ingestion & cone search
* [ ] Relativistic light‑time & gravitational deflection
* [ ] Batch orbit determination helpers (LSQ & EKF)
* [ ] GPU acceleration via `wgpu` (experiment)

---

## Contributing

Contributions of algorithms, bug fixes or docs are welcome! Please:

1. Fork & clone (`git clone`)
2. Create a feature branch
3. Run **all** tests & clippy (`cargo test && cargo clippy -- -D warnings`)
4. Open a PR with a clear description

By participating you agree to follow the [Rust Code of Conduct](https://www.rust-lang.org/policies/code-of-conduct).

---

## License

Copyright (C) 2026 Vallés Puig, Ramon

This project is licensed under the **GNU Affero General Public License v3.0 (AGPL-3.0)**.
The AGPL-3.0 ensures that:

- Any modifications and redistributions of the code (including as a network service) remain free and open.
- End users have access to the full source code, including any improvements or extensions made.

> **Note for commercial or proprietary use:**
> If you wish to incorporate this code into a closed-source or otherwise differently licensed project, a **dual-licensing** arrangement can be negotiated. Please contact the authors to discuss terms and conditions for a commercial or proprietary license that suits your needs.

---

## Acknowledgments

Big thanks to **Màrius Montón** ([@mariusmm](https://github.com/mariusmm)) for inviting me to his three-week Rust intro course at the Universitat Autònoma de Barcelona (UAB) in 2024 that nudge set this project in motion.
