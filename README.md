# Siderust

[![Crates.io](https://img.shields.io/crates/v/siderust.svg)](https://crates.io/crates/siderust)
[![Docs.rs](https://docs.rs/siderust/badge.svg)](https://docs.rs/siderust)
[![CI](https://github.com/Siderust/siderust/actions/workflows/github-ci.yml/badge.svg)](https://github.com/Siderust/siderust/actions/workflows/github-ci.yml)
[![License: AGPL-3.0-only](https://img.shields.io/badge/license-AGPL--3.0--only-blue)](LICENSE)


> **Typed astronomy & satellite mechanics in safe Rust.**

Siderust provides ephemerides, coordinate transforms, time-scale handling, and orbit-analysis building blocks for scientific applications. The primary Rust implementation avoids `unsafe`; allocation behavior depends on the enabled subsystems and is documented at the module level. External reference fixtures are being expanded, so scientific claims should be read against the validation tests that ship with the crate.

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

| Feature        | Default | What it enables |
|----------------|---------|-----------------|
| `serde`        | ✔       | `Serialize` / `Deserialize` on public types (default) |
| *(base)*       |         | VSOP87 + ELP2000-82B analytical ephemerides, full coordinate/altitude API |
| `atmosphere`   |         | Atmospheric tables and radiative transfer helpers |
| `photometry`   |         | Photometric passbands and throughput unit (Johnson–Cousins UBVRI) |
| `spice`        |         | High-level SPICE kernel context (`SpiceContext`, `KernelSet`) |
| `pod`          |         | Precise Orbit Determination toolkit (WLS, EKF, force models, I/O) |
| `pod-parquet`  |         | Parquet residuals writer (implies `pod`) |
| `pod-doris`    |         | DORIS RINEX observation parser (implies `pod`) |
| `runtime-data` |         | Runtime dataset-loading helpers via `siderust-archive` |

> **Note:** `no_std` and `f128` quad‑precision are **not** supported today.
> The crate depends on `std`‑only libraries such as `chrono`.
> Sub‑crates `qtty` and `qtty-core` do offer `no_std` support independently.
> C ABI bindings live in the separate `siderust-ffi` crate rather than a `siderust` feature flag.

---

## Features

| Category                | What you get                                                                                                                                                                  |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Coordinate Systems**  | `Position` and spherical `Direction` types parameterised by `ReferenceCenter`, `ReferenceFrame`, and `Unit`. Compile‑time guarantees prevent mixing frames by accident.       |
| **Target Tracking**     | `Target<T>` couples any coordinate with an observation epoch and optional `ProperMotion`, enabling extrapolation & filtering.                                                  |
| **Physical Units**      | Strongly typed `Mass`, `Length`, `Angle`, `Velocity`, `Duration` & more via the [`qtty`](https://crates.io/crates/qtty) crate, dimensional correctness at compile time.       |
| **Celestial Mechanics** | VSOP87 & ELP2000 theories, Pluto (Meeus/Williams), light‑time & aberration, nutation & precession, apparent Sun & Moon, culmination searches, SGP4/TLE propagation.         |
| **Ephemeris Backends**  | Pluggable `Ephemeris` / `DynEphemeris` traits: `Vsop87Ephemeris` (always available) and [`RuntimeEphemeris`](https://docs.rs/siderust/latest/siderust/ephemeris/struct.RuntimeEphemeris.html) for JPL DE4xx BSP files at runtime. |
| **Altitude API**        | Unified `AltitudePeriodsProvider` trait for Sun, Moon, stars, and arbitrary ICRS directions, find crossings, culminations, altitude ranges, and above/below‑threshold periods.|
| **Catalogs & Bodies**   | Built‑in Sun→Neptune, asteroids (Ceres, Bennu, Apophis), comets (Halley, Encke, Hale-Bopp), a starter star catalog, + helpers for custom datasets.                           |
| **Observatories**       | Predefined sites (Roque de los Muchachos, El Paranal, Mauna Kea, La Silla) with `ObserverSite` for topocentric transforms.                                                  |

Coordinate algebra and reusable conic geometry are provided by [`affn`](https://crates.io/crates/affn); Kepler-equation solving and domain-neutral conic propagation live in [`keplerian`](https://crates.io/crates/keplerian); `siderust` adds astronomy-specific time, frame transforms, ephemeris backends, and body/observer orchestration on top.

### API Design Pillars

Siderust is built on two cross-cutting principles documented in
[`doc/conventions.md`](doc/conventions.md):

1. **Typed quantities everywhere.** Every scalar that has physical meaning —
   pressures, scale heights, optical depths, airmasses, albedos, illumination
   fractions, CIP coordinates — is a `qtty` newtype. Passing a raw `f64` where
   a `Hectopascals` or `Kilometers` is expected is a compile-time error.

2. **Phantom-typed model selection.** Algorithm variants (e.g. nutation models)
   are selected at the call site via zero-sized phantom type parameters such as
   `to_frame_as::<EquatorialFrame, Iau2006A>(jd)`. There are no runtime enums
   to match on. Dispatch is fully monomorphised.

### Astrometry Compliance Note

- Stellar aberration uses the full special-relativistic (Lorentz) formula per IERS Conventions (2020, §7.2); annual uses VSOP87E barycentric Earth velocity and topocentric adds a diurnal `ω×r` term.
- The Earth-orientation chain now exposes public frame transforms for `GCRS ↔ CIRS ↔ TIRS ↔ ITRF/ECEF`, plus the usual inertial and operational frames (`ICRS`, `ICRF`, `EME2000`, `TEME`, `Galactic`, planetary body-fixed).
- Local orbital frames (`RTN`, `LVLH`, `VNC`) and covariance transport in those frames are first-class via [`astro::dynamics::frames::LocalOrbitalFrame`] and [`astro::dynamics::covariance::StateCovariance`], built from a typed [`OrbitState`].

---

## Scientific Archive (submodule)

Large scientific datasets — VSOP87, IAU 2000A nutation, ELP2000-82B,
Sun-Earth Lagrange Chebyshev kernels, SPICE time/frame/constants kernels,
and dataset generators/validators — live in the separate
[`Siderust/archive`](https://github.com/Siderust/archive) repository, attached
here as the `archive/` git submodule.

After cloning Siderust, initialise the submodule:

```sh
git submodule update --init --recursive
```

All archive metadata uses **TOML** (`MANIFEST.toml`, per-family
`manifest.toml`). Binary payloads use their authoritative formats (SPICE
`.bsp`, Siderust Chebyshev Kernel `.sck`, raw `.dat`). The archive layout,
manifest schema, and regeneration recipes are documented in
[`archive/README.md`](archive/README.md) and
[`archive/schema/archive-manifest-v1.md`](archive/schema/archive-manifest-v1.md).

The default `siderust` build does **not** require a separate archive checkout.
Large scientific datasets (VSOP87, ELP2000, nutation, gravity, atmosphere, Pluto)
are embedded via the [`siderust-archive`](https://crates.io/crates/siderust-archive)
crate which is a regular Cargo dependency.  JPL DE4xx kernels are resolved at
runtime from the local filesystem or downloaded on demand via `runtime-data`.

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
siderust = "0.9"
```

VSOP87/ELP2000 coefficients, nutation tables, and EOP data are provided by
[`siderust-archive`](https://crates.io/crates/siderust-archive) (scientific
datasets, manifests, checksums, provenance) and
[`tempoch`](https://crates.io/crates/tempoch) (UTC/TAI/TT/UT1/TDB time scales,
ΔT, and EOP freshness). Optional JPL kernels are downloaded on demand when the
corresponding feature is enabled; see `doc/datasets.md`.

### Ephemeris Backends

Siderust always includes [`Vsop87Ephemeris`](https://docs.rs/siderust/latest/siderust/ephemeris/struct.Vsop87Ephemeris.html) (VSOP87 + ELP2000-82B). Compile-time [`DefaultEphemeris`](https://docs.rs/siderust/latest/siderust/coordinates/transform/context/type.DefaultEphemeris.html) is an alias to VSOP87.

JPL DE440/DE441 (and other NAIF BSP kernels) are provided at **runtime** via [`RuntimeEphemeris`](https://docs.rs/siderust/latest/siderust/ephemeris/struct.RuntimeEphemeris.html):

```toml
[dependencies]
siderust = { version = "0.9", features = ["runtime-data"] }
```

```rust
use siderust::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::time::JulianDate;

let eph = RuntimeEphemeris::from_bsp("/path/to/de440.bsp")?;
let earth = eph.earth_heliocentric(JulianDate::J2000);
```

With `runtime-data`, [`siderust_archive::jpl::DatasetManager`](https://docs.rs/siderust-archive/latest/siderust_archive/jpl/struct.DatasetManager.html) can download kernels on first use (see example `12_runtime_ephemeris`).

VSOP87-only (explicit):

```toml
[dependencies]
siderust = { version = "0.9", default-features = false }
```

Combine with other features (for example `serde`):

```toml
[dependencies]
siderust = { version = "0.9", features = ["runtime-data", "serde"] }
```

### JPL datasets (runtime, not compile-time features)

DE kernels are **not** embedded at compile time. Typical files:

- `de440.bsp` (~120 MB)
- `de441_part-2.bsp` (~1.65 GB)

Prefetch into a persistent cache (optional):

```bash
export SIDERUST_DATASETS_DIR="$HOME/.cache/siderust"
# from an archive checkout:
../archive/scripts/prefetch_datasets.sh --de440
```

`cargo test --all-features` does **not** download JPL kernels by itself. Use `cargo test --features runtime-data` when exercising download paths, or point `SIDERUST_BSP_PATH` at a local BSP for integration tests (see `tests/test_jpl_real_backend.rs`).

CI sets `SIDERUST_JPL_STUB=all` for deterministic `--all-features` jobs. For local parity:

```bash
SIDERUST_JPL_STUB=all cargo test --all-features
```

Real-kernel validation (manual / optional workflow):

```bash
SIDERUST_BSP_PATH=/path/to/de440.bsp cargo test --test test_jpl_real_backend
SIDERUST_BSP_PATH=/path/to/de440.bsp cargo bench --bench de441
```

Coverage:
- Fallible JPL APIs such as `try_position`, `try_velocity`, and `try_position_velocity` return an error outside the Chebyshev segment coverage.
- The legacy infallible JPL APIs are retained for compatibility and panic with an explicit out-of-range message instead of silently extrapolating.
- Earth-orientation lookups are keyed by UTC/MJD. The default IERS provider uses `tempoch`'s bundled EOP data and reports missing coverage rather than fabricating zero EOP values. Use `NullEop` only when a documented zero-EOP approximation is intended.

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

The public frame set includes:
- inertial / catalogue frames: `ICRS`, `ICRF`, `EquatorialMeanJ2000`, `EME2000`, `EquatorialMeanOfDate`, `EquatorialTrueOfDate`, `FK4B1950`, `Galactic`
- Earth-rotation chain frames: `GCRS`, `CIRS`, `TIRS`, `ITRF`, `ECEF`
- operational / mission frames: `TEME`, `Horizontal`, and the planetary body-fixed frames

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
- Run one: `cargo run --example 01_basic_coordinates`

Feature-gated examples:

```bash
# Runtime JPL ephemeris (BSP path argument or runtime-data download)
cargo run --example 12_runtime_ephemeris -- /path/to/de440.bsp
cargo run --features runtime-data --example 12_runtime_ephemeris

# Serde
cargo run --example 11_serde_serialization --features serde
```

## Crate Layout

```
├─ astro/         # Aberration, nutation, precession, sidereal time
├─ bodies/        # Planet, Star, Satellite, Asteroid, Comet + built-in catalogs
├─ calculus/
│   ├─ altitude/     # Unified altitude API (AltitudePeriodsProvider trait)
│   ├─ ephemeris/    # Ephemeris trait + VSOP87 + runtime JPL backends
│   ├─ jpl/          # Shared JPL DE4xx infrastructure (Chebyshev evaluation)
│   ├─ math_core/    # Root-finding (Brent/bisection), extrema, interval assembly
│   ├─ solar/        # Sun altitude, night/day/twilight periods
│   ├─ lunar/        # Moon altitude with topocentric parallax
│   ├─ stellar/      # Analytical star altitude engine
│   ├─ vsop87/       # VSOP87 planetary theory
│   ├─ elp2000/      # ELP2000-82B lunar theory
│   ├─ conic_equations.rs # Siderust wrappers over keplerian solvers
│   └─ pluto         # Meeus/Williams Pluto ephemeris
├─ coordinates/   # Cartesian/Spherical types, frames, centers, transforms
├─ mission/       # Mission-analysis building blocks
│   ├─ geometry/     # AzElRange, Fov, TerrainMask, eclipse, orbit-relative geometry
│   ├─ context.rs    # MissionContext — runtime aggregation of instruments and sites
│   └─ site.rs       # Location — ground-station / observing-site metadata
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
