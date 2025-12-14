# Siderust

## Development

### Running CI Checks Locally

```bash
# Run all CI checks (check, fmt, clippy, tests)
./ci-local.sh

# Run coverage checks (matches GitHub CI workflow)
./ci-local-coverage.sh
```

**Note**: Coverage requires the nightly toolchain. Install it with:
```bash
rustup toolchain install nightly
rustup component add llvm-tools-preview --toolchain nightly
cargo install cargo-llvm-cov
```

[![Crates.io](https://img.shields.io/crates/v/siderust.svg)](https://crates.io/crates/siderust)
[![Docs.rs](https://docs.rs/siderust/badge.svg)](https://docs.rs/siderust)

> **Precision astronomy & satellite mechanics in safe, fast Rust.**

Siderust aims to be the reference ephemeris and orbit‑analysis library for embedded flight‑software as well as research‐grade pipelines. Every algorithm ships with validation tests against authoritative data (JPL Horizons, IMCCE, SOFA). No unsafe blocks, no hidden allocations.

---

## Table of Contents

1. [Features](#features)
2. [Installation](#installation)
3. [Coordinate Systems](#coordinate-systems)
4. [Units & Physical Quantities](#units--physical-quantities)
5. [Quick Start](#quick-start)
6. [Accuracy & Benchmarks](#accuracy--benchmarks)
7. [Crate Layout](#crate-layout)
8. [Roadmap](#roadmap)
9. [Contributing](#contributing)
10. [License](#license)
11. [Acknowledgments](#acknowledgments)

---

## Features

| Category                | What you get                                                                                                                                                                                                         |
| ----------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Coordinate Systems**  | `Vector` & `SphericalCoord` parametrised by `ReferenceCenter` (Helio, Geo, Bary, …) and `ReferenceFrame` (ICRS, Ecliptic, Equatorial, Topocentric, etc.). Compile‑time guarantees ensure you never mix frames by accident. |
| **Target Tracking**     | `Target<T>` couples any coordinate with an observation epoch and optional `ProperMotion`, enabling extrapolation & filtering pipelines.                                                                              |
| **Physical Units**      | Strongly typed `Mass`, `Length`, `Angle`, `Velocity`, `Time` & more; operator overloading makes math look natural while the compiler guards dimensional correctness.                                             |
| **Celestial Mechanics** | Kepler solvers, VSOP87 & ELP2000 planetary/lunar theories, light‑time & aberration, nutation & precession matrices, apparent Sun & Moon, culmination searches.                                                       |
| **Catalogs & Bodies**   | Built‑in Sun→Neptune, major moons, a starter star catalog, + helper builders to load *Gaia*, *Hipparcos* or custom datasets.                                                                                         |

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
siderust = "0.1"
```

---

## Coordinate Systems

Siderust encodes both the **origin** and the **orientation** of every vector at the type level:

```rust
use siderust::coordinates::{Vector, centers::*, frames::*};

// Position of Mars in the Heliocentric Ecliptic frame
let mars_helio = Vector::<Heliocentric, Ecliptic>::new(x, y, z);

// Convert to Geocentric Ecliptic Cartesian coordinates
let mars_geo: Vector::<Geocentric, Ecliptic> = mars_helio.transform(jd);
```

Impossible states (e.g. adding heliocentric and geocentric vectors) simply do not compile.

---

## Units & Physical Quantities

```rust
use siderust::units::{AU, KM, DEG, DAY};
use siderust::units::*;

let distance = 1.523 * AU; // Mars semi‑major axis
let period   = 686.97 * DAY;
```

The compiler will refuse `distance + period` – dimensional analysis at compile time.

---

## Quick Start

```rust
use siderust::{
    bodies::Mars,
    astro::JulianDate,
};
use chrono::prelude::*;

// 1. Select an epoch (UTC now to JD)
let jd = JulianDate::from_utc(Utc::now());

// 2. Compute barycentric ecliptic coordinates via VSOP87
let mars = Mars::vsop87e(jd);

// 3. Print mars

println!("{}", mars.position);
```

---

## Accuracy & Benchmarks

All numeric kernels are cross‑checked against JPL Horizons (see **siderust-py** #TODO):<br/>
`|Δα|, |Δδ| < #TBD mas` for planets (1800–2200 CE). Full tables in `#TBD`.

| Routine              | Mean time (ns) | Note                         | HW            |
| -------------------- | -------------- | ---------------------------- | ------------- |
| VSOP87 planet        | **120**        | SIMD auto‑vectorised by LLVM | Ryzen 7 5800X |
| ELP2000 Moon         | **310**        |                              |               |
| Coordinate transform | **<50**        | center+frame change          |               |

--  This is just a placeholder to put the real data --

---


## Crate Layout

```
├─ astro/         # Astronomical properties (nutation, precession, …)
├─ bodies/        # Data structures for celestial bodies (planet, star, satellite, …)
├─ calculus/      # Numerical kernels (kepler, vsop87, …)
├─ coordinates/   # Coordinate types & transforms
├─ observatories/ # Ground stations & observer helpers
├─ targets/       # Target<T> & ProperMotion
└─ units/         # Dimensional quantities
```

---

## Roadmap

* [ ] Custom dynamic reference centers (topocentric)
* [ ] DE440/441 ephemerides (barycentric)
* [ ] Gaia DR3 star ingestion & cone search
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

By participating you agree to follow the [Rust Code of Conduct](https://www.rust-lang.org/policies/code-of-conduct).

---

## License

This project is licensed under the **GNU Affero General Public License v3.0 (AGPL-3.0)**.
The AGPL-3.0 ensures that:

- Any modifications and redistributions of the code (including as a network service) remain free and open.
- End users have access to the full source code, including any improvements or extensions made.

> **Note for commercial or proprietary use:**
> If you wish to incorporate this code into a closed-source or otherwise differently licensed project, a **dual-licensing** arrangement can be negotiated. Please contact the authors to discuss terms and conditions for a commercial or proprietary license that suits your needs.


---

## Acknowledgments

Big thanks to **Màrius Montón** ([@mariusmm](https://github.com/mariusmm)) for inviting me to his three-week Rust intro course at the Universitat Autònoma de Barcelona (UAB) in 2024 that nudge set this project in motion.