// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar surface Alt/Az quickstart (no downloads)
//!
//! This is the smallest useful body-surface example. It answers:
//! “From a site on the Moon, where is a fixed sky direction right now?”
//!
//! ## Concepts (in order)
//!
//! 1. **Observer site** — a longitude, latitude, and height on the Moon’s
//!    reference sphere (`SurfaceSite::<Moon>`). Longitude is east-positive;
//!    latitude is planetocentric on the IAU 1737.4 km sphere.
//! 2. **ICRS** — the inertial celestial reference system (quasar-fixed axes).
//!    A star catalogue direction is naturally an ICRS unit vector.
//! 3. **Direction** — a unit vector with **no distance**. Parallax is therefore
//!    not invented; only the orientation of the Moon matters.
//! 4. **Epoch** — the observation time as a TT Julian Date (`J2000` here).
//!    Lunar orientation is evaluated in TDB via tempoch’s TT→TDB conversion.
//! 5. **Altitude / azimuth** — geometric horizon angles at the site
//!    (alt = 0° on the horizon, +90° at zenith; az = 0° north, 90° east).
//! 6. **`AstroContext`** — carries ephemeris / EOP providers used by frame
//!    rotations. `Default` is enough for this geometric demo.
//!
//! Orientation model: NAIF `IAU_MOON` analytical Mean-Earth
//! (`MoonIauMeanEarth`), **not** a JPL DE `MOON_PA` kernel.
//! Results are **geometric**, not apparent (no refraction, aberration,
//! light-time, or terrain).
//!
//! ```text
//! cargo run --example 22_moon_surface_quickstart
//! ```
//!
//! Representative output (deterministic for this fixed site, direction, and
//! J2000 TT epoch; values change only if the lunar orientation model changes):
//! ```text
//! Lunar site: lon=20 °, lat=-3.5 °
//! Geometric Alt/Az at J2000 TT: alt=17.765928422856998 °, az=12.045094675187446 ° (N=0°, E=90°)
//! ```
#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::Moon;
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::frames::ICRS;
use siderust::coordinates::surface::SurfaceSite;
use siderust::coordinates::transform::{AstroContext, SurfaceHorizontalExt};
use siderust::qtty::{DEG, M};
use siderust::J2000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Observer on the lunar sphere (lon east, lat planetocentric, height 0).
    let site = SurfaceSite::<Moon>::new(20.0 * DEG, (-3.5) * DEG, 0.0 * M)?;

    // 2–3. Fixed ICRS direction toward the north celestial pole (unit vector).
    // `new_unchecked` is safe here: (0,0,1) is already a unit vector.
    let star = Direction::<ICRS>::new_unchecked(0.0, 0.0, 1.0);

    // 4–6. Evaluate geometric Alt/Az at J2000 TT with the default context.
    let ctx: AstroContext = AstroContext::default();
    let observed = site.observe_direction(&star, &J2000, &ctx)?;

    println!("Lunar site: lon={}, lat={}", site.lon(), site.lat());
    println!(
        "Geometric Alt/Az at J2000 TT: alt={}, az={} (N=0°, E=90°)",
        observed.alt(),
        observed.az()
    );
    Ok(())
}
