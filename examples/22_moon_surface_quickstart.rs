// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Minimal lunar surface Alt/Az example, with no external ephemeris files.
//!
//! Run: cargo run --example 22_moon_surface_quickstart
//! Azimuth: north through east; altitude: geometric horizon through zenith.
//! Longitude positive east; planetocentric latitude on a 1737.4-km sphere.
//! Local up is radial. Orientation: NAIF IAU_MOON analytical Mean Earth,
//! exposed by the legacy MoonPrincipalAxes frame marker (not DE-specific PA).
//! J2000 is interpreted as a TT date, approximating the model's TDB input.
//! Alt/Az is geometric, not apparent; no terrain occultation, refraction,
//! aberration, or light-time.
#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::Moon;
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::frames::ICRS;
use siderust::coordinates::surface::SurfaceSite;
use siderust::coordinates::transform::{AstroContext, SurfaceHorizontalExt};
use siderust::qtty::{DEG, M};
use siderust::J2000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Longitude east, latitude planetocentric; height above a 1737.4-km sphere.
    let site = SurfaceSite::<Moon>::new(20.0 * DEG, (-3.5) * DEG, 0.0 * M)?;

    // A fixed, unit-length ICRS direction toward the north celestial pole.
    // A direction has no distance, so parallax is not fabricated.
    let star = Direction::<ICRS>::new_unchecked(0.0, 0.0, 1.0);

    // The Moon's body-fixed orientation is evaluated at this TT epoch.
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
