// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic body-surface geometric Alt/Az (Moon, Mars, Earth).
//!
//! Demonstrates the shared [`SurfaceSite`] / ENU pipeline for observers on
//! the Moon, Mars, and Earth. Results are **geometric** (no refraction,
//! aberration, or light-time).
//!
//! Conventions:
//! - Azimuth: north = 0°, east = 90° (north through east)
//! - Altitude: geometric horizon = 0°, zenith = +90°
//! - Moon frame: IAU 2015 analytical `MoonPrincipalAxes` (not JPL DE PA / ME)
//! - Moon surface: IAU mean sphere R = 1737.4 km
//! - Mars surface: IAU ellipsoid a = 3396.19 km, planetographic latitude
//! - Earth surface: WGS84 / ECEF
//!
//! Run with: `cargo run --example 21_body_surface_horizontal`
#![allow(clippy::print_stdout)]

use siderust::coordinates::cartesian::{Direction, Position};
use siderust::coordinates::centers::{Geocentric, Selenocentric};
use siderust::coordinates::frames::{EclipticMeanJ2000, ICRS};
use siderust::coordinates::surface::{
    EarthSurface, MarsSurface, MoonSurface, SurfaceSite, EARTH_WGS84_SURFACE, MARS_IAU_ELLIPSOID,
    MOON_IAU_SPHERE,
};
use siderust::coordinates::transform::{
    direction_to_horizontal, horizontal_to_direction, position_to_horizontal, AstroContext,
    DirectionAstroExt, SurfaceHorizontalExt,
};
use siderust::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::qtty::*;
use siderust::time::JulianDate;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let jd = JulianDate::new(2_460_000.5); // TT
    let ctx: AstroContext = AstroContext::default();

    println!("Geometric body-surface Alt/Az");
    println!("Epoch (TT JD): {jd:.6}");
    println!("Azimuth: N=0° → E=90°; Altitude: horizon=0°, zenith=+90°");
    println!("Corrections: none (geometric only)\n");

    // ── Moon: observe Earth from a lunar site ─────────────────────────────
    let moon_site = SurfaceSite::<MoonSurface>::new(20.0 * DEG, (-3.5) * DEG, 0.0 * M)?;
    println!("Moon observer");
    println!(
        "  site: lon={:.3}, lat={:.3}, h={}, surface=IAU sphere R={} km",
        moon_site.lon(),
        moon_site.lat(),
        moon_site.height(),
        MOON_IAU_SPHERE
            .equatorial_radius()
            .to::<Kilometer>()
            .value()
    );
    println!(
        "  frame: MoonPrincipalAxes (IAU 2015 α0,δ0,W); lat={}",
        moon_site.latitude_convention().as_str()
    );

    // Earth as seen from the Moon centre = −Moon(geocentric) in the ecliptic frame.
    let moon_geo = Vsop87Ephemeris::moon_geocentric(jd);
    let earth_sel = Position::<Selenocentric, EclipticMeanJ2000, Kilometer>::new(
        -moon_geo.x(),
        -moon_geo.y(),
        -moon_geo.z(),
    );
    let lunar_obs = position_to_horizontal(&earth_sel, &moon_site, &jd, &ctx)?;
    println!(
        "  target: Earth (selenocentric, finite)  alt={:.4}  az={:.4}  range={}",
        lunar_obs.alt(),
        lunar_obs.az(),
        lunar_obs.range().to::<Kilometer>()
    );

    // Direction-only Sun (no parallax invented).
    let sun_bary = Vsop87Ephemeris::sun_barycentric(jd);
    let r = (sun_bary.x().value().powi(2)
        + sun_bary.y().value().powi(2)
        + sun_bary.z().value().powi(2))
    .sqrt();
    let sun_ecl = Direction::<EclipticMeanJ2000>::new_unchecked(
        sun_bary.x().value() / r,
        sun_bary.y().value() / r,
        sun_bary.z().value() / r,
    );
    let sun_icrs: Direction<ICRS> = sun_ecl.to_frame(&jd);
    let sun_hz = moon_site.observe_direction(&sun_icrs, &jd, &ctx)?;
    let sun_sph = sun_hz.to_spherical();
    println!(
        "  target: Sun (direction-only)          alt={:.4}  az={:.4}",
        sun_sph.alt(),
        sun_sph.az()
    );

    let back: Direction<ICRS> = horizontal_to_direction(&sun_hz, &moon_site, &jd, &ctx)?;
    println!(
        "  inverse check (Sun): angular error = {:.3e} rad\n",
        sun_icrs.angle_to(&back).value()
    );

    // ── Mars: observe a fixed ICRS direction ──────────────────────────────
    let mars_site = SurfaceSite::<MarsSurface>::new(137.4 * DEG, (-4.6) * DEG, 0.0 * M)?;
    println!("Mars observer (Gale-crater-like)");
    println!(
        "  site: lon={:.3}, lat={:.3} ({}), surface=IAU ellipsoid a={} km",
        mars_site.lon(),
        mars_site.lat(),
        mars_site.latitude_convention().as_str(),
        MARS_IAU_ELLIPSOID
            .equatorial_radius()
            .to::<Kilometer>()
            .value()
    );
    println!("  frame: MarsFixed (IAU 2015 α0,δ0,W)");

    let n = (0.1_f64.powi(2) + 0.7_f64.powi(2) + 0.7_f64.powi(2)).sqrt();
    let star = Direction::<ICRS>::new_unchecked(0.1 / n, 0.7 / n, 0.7 / n);
    let mars_hz = direction_to_horizontal(&star, &mars_site, &jd, &ctx)?;
    let mars_sph = mars_hz.to_spherical();
    println!(
        "  target: fixed ICRS direction          alt={:.4}  az={:.4}",
        mars_sph.alt(),
        mars_sph.az()
    );
    let back_m: Direction<ICRS> = horizontal_to_direction(&mars_hz, &mars_site, &jd, &ctx)?;
    println!(
        "  inverse check: angular error = {:.3e} rad\n",
        star.angle_to(&back_m).value()
    );

    // ── Earth: same generic API with WGS84 ────────────────────────────────
    let earth_site = SurfaceSite::<EarthSurface>::new((-17.89) * DEG, 28.75 * DEG, 2326.0 * M)?;
    println!("Earth observer (La Palma-like)");
    println!(
        "  site: lon={:.3}, lat={:.3}, h={}, surface=WGS84 a={:.0} m",
        earth_site.lon(),
        earth_site.lat(),
        earth_site.height(),
        EARTH_WGS84_SURFACE.equatorial_radius().value()
    );
    println!("  frame: ECEF (ITRF orientation chain via ICRS hub)");

    let moon_target = Position::<Geocentric, EclipticMeanJ2000, Kilometer>::new(
        moon_geo.x(),
        moon_geo.y(),
        moon_geo.z(),
    );
    let earth_obs = earth_site.observe_position(&moon_target, &jd, &ctx)?;
    println!(
        "  target: Moon (geocentric, finite)     alt={:.4}  az={:.4}  range={}",
        earth_obs.alt(),
        earth_obs.az(),
        earth_obs.range()
    );

    Ok(())
}
