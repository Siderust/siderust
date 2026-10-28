// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Body-surface geometric Alt/Az (Moon → Mars → Earth)
//!
//! Progressive workflow after `22_moon_surface_quickstart`:
//!
//! 1. **Moon, finite target** — Earth and Sun as selenocentric positions
//!    (parallax included because range is known).
//! 2. **Mars, direction-only** — a fixed ICRS star direction (no parallax).
//! 3. **Earth, generic API** — same `SurfaceSite` pipeline with WGS84 / ECEF.
//!
//! Conventions:
//! - Azimuth: north = 0°, east = 90° (north through east)
//! - Altitude: geometric horizon = 0°, zenith = +90°
//! - Longitude: always east-positive in Siderust. Mars SPICE planetographic
//!   longitudes are west-positive — convert with
//!   `Mars::east_positive_longitude` / `SurfaceSite::from_planetographic_longitude`.
//! - Moon frame: NAIF `IAU_MOON` analytical Mean-Earth (`MoonIauMeanEarth`)
//! - Epoch is TT; lunar orientation converts TT → TDB through tempoch
//! - Geometric only: no refraction, aberration, light-time, or terrain
//!
//! ```text
//! cargo run --example 21_body_surface_horizontal
//! ```
#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::{Earth, Mars, Moon};
use siderust::coordinates::cartesian::{Direction, Position};
use siderust::coordinates::centers::{Geocentric, Selenocentric};
use siderust::coordinates::frames::{EclipticMeanJ2000, ICRS};
use siderust::coordinates::surface::{
    BodySurfaceModel, SurfaceSite, EARTH_WGS84_SURFACE, MARS_IAU_ELLIPSOID, MOON_IAU_SPHERE,
};
use siderust::coordinates::transform::{
    direction_to_horizontal, horizontal_to_direction, position_to_horizontal, AstroContext,
    SurfaceHorizontalExt, TransformCenter,
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

    // ── 1. Moon: finite-distance Earth and Sun (parallax applied) ─────────
    let moon_site = SurfaceSite::<Moon>::new(20.0 * DEG, (-3.5) * DEG, 0.0 * M)?;
    println!("Moon observer");
    println!(
        "  site: lon={:.3} east, lat={:.3}, h={}, surface=IAU sphere R={} km",
        moon_site.lon(),
        moon_site.lat(),
        moon_site.height(),
        MOON_IAU_SPHERE
            .equatorial_radius()
            .to::<Kilometer>()
            .value()
    );
    println!(
        "  frame: MoonIauMeanEarth (NAIF IAU_MOON analytical); lat={}",
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

    // Convert the Sun's barycentric position to a Moon-centred position first.
    // Normalising sun_bary directly would yield the barycentre→Sun direction.
    let sun_bary = Vsop87Ephemeris::sun_barycentric(jd);
    let sun_sel: Position<Selenocentric, EclipticMeanJ2000, AstronomicalUnit> =
        TransformCenter::to_center_with(&sun_bary, (), jd, &ctx);
    let sun_obs = moon_site.observe_position(&sun_sel, &jd, &ctx)?;
    println!(
        "  target: Sun (selenocentric, finite)   alt={:.4}  az={:.4}  range={}\n",
        sun_obs.alt(),
        sun_obs.az(),
        sun_obs.range().to::<AstronomicalUnit>()
    );

    // ── 2. Mars: direction-only ICRS target (no invented parallax) ────────
    // Gale crater ≈ 137.4°E east-positive. If you have a SPICE planetographic
    // west longitude W, use Mars::east_positive_longitude(W) instead.
    let mars_site = SurfaceSite::<Mars>::new(137.4 * DEG, (-4.6) * DEG, 0.0 * M)?;
    println!("Mars observer (Gale-crater-like, lon east-positive)");
    let pgr_west = mars_site
        .planetographic_longitude()
        .value()
        .rem_euclid(360.0);
    println!(
        "  site: lon={} east (= {:.3}° west planetographic), lat={} ({}), a={} km",
        mars_site.lon(),
        pgr_west,
        mars_site.lat(),
        mars_site.latitude_convention().as_str(),
        MARS_IAU_ELLIPSOID
            .equatorial_radius()
            .to::<Kilometer>()
            .value()
    );
    println!(
        "  frame: MarsFixed (linear IAU α0,δ0,W); pgr lon sense={}",
        Mars::planetographic_longitude_sense().as_str()
    );

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

    // ── 3. Earth: same generic API with WGS84 (coexists with Topocentric) ─
    let earth_site = SurfaceSite::<Earth>::new((-17.89) * DEG, 28.75 * DEG, 2326.0 * M)?;
    println!("Earth observer (La Palma-like)");
    println!(
        "  site: lon={:.3}, lat={:.3}, h={}, surface=WGS84 a={:.0} m",
        earth_site.lon(),
        earth_site.lat(),
        earth_site.height(),
        EARTH_WGS84_SURFACE.equatorial_radius().value()
    );
    println!("  frame: ECEF (ITRF orientation chain via ICRS hub)");
    println!("  note: Earth-only GMST/GAST Topocentric APIs remain available");

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
