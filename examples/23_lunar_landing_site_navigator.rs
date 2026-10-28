// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar landing site navigator (NASA Space Apps 2026)
//!
//! Observe the Sun and Earth from a candidate lunar landing site, every six
//! hours over two days. Uses typed, finite-distance positions (parallax included).
//!
//! ```text
//! cargo run --example 23_lunar_landing_site_navigator
//! cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5
//! ```
//!
//! Arguments: planetocentric latitude (degrees), east-positive longitude
//! (degrees), and start Julian Date in **TT** (not UTC).
//!
//! Alt/Az uses the IAU lunar reference sphere and analytical IAU_MOON frame.
//! "Above horizon" refers to the **centre** of the Sun or Earth above the
//! geometric horizon. Terrain, solar-disc extent, aberration, light-time and
//! communications constraints are not modelled.

#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::Moon;
use siderust::coordinates::cartesian::position;
use siderust::coordinates::centers::Selenocentric;
use siderust::coordinates::surface::{SurfaceError, SurfaceSite};
use siderust::coordinates::transform::{
    AstroContext, HorizontalObservation, SurfaceHorizontalExt, TransformCenter,
};
use siderust::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::qtty::{AstronomicalUnit, DEG, M};
use siderust::time::JulianDate;

// Reuse Siderust/affn's typed positions instead of building Cartesian vectors.
type LunarTarget = position::EclipticMeanJ2000<AstronomicalUnit, Selenocentric>;
type LunarObservation = HorizontalObservation<Moon, AstronomicalUnit>;

fn observe(
    site: &SurfaceSite<Moon>,
    jd: JulianDate,
    ctx: &AstroContext,
) -> Result<(LunarObservation, LunarObservation), SurfaceError> {
    // Barycentric ephemerides -> Moon-centred positions; no manual XYZ arithmetic.
    let earth: LunarTarget =
        TransformCenter::to_center_with(&Vsop87Ephemeris::earth_barycentric(jd), (), jd, ctx);
    let sun: LunarTarget =
        TransformCenter::to_center_with(&Vsop87Ephemeris::sun_barycentric(jd), (), jd, ctx);

    Ok((
        site.observe_position(&sun, &jd, ctx)?,
        site.observe_position(&earth, &jd, ctx)?,
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (latitude, longitude, start_jd_tt) = match args.as_slice() {
        [] => (-89.0, 0.0, 2_461_358.5),
        [lat, lon, jd] => (lat.parse::<f64>()?, lon.parse::<f64>()?, jd.parse::<f64>()?),
        _ => return Err("Usage: [LAT_DEG LON_EAST_DEG START_JD_TT]".into()),
    };
    if !start_jd_tt.is_finite() || start_jd_tt <= 0.0 {
        return Err("Start JD(TT) must be a finite, positive Julian date".into());
    }

    let site = SurfaceSite::<Moon>::new(longitude * DEG, latitude * DEG, 0.0 * M)?;
    let ctx = AstroContext::default();

    println!("Lunar site: lat={}, east lon={} (JD in TT)", site.lat(), site.lon());
    println!("Geometric horizon: azimuth N=0°, E=90°; altitude > 0° means above.");

    for hours in (0..=48).step_by(6) {
        let jd = JulianDate::new(start_jd_tt + f64::from(hours) / 24.0);
        let (sun, earth) = observe(&site, jd, &ctx)?;

        println!(
            "T+{hours:02}h (JD {jd:.5}): Sun alt={:.2}, az={:.2}, above={}; Earth alt={:.2}, az={:.2}, above={}",
            sun.alt(), sun.az(), sun.alt() > 0.0 * DEG,
            earth.alt(), earth.az(), earth.alt() > 0.0 * DEG,
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_sun_and_earth_observations() {
        let site = SurfaceSite::<Moon>::new(0.0 * DEG, -89.0 * DEG, 0.0 * M).unwrap();
        let (sun, earth) = observe(
            &site,
            JulianDate::new(2_461_358.5),
            &AstroContext::default(),
        )
        .unwrap();

        for observation in [sun, earth] {
            assert!((-90.0 * DEG..=90.0 * DEG).contains(&observation.alt()));
            assert!((0.0 * DEG..360.0 * DEG).contains(&observation.az()));
        }
    }
}
