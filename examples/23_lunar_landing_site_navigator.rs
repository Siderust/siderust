// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar landing site navigator: geometric Sun/Earth horizon coordinates
//!
//! This self-contained CLI example supports the NASA Space Apps Challenge 2026
//! "Lunar Landing Site Navigator" use case. It samples the Sun and Earth above
//! a selected lunar site, reporting azimuth, altitude, and a geometric
//! above-horizon flag. Use CSV output to plot or compare sites and dates.
//!
//! Coordinates: planetocentric latitude in degrees, east-positive longitude
//! in degrees, and height zero on the IAU mean lunar sphere (1737.4 km).
//! Times are Julian dates in TT, NOT UTC. Azimuth increases north through east.
//!
//! Quick start (no downloads or additional features):
//!
//!   cargo run --example 23_lunar_landing_site_navigator
//!   cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5 48 60
//!   cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5 48 60 --csv
//!
//! This is a GEOMETRIC reference-sphere screening tool, not an operational
//! illumination or communications prediction. It omits local terrain, lunar
//! topography, antenna masks, solar-disc size, light-time and aberration.
//! Sun above horizon does not directly imply usable solar power; Earth above
//! horizon does not guarantee a communications link.

#![allow(clippy::print_stdout, clippy::print_stderr)]

use siderust::bodies::solar_system::Moon;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Selenocentric;
use siderust::coordinates::frames::EclipticMeanJ2000;
use siderust::coordinates::surface::SurfaceSite;
use siderust::coordinates::transform::{AstroContext, SurfaceHorizontalExt, TransformCenter};
use siderust::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::qtty::{AstronomicalUnit, Kilometer, DEG, M};
use siderust::time::JulianDate;

const DEFAULT_LAT_DEG: f64 = -89.0;
const DEFAULT_LON_DEG: f64 = 0.0;
const DEFAULT_START_JD_TT: f64 = 2_461_358.5;
const DEFAULT_HOURS: f64 = 48.0;
const DEFAULT_STEP_MINUTES: f64 = 60.0;
const MAX_INTERVALS: usize = 10_000;

#[derive(Debug, Clone, Copy)]
struct Options {
    lat_deg: f64,
    lon_deg: f64,
    start_jd_tt: f64,
    duration_minutes: f64,
    step_minutes: f64,
    intervals: usize,
    csv: bool,
}

#[derive(Debug, Clone, Copy)]
struct HorizonAngles {
    az_deg: f64,
    alt_deg: f64,
}

fn usage() -> &'static str {
    "Usage: cargo run --example 23_lunar_landing_site_navigator -- [LAT_DEG LON_EAST_DEG START_JD_TT [HOURS STEP_MINUTES]] [--csv]\n\
     Defaults: lat=-89, lon=0, JD(TT)=2461358.5, hours=48, step=60 minutes.\n\
     Latitude: -90..90; longitude: -180..180, positive east. JD is TT, NOT UTC.\n\
     --csv emits only a machine-readable header and rows; --help shows this message."
}

fn finite_number(input: &str, name: &str) -> Result<f64, String> {
    let value = input
        .parse::<f64>()
        .map_err(|_| format!("{name} must be a number"))?;
    if !value.is_finite() {
        return Err(format!("{name} must be finite"));
    }
    Ok(value)
}

fn parse_args(args: &[String]) -> Result<Options, String> {
    let csv = args.last().is_some_and(|arg| arg == "--csv");
    let positional = if csv { &args[..args.len() - 1] } else { args };
    let (lat_deg, lon_deg, start_jd_tt, hours, step_minutes) = match positional {
        [] => (
            DEFAULT_LAT_DEG,
            DEFAULT_LON_DEG,
            DEFAULT_START_JD_TT,
            DEFAULT_HOURS,
            DEFAULT_STEP_MINUTES,
        ),
        [lat, lon, start] => (
            finite_number(lat, "latitude")?,
            finite_number(lon, "longitude")?,
            finite_number(start, "start JD(TT)")?,
            DEFAULT_HOURS,
            DEFAULT_STEP_MINUTES,
        ),
        [lat, lon, start, hours, step] => (
            finite_number(lat, "latitude")?,
            finite_number(lon, "longitude")?,
            finite_number(start, "start JD(TT)")?,
            finite_number(hours, "hours")?,
            finite_number(step, "step minutes")?,
        ),
        _ => return Err(usage().to_owned()),
    };

    if !(-90.0..=90.0).contains(&lat_deg) {
        return Err("latitude must lie in [-90, 90] degrees".to_owned());
    }
    if !(-180.0..=180.0).contains(&lon_deg) {
        return Err("longitude must lie in [-180, 180] degrees (positive east)".to_owned());
    }
    if start_jd_tt <= 0.0 {
        return Err("start JD(TT) must be positive".to_owned());
    }
    if hours <= 0.0 || step_minutes <= 0.0 {
        return Err("hours and step minutes must be positive".to_owned());
    }

    let duration_minutes = hours * 60.0;
    let interval_count = (duration_minutes / step_minutes).ceil();
    if !duration_minutes.is_finite()
        || !interval_count.is_finite()
        || interval_count > MAX_INTERVALS as f64
    {
        return Err(format!(
            "requested time range is too large (maximum {MAX_INTERVALS} intervals)"
        ));
    }
    if !(start_jd_tt + duration_minutes / 1440.0).is_finite() {
        return Err("end JD(TT) must be finite".to_owned());
    }

    Ok(Options {
        lat_deg,
        lon_deg,
        start_jd_tt,
        duration_minutes,
        step_minutes,
        intervals: interval_count as usize,
        csv,
    })
}

fn observe_sun_and_earth(
    site: &SurfaceSite<Moon>,
    jd: JulianDate,
    ctx: &AstroContext,
) -> Result<(HorizonAngles, HorizonAngles), Box<dyn std::error::Error>> {
    // The Earth viewed from the lunar centre is the negative of the
    // geocentric Moon vector. Keep the finite range for lunar parallax.
    let moon_geocentric = Vsop87Ephemeris::moon_geocentric(jd);
    let earth_selenocentric = Position::<Selenocentric, EclipticMeanJ2000, Kilometer>::new(
        -moon_geocentric.x(),
        -moon_geocentric.y(),
        -moon_geocentric.z(),
    );
    let earth = site.observe_position(&earth_selenocentric, &jd, ctx)?;

    // Convert Sun barycentric -> selenocentric before forming the sightline.
    let sun_barycentric = Vsop87Ephemeris::sun_barycentric(jd);
    let sun_selenocentric: Position<Selenocentric, EclipticMeanJ2000, AstronomicalUnit> =
        TransformCenter::to_center_with(&sun_barycentric, (), jd, ctx);
    let sun = site.observe_position(&sun_selenocentric, &jd, ctx)?;

    Ok((
        HorizonAngles {
            az_deg: sun.az().value(),
            alt_deg: sun.alt().value(),
        },
        HorizonAngles {
            az_deg: earth.az().value(),
            alt_deg: earth.alt().value(),
        },
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!("{}", usage());
        return Ok(());
    }
    let options = parse_args(&args)
        .map_err(|message| std::io::Error::new(std::io::ErrorKind::InvalidInput, message))?;
    let site = SurfaceSite::<Moon>::new(
        options.lon_deg * DEG,
        options.lat_deg * DEG,
        0.0 * M,
    )?;
    let ctx = AstroContext::default();

    if options.csv {
        println!(
            "jd_tt,sun_az_deg,sun_alt_deg,sun_above_horizon,earth_az_deg,earth_alt_deg,earth_above_horizon"
        );
    } else {
        println!("Lunar landing site navigator (geometric reference horizon)");
        println!(
            "Site: latitude {:.4} deg, east longitude {:.4} deg, height 0 m",
            options.lat_deg, options.lon_deg
        );
        println!(
            "Start JD(TT): {:.6}; duration: {:.2} h; step: {:.2} min",
            options.start_jd_tt,
            options.duration_minutes / 60.0,
            options.step_minutes
        );
        println!("Azimuth: N=0, E=90; altitude: horizon=0, zenith=+90 degrees.");
        println!("  JD(TT)       Sun az   Sun alt  Sun>0  Earth az Earth alt Earth>0");
    }

    let mut sun_above_samples = 0usize;
    let mut earth_above_samples = 0usize;
    // Include both endpoints. The final interval may be shorter than step_minutes.
    for i in 0..=options.intervals {
        let elapsed_minutes =
            (i as f64 * options.step_minutes).min(options.duration_minutes);
        let jd_tt = options.start_jd_tt + elapsed_minutes / 1440.0;
        let (sun, earth) = observe_sun_and_earth(&site, JulianDate::new(jd_tt), &ctx)?;
        let sun_above = sun.alt_deg > 0.0;
        let earth_above = earth.alt_deg > 0.0;
        sun_above_samples += usize::from(sun_above);
        earth_above_samples += usize::from(earth_above);

        if options.csv {
            println!(
                "{jd_tt:.8},{:.6},{:.6},{},{:.6},{:.6},{}",
                sun.az_deg,
                sun.alt_deg,
                u8::from(sun_above),
                earth.az_deg,
                earth.alt_deg,
                u8::from(earth_above)
            );
        } else {
            println!(
                "{jd_tt:>13.6} {:>8.2} {:>9.2} {:>6} {:>9.2} {:>9.2} {:>7}",
                sun.az_deg,
                sun.alt_deg,
                if sun_above { "yes" } else { "no" },
                earth.az_deg,
                earth.alt_deg,
                if earth_above { "yes" } else { "no" }
            );
        }
    }

    if !options.csv {
        let samples = options.intervals + 1;
        println!(
            "Samples above geometric horizon: Sun {sun_above_samples}/{samples}, Earth {earth_above_samples}/{samples}."
        );
        println!(
            "These sample counts are not duration-weighted windows or mission feasibility metrics."
        );
        println!(
            "No topographic horizon, partial solar disc, antenna mask, light-time, or aberration."
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options(input: &[&str]) -> Result<Options, String> {
        parse_args(&input.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }

    #[test]
    fn defaults_are_south_polar_and_bounded() {
        let cfg = options(&[]).unwrap();
        assert_eq!(cfg.lat_deg, -89.0);
        assert_eq!(cfg.intervals, 48);
        assert!(!cfg.csv);
    }

    #[test]
    fn custom_sampling_includes_a_short_final_interval() {
        let cfg = options(&["-89.5", "0", "2461358.5", "1.5", "60", "--csv"]).unwrap();
        assert_eq!(cfg.intervals, 2);
        assert_eq!(cfg.duration_minutes, 90.0);
        assert!(cfg.csv);
        assert_eq!(
            (cfg.intervals as f64 * cfg.step_minutes).min(cfg.duration_minutes),
            90.0
        );
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(options(&["-91", "0", "2461358.5"]).is_err());
        assert!(options(&["0", "181", "2461358.5"]).is_err());
        assert!(options(&["0", "0", "NaN"]).is_err());
        assert!(options(&["0", "0", "2461358.5", "0", "60"]).is_err());
        assert!(options(&["0", "0", "2461358.5", "48", "0"]).is_err());
        assert!(options(&["0", "0", "2461358.5", "48", "0.00001"]).is_err());
        assert!(options(&["--csv", "0", "0"]).is_err());
    }

    #[test]
    fn real_sun_and_earth_observations_are_finite() {
        let site = SurfaceSite::<Moon>::new(0.0 * DEG, -89.0 * DEG, 0.0 * M).unwrap();
        let (sun, earth) = observe_sun_and_earth(
            &site,
            JulianDate::new(DEFAULT_START_JD_TT),
            &AstroContext::default(),
        )
        .unwrap();
        for angles in [sun, earth] {
            assert!(angles.alt_deg.is_finite());
            assert!((-90.0..=90.0).contains(&angles.alt_deg));
            assert!((0.0..360.0).contains(&angles.az_deg));
        }
    }
}
