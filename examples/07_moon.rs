// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Moon phase quick examples.
//!
//! Shows:
//! 1) how to get moon phase properties at a given instant,
//! 2) how to find windows where illumination is in a given range.
//!
//! Run with:
//! `cargo run --example 07_moon_phase -- [YYYY-MM-DD] [lat_deg] [lon_deg] [height_m]`

use chrono::{NaiveDate, NaiveDateTime, NaiveTime, TimeZone, Utc};
use qtty::{Days, Degree, Degrees, Meter, Quantity};
use siderust::calculus::ephemeris::Vsop87Ephemeris;
use siderust::calculus::lunar::phase::{
    find_phase_events, illumination_range, moon_phase_geocentric, moon_phase_topocentric,
    PhaseSearchOpts,
};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::{JulianDate, ModifiedJulianDate, Period, MJD};

fn print_periods(label: &str, periods: &[Period<MJD>]) {
    println!("\n{label}: {} period(s)", periods.len());
    for p in periods {
        let dur_h = p.duration_days().to::<qtty::Hour>();
        match (p.start.to_utc(), p.end.to_utc()) {
            (Some(s), Some(e)) => println!(
                "  - {} -> {} ({dur_h})",
                s.format("%Y-%m-%d %H:%M UTC"),
                e.format("%Y-%m-%d %H:%M UTC")
            ),
            _ => println!("  - MJD {} -> {} ({dur_h})", p.start, p.end),
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let today = Utc::now().date_naive();

    let start_date = args
        .get(1)
        .and_then(|s| NaiveDate::parse_from_str(s, "%Y-%m-%d").ok())
        .unwrap_or(today);
    let lat = args.get(2).and_then(|s| s.parse::<f64>().ok()).unwrap_or(28.762);
    let lon = args
        .get(3)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(-17.892);
    let h_m = args
        .get(4)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(2396.0);

    let site = Geodetic::<ECEF>::new(
        Degrees::new(lon),
        Degrees::new(lat),
        Quantity::<Meter>::new(h_m),
    );

    let midnight = NaiveDateTime::new(
        start_date,
        NaiveTime::from_hms_opt(0, 0, 0).expect("00:00:00 must be valid"),
    );
    let jd = JulianDate::from_utc(Utc.from_utc_datetime(&midnight));
    let mjd = ModifiedJulianDate::from(jd);
    let window = Period::new(mjd, mjd + Days::new(35.0));
    let opts = PhaseSearchOpts::default();

    // 1) Point-in-time phase properties.
    let geo = moon_phase_geocentric::<Vsop87Ephemeris>(jd);
    let topo = moon_phase_topocentric::<Vsop87Ephemeris>(jd, site);

    println!("Moon phase at {} 00:00 UTC", start_date);
    println!("==================================");
    println!("Site: lat={lat:.4} deg, lon={lon:.4} deg, h={h_m:.0} m");
    println!("\nGeocentric:");
    println!("  label                 : {}", geo.label());
    println!("  illuminated fraction  : {:.4}", geo.illuminated_fraction);
    println!("  illuminated percent   : {:.2} %", geo.illuminated_percent());
    println!("  phase angle           : {}", geo.phase_angle.to::<Degree>());
    println!("  elongation            : {}", geo.elongation.to::<Degree>());
    println!("  waxing                : {}", geo.waxing);

    println!("\nTopocentric:");
    println!("  label                 : {}", topo.label());
    println!("  illuminated fraction  : {:.4}", topo.illuminated_fraction);
    println!(
        "  illumination delta    : {:+.4} %",
        (topo.illuminated_fraction - geo.illuminated_fraction) * 100.0
    );
    println!("  elongation            : {}", topo.elongation.to::<Degree>());

    // 2) Principal phase events (optional but useful alongside range searches).
    let events = find_phase_events::<Vsop87Ephemeris>(window, opts);
    println!("\nPrincipal phase events in next 35 days: {}", events.len());
    for ev in &events {
        if let Some(utc) = ev.mjd.to_utc() {
            println!("  - {:>13} at {}", ev.kind, utc.format("%Y-%m-%d %H:%M UTC"));
        } else {
            println!("  - {:>13} at MJD {}", ev.kind, ev.mjd);
        }
    }

    // 3) Find periods where Moon illumination is in requested phase ranges.
    let crescent = illumination_range::<Vsop87Ephemeris>(window, 0.05, 0.35, opts);
    let quarterish = illumination_range::<Vsop87Ephemeris>(window, 0.45, 0.55, opts);
    let gibbous = illumination_range::<Vsop87Ephemeris>(window, 0.65, 0.95, opts);

    print_periods("Crescent-like range (5%-35%)", &crescent);
    print_periods("Quarter-like range (45%-55%)", &quarterish);
    print_periods("Gibbous-like range (65%-95%)", &gibbous);
}
