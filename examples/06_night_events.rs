// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Night Events Example (prefixed)
//!
//! Shows how to spot "night-type" crossing events and night periods in a
//! one-week window using civil/nautical/astronomical/horizon thresholds.
//!
//! Run with:
//! `cargo run --example 06_night_events -- [YYYY-MM-DD] [lat_deg] [lon_deg] [height_m]`

use chrono::{NaiveDate, NaiveDateTime, TimeZone, Utc};
use qtty::{Days, Degrees, Meter, Quantity};
use siderust::bodies::Sun;
use siderust::calculus::altitude::{below_threshold, crossings, CrossingDirection, SearchOpts};
use siderust::calculus::solar::night_types::{twilight, Twilight};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::ModifiedJulianDate;
use siderust::time::{Period, MJD};

fn week_period_from_date(start_date: NaiveDate) -> Period<MJD> {
    let start_dt = Utc.from_utc_datetime(&NaiveDateTime::new(
        start_date,
        chrono::NaiveTime::from_hms_opt(0, 0, 0).expect("00:00:00 is valid"),
    ));
    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = mjd_start + Days::new(7.0);
    Period::new(mjd_start, mjd_end)
}

fn print_events_for_type(site: &Geodetic<ECEF>, week: Period<MJD>, name: &str, threshold: Degrees) {
    let events = crossings(&Sun, site, week, threshold, SearchOpts::default());
    let mut downs = 0usize;
    let mut raises = 0usize;

    println!(
        "{name:18} threshold {threshold:>8} -> {:2} crossing(s)",
        events.len()
    );

    for ev in &events {
        let label = match ev.direction {
            // Sun going below threshold: night-type starts/intensifies.
            CrossingDirection::Setting => {
                downs += 1;
                "night-type down (Sun setting below threshold)"
            }
            // Sun going above threshold: night-type ends/weakens.
            CrossingDirection::Rising => {
                raises += 1;
                "night-type raise (Sun rising above threshold)"
            }
        };
        if let Some(t_utc) = ev.mjd.to_utc() {
            println!("  - {} at {}", label, t_utc.format("%Y-%m-%dT%H:%M:%S"));
        }
    }

    println!("  summary: down={} raise={}", downs, raises);
}

fn print_periods_for_type(
    site: &Geodetic<ECEF>,
    week: Period<MJD>,
    name: &str,
    threshold: Degrees,
) {
    let periods = below_threshold(&Sun, site, week, threshold, SearchOpts::default());
    println!(
        "{name:18} night periods (Sun < {threshold}): {}",
        periods.len()
    );

    for p in periods {
        if let (Some(s), Some(e)) = (p.start.to_utc(), p.end.to_utc()) {
            println!(
                "  - {} -> {} ({:.1} h)",
                s.format("%Y-%m-%dT%H:%M:%S"),
                e.format("%Y-%m-%dT%H:%M:%S"),
                p.duration_days().to::<qtty::Hour>()
            );
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

    let lat_deg = args
        .get(2)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(51.4769);
    let lon_deg = args
        .get(3)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let height_m = args
        .get(4)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);

    let site = Geodetic::<ECEF>::new(
        Degrees::new(lon_deg),
        Degrees::new(lat_deg),
        Quantity::<Meter>::new(height_m),
    );
    let week = week_period_from_date(start_date);

    let night_types = [
        ("Horizon", Degrees::from(Twilight::Horizon)),
        ("Apparent Horizon", Degrees::from(Twilight::ApparentHorizon)),
        ("Civil", twilight::CIVIL),
        ("Nautical", twilight::NAUTICAL),
        ("Astronomical", twilight::ASTRONOMICAL),
    ];

    println!("Night events over one week");
    println!("==========================");
    println!(
        "Site: lat={} lon={} height={}",
        site.lat, site.lon, site.height
    );
    println!("Week start: {} UTC\n", start_date);

    println!("1) Night-type crossing events");
    for (name, thr) in night_types {
        print_events_for_type(&site, week, name, thr);
    }

    println!("\n2) Night periods per night type");
    for (name, thr) in night_types {
        print_periods_for_type(&site, week, name, thr);
    }
}
