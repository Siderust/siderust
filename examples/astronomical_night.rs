// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical Night Example
//!
//! Demonstrates finding astronomical night periods using the siderust library.
//!
//! Astronomical night is defined as the period when the Sun's center is
//! more than 18° below the horizon (altitude < -18°).
//!
//! ## Usage
//!
//! ```bash
//! cargo run --example astronomical_night -- [YYYY-MM-DD] [lat_deg] [lon_deg] [height_m]
//! ```
//!
//! ## Defaults
//! - Start date: today UTC
//! - Location: Greenwich Observatory (51.4769°N, 0°E)
//! - Search period: 7 days

use chrono::{NaiveDate, NaiveDateTime, TimeZone, Utc};
use qtty::{Degrees, Meter, Quantity};
use siderust::time::ModifiedJulianDate;
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::night_types::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::time::Period;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // Parse command-line arguments with defaults
    let today = Utc::now().date_naive();
    let start_date = args
        .get(1)
        .map(|s| NaiveDate::parse_from_str(s, "%Y-%m-%d").expect("Invalid date format"))
        .unwrap_or(today);

    let lat_deg: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(51.4769);
    let lon_deg: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let height_m: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);

    // Create observer site
    let site = ObserverSite::new(
        Degrees::new(lon_deg),
        Degrees::new(lat_deg),
        Quantity::<Meter>::new(height_m),
    );

    // Convert start date to MJD and compute 7-day interval
    let datetime_start = Utc.from_utc_datetime(&NaiveDateTime::new(
        start_date,
        chrono::NaiveTime::from_hms_opt(0, 0, 0).unwrap(),
    ));
    let mjd_start = ModifiedJulianDate::from_utc(datetime_start);
    let mjd_end = ModifiedJulianDate::new(mjd_start.value() + 7.0);

    // Find astronomical night periods (Sun altitude < -18°)
    let period = Period::new(mjd_start, mjd_end);
    let nights = Sun.below_threshold(site, period, twilight::ASTRONOMICAL);

    // Print results
    println!("Astronomical Night Periods (Sun altitude < -18°)");
    println!("================================================");
    println!(
        "Observer: lat = {:.4}°, lon = {:.4}°, height = {} m",
        lat_deg, lon_deg, height_m
    );
    println!("Week starting: {} UTC", start_date);
    println!();

    if !nights.is_empty() {
        for period in nights {
            let start_utc = period.start.to_utc();
            let end_utc = period.end.to_utc();

            if let (Some(s), Some(e)) = (start_utc, end_utc) {
                let duration_mins = (period.duration_days().value() * 24.0 * 60.0).round() as i64;
                println!(
                    "{} → {}  ({} min)",
                    s.format("%Y-%m-%dT%H:%M:%S"),
                    e.format("%Y-%m-%dT%H:%M:%S"),
                    duration_mins
                );
            }
        }
    } else {
        println!("No astronomical night periods found in this week.");
        println!("(This can happen at high latitudes during summer.)");
    }
}
