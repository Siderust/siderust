// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Culmination-Based Solar Altitude Example
//! Demonstrates how to use `find_sun_altitude_periods_via_culminations` to
//! compute astronomical night (Sun altitude < -18°) for a whole year at the
//! Roque de los Muchachos observatory.
//!
//! ## Usage
//! ```
//! cargo run --example solar_altitude_culminations
//! ```

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use siderust::calculus::solar::altitude_periods::find_night_periods;
use siderust::calculus::solar::night_types::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};

fn main() {
    const OBSERVED_YEAR: i32 = 2026;

    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let start_naive = NaiveDate::from_ymd_opt(OBSERVED_YEAR, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = NaiveDate::from_ymd_opt(OBSERVED_YEAR + 1, 1, 1)
        .unwrap()
        .and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);

    let period = Period::new(mjd_start, mjd_end);

    println!("Roque de los Muchachos (La Palma) – astronomical night 2026");
    println!("Interval: {} → {} UTC", start_dt.format("%Y-%m-%d"), end_dt.format("%Y-%m-%d"));

    let periods = find_night_periods(site, period, twilight::ASTRONOMICAL);
    if !periods.is_empty() {
        let total_days: f64 = periods.iter().map(|p| p.duration_days()).sum();
        println!("Found {} astronomical-night windows spanning {:.1} hours", periods.len(), total_days * 24.0);

        for (idx, period) in periods.iter().take(5).enumerate() {
            let start = period.start.to_utc().unwrap();
            let end = period.end.to_utc().unwrap();
            println!(
                "  #{:<2} {:<19} → {:<19} ({:.1} h)",
                idx + 1,
                start.format("%Y-%m-%dT%H:%M:%S"),
                end.format("%Y-%m-%dT%H:%M:%S"),
                period.duration_days() * 24.0
            );
        }

        if periods.len() > 5 {
            println!("  ... ({} more nights)", periods.len() - 5);
        }
    } else {
        println!("No astronomical night windows found on the selected interval.");
    }
}
