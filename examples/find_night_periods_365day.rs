// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example: Find Night Periods for a Full Year
//!
//! Runs `find_night_periods` for a 365-day horizon and prints the
//! astronomical night periods (Sun altitude < -18°) for the chosen site.
//!
//! Usage:
//!
//! ```bash
//! cargo run --example find_night_periods_365day -- [YYYY-MM-DD]
//! ```
//!
//! Defaults:
//! - Start date: 2026-01-01 UTC
//! - Site: Roque de los Muchachos Observatory

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::night_types::twilight;
use siderust::coordinates::centers::ObserverSite;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};

fn build_period(start_date: NaiveDate, days: u32) -> Period<ModifiedJulianDate> {
    let start_naive = start_date.and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let end_naive = start_naive + chrono::Duration::days(days as i64);

    let start_dt = Utc.from_utc_datetime(&start_naive);
    let end_dt = Utc.from_utc_datetime(&end_naive);

    let mjd_start = ModifiedJulianDate::from_utc(start_dt);
    let mjd_end = ModifiedJulianDate::from_utc(end_dt);

    Period::new(mjd_start, mjd_end)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let start_date = args
        .get(1)
        .and_then(|s| NaiveDate::parse_from_str(s, "%Y-%m-%d").ok())
        .unwrap_or_else(|| NaiveDate::from_ymd_opt(2026, 1, 1).unwrap());

    println!(
        "Find astronomical night periods for 365 days starting {} UTC",
        start_date
    );
    println!("Observer: Roque de los Muchachos Observatory (La Palma)");

    let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);

    let period = build_period(start_date, 365);

    let nights = Sun.below_threshold(site, period, twilight::ASTRONOMICAL);

    if !nights.is_empty() {
        println!("Found {} night periods:\n", nights.len());
        for p in nights {
            if let (Some(s), Some(e)) = (p.start.to_utc(), p.end.to_utc()) {
                let mins = (p.duration_days() * 24.0 * 60.0).round() as i64;
                println!(
                    "{} → {}  ({} min)",
                    s.format("%Y-%m-%dT%H:%M:%S"),
                    e.format("%Y-%m-%dT%H:%M:%S"),
                    mins
                );
            }
        }
    } else {
        println!("No astronomical night periods found for this year at this site.");
    }
}
