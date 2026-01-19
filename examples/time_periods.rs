// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example demonstrating the generic Period<T> with different time types
//!
//! This example shows how to use Period with JulianDate, ModifiedJulianDate,
//! and DateTime<Utc>.

use chrono::DateTime;
use siderust::time::{JulianDate, ModifiedJulianDate, Period};

fn main() {
    println!("Generic Time Period Examples");
    println!("============================\n");

    // Example 1: Period with JulianDate
    println!("1. Period with JulianDate:");
    let jd_start = JulianDate::new(2451545.0); // J2000.0
    let jd_end = JulianDate::new(2451546.5); // 1.5 days later
    let jd_period = Period::new(jd_start, jd_end);
    println!("   Start: JD {}", jd_start.value());
    println!("   End:   JD {}", jd_end.value());
    println!("   Duration: {} days\n", jd_period.duration_days());

    // Example 2: Period with ModifiedJulianDate
    println!("2. Period with ModifiedJulianDate:");
    let mjd_start = ModifiedJulianDate::new(59000.0);
    let mjd_end = ModifiedJulianDate::new(59002.5);
    let mjd_period = Period::new(mjd_start, mjd_end);
    println!("   Start: MJD {}", mjd_start.value());
    println!("   End:   MJD {}", mjd_end.value());
    println!("   Duration: {} days\n", mjd_period.duration_days());

    // Example 3: Period with UTC DateTime
    println!("3. Period with UTC DateTime:");
    let utc_start = DateTime::from_timestamp(0, 0).unwrap(); // Unix epoch
    let utc_end = DateTime::from_timestamp(86400 * 2, 0).unwrap(); // 2 days later
    let utc_period = Period::new(utc_start, utc_end);
    println!("   Start: {}", utc_start);
    println!("   End:   {}", utc_end);
    println!("   Duration: {} days", utc_period.duration_days());
    println!("   Duration: {} seconds\n", utc_period.duration_seconds());

    // Example 4: Converting between time systems
    println!("4. Converting between time systems:");
    let mjd = ModifiedJulianDate::new(51544.5); // MJD at J2000.0
    let jd = mjd.to_julian_day();
    let utc = mjd.to_utc().unwrap();
    println!("   MJD: {}", mjd.value());
    println!("   JD:  {}", jd.value());
    println!("   UTC: {}\n", utc);

    // Example 5: Backward compatibility with AltitudePeriod type alias
    println!("5. AltitudePeriod (type alias for Period<ModifiedJulianDate>):");
    use siderust::calculus::events::AltitudePeriod;
    let night_period = AltitudePeriod::new(
        ModifiedJulianDate::new(59000.0),
        ModifiedJulianDate::new(59000.5),
    );
    println!("   Start: {}", night_period.start);
    println!("   End:   {}", night_period.end);
    println!("   Duration: {} hours", night_period.duration_days() * 24.0);
}
