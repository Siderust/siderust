// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Periods Example (prefixed)
//!
//! Run with: `cargo run --example 05_time_periods`

use chrono::DateTime;
use siderust::time::{JulianDate, ModifiedJulianDate, Period, UtcPeriod, MJD};

fn main() {
    println!("Generic Time Period Examples");
    println!("============================\n");

    // Example 1: Period with JulianDate
    println!("1. Period with JulianDate:");
    let jd_start = JulianDate::new(2451545.0); // J2000.0
    let jd_end = JulianDate::new(2451546.5); // 1.5 days later
    let jd_period = Period::new(jd_start, jd_end);
    println!("   Start: {}", jd_start);
    println!("   End:   {}", jd_end);
    println!("   Duration: {} days\n", jd_period.duration_days());

    // ... (rest unchanged)
}
