// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Find Night Periods (365 days) Example (prefixed)
//!
//! Run with: `cargo run --example 07_find_night_periods_365day`

use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::night_types::twilight;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period, MJD};

fn main() {
    println!("Find astronomical night periods (prefixed example)");
}
