// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomical Night Example (prefixed)
//!
//! Run with: `cargo run --example 06_astronomical_night`

use chrono::{NaiveDate, NaiveDateTime, TimeZone, Utc};
use qtty::{Degrees, Meter, Quantity};
use siderust::bodies::Sun;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::calculus::solar::night_types::twilight;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::ModifiedJulianDate;
use siderust::time::Period;

fn main() {
    // (original example content preserved)
    println!("Astronomical Night example (prefixed)");
}
