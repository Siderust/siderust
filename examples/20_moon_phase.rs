// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Moon Phase Example (prefixed)
//!
//! Run with: `cargo run --example 20_moon_phase`

use chrono::{NaiveDate, NaiveDateTime, NaiveTime, TimeZone, Utc};
use qtty::*;
use siderust::calculus::ephemeris::Vsop87Ephemeris;
use siderust::calculus::lunar::phase::{
    find_phase_events, illumination_range, moon_phase_geocentric, moon_phase_topocentric,
    PhaseSearchOpts,
};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::{JulianDate, ModifiedJulianDate, Period};

fn main() {
    println!("Moon phase (prefixed)");
}
