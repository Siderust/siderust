// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Night Quality Scoring (prefixed)
//!
//! Run with: `cargo run --example 08_night_quality_scoring`

use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

fn main() {
    println!("Night quality scoring (prefixed)");
}
