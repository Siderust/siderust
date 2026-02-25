// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Altitude Periods — Trait-Based API Example (prefixed)
//!
//! Run with: `cargo run --example 10_altitude_periods_trait`

use siderust::bodies::catalog::{POLARIS, SIRIUS, VEGA};
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::spherical::direction;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

fn main() {
    println!("=== Altitude Periods API Examples (prefixed) ===\n");
    // (rest of example preserved)
}
