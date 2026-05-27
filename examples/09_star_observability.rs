// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Star observability in altitude + azimuth ranges.
//!
//! Run with: `cargo run --example 09_star_observability`
#![allow(clippy::print_stdout)]

use siderust::bodies::catalog::SIRIUS;
use siderust::event::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::event::azimuth::{AzimuthProvider, AzimuthQuery};
use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::intersect_periods;
use siderust::time::{Interval, ModifiedJulianDate};

use siderust::qtty::*;

fn main() {
    println!("Star observability: altitude + azimuth constraints\n");

    let observer = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let target = SIRIUS;

    // One-night search window (MJD TT).
    let t_0 = ModifiedJulianDate::new(60000.0);
    let window = Interval::new(t_0, t_0 + Days::new(1.0));

    // Constraint 1: altitude between 25° and 65°.
    let altitude_query = AltitudeQuery {
        observer,
        window,
        min_altitude: Degrees::new(25.0),
        max_altitude: Degrees::new(65.0),
        correction_policy: siderust::astro::apparent::CorrectionPolicy::APPARENT,
    };
    let altitude_periods = target.altitude_periods(&altitude_query);

    // Constraint 2: azimuth between 110° and 220° (ESE -> SW sector).
    let azimuth_query = AzimuthQuery {
        observer,
        window,
        min_azimuth: Degrees::new(110.0),
        max_azimuth: Degrees::new(220.0),
        opts: siderust::event::azimuth::SearchOpts::default(),
        correction_policy: siderust::astro::apparent::CorrectionPolicy::APPARENT,
    };
    let azimuth_periods = target.azimuth_periods(&azimuth_query);

    // Final observability: periods satisfying both constraints simultaneously.
    let observable_periods = intersect_periods(&altitude_periods, &azimuth_periods);

    println!("Observer: {}", ROQUE_DE_LOS_MUCHACHOS.name);
    println!("Target: Sirius");
    println!("Window: {}\n", window);

    println!(
        "Altitude range: {}..{}",
        altitude_query.min_altitude, altitude_query.max_altitude
    );
    println!(
        "Azimuth range:  {}..{}\n",
        azimuth_query.min_azimuth, azimuth_query.max_azimuth
    );

    println!("Matched periods: {}", observable_periods.len());
    for (idx, period) in observable_periods.iter().enumerate() {
        let hours = period.length().to::<Hour>();
        println!(
            "  {}. {} -> {}  ({})",
            idx + 1,
            period.start,
            period.end,
            hours
        );
    }

    let total_hours = observable_periods
        .iter()
        .fold(Hours::new(0.0), |acc, p| acc + p.length().to::<Hour>());
    println!("\nTotal observable time in both ranges: {}", total_hours);
}
