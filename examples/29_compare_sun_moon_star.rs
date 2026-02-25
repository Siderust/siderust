// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Side-by-Side Comparison: Sun, Moon, and Star
//!
//! Demonstrates that the `AltitudePeriodsProvider` trait provides a truly
//! unified interface — the same code works for any celestial body.
//!
//! Run with: `cargo run --example compare_sun_moon_star`

use siderust::bodies::catalog::SIRIUS;
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::coordinates::spherical::direction;
use siderust::time::*;

use qtty::*;

/// Generic function that works for ANY body implementing AltitudePeriodsProvider
fn analyze_body<B: AltitudePeriodsProvider>(
    body: &B,
    name: &str,
    observer: Geodetic<ECEF>,
    window: Period<MJD>,
) {
    println!("─────────────────────────────────────────────────");
    println!("  {}", name);
    println!("─────────────────────────────────────────────────");

    // Above horizon
    let up_periods = body.above_threshold(observer, window, Degrees::new(0.0));
    let total_up = up_periods.iter().fold(Hours::new(0.0), |acc, p| {
        acc + p.duration_days().to::<Hour>()
    });
    println!(
        "  Above horizon: {} in {} periods",
        total_up,
        up_periods.len()
    );

    // High altitude (> 45°)
    let high_periods = body.above_threshold(observer, window, Degrees::new(45.0));
    let total_high = high_periods.iter().fold(Hours::new(0.0), |acc, p| {
        acc + p.duration_days().to::<Hour>()
    });
    println!(
        "  Above 45°:     {} in {} periods",
        total_high,
        high_periods.len()
    );

    // Twilight/low band (0° to 10°)
    let query = AltitudeQuery {
        observer,
        window,
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(10.0),
    };
    let low_periods = body.altitude_periods(&query);
    let total_low = low_periods.iter().fold(Hours::new(0.0), |acc, p| {
        acc + p.duration_days().to::<Hour>()
    });
    println!(
        "  Low (0-10°):   {} in {} periods",
        total_low,
        low_periods.len()
    );

    // Peak altitude during window
    let mid_mjd =
        siderust::time::ModifiedJulianDate::new((window.start.value() + window.end.value()) / 2.0);
    let alt = body.altitude_at(&observer, mid_mjd).to::<Degree>();
    println!("  Altitude at window midpoint: {}", alt);
    println!();
}

fn main() {
    println!("╔════════════════════════════════════════════════════════╗");
    println!("║    Generic Altitude Analysis: Sun, Moon, Star         ║");
    println!("╚════════════════════════════════════════════════════════╝\n");

    let observatory =
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0));
    println!(
        "Observatory: Greenwich ({}, {})",
        observatory.lat, observatory.lon
    );

    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    println!("Window: 24 hours ({} to {})\n", window.start, window.end);

    println!("═════════════════════════════════════════════════════════");

    // Same analysis for three different body types using ONE generic function
    analyze_body(&Sun, "Sun", observatory, window);
    analyze_body(&Moon, "Moon", observatory, window);
    analyze_body(&SIRIUS, "Sirius (α CMa)", observatory, window);
    analyze_body(
        &direction::ICRS::new(Degrees::new(279.23), Degrees::new(38.78)),
        "Vega (via direction::ICRS)",
        observatory,
        window,
    );

    println!("═════════════════════════════════════════════════════════");
    println!("\nKey Insight:");
    println!("  The SAME generic function analyzed all four bodies.");
    println!("  Sun/Moon delegate to VSOP87/ELP2000 engines.");
    println!("  Stars use the analytical sinusoidal model.");
    println!("  direction::ICRS is the lightweight path (RA/Dec only).");
    println!();
    println!("  All share the AltitudePeriodsProvider trait interface!");
    println!("═════════════════════════════════════════════════════════");
}
