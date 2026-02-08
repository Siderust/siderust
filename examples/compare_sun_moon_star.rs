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
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::direction;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

/// Generic function that works for ANY body implementing AltitudePeriodsProvider
fn analyze_body<B: AltitudePeriodsProvider>(
    body: &B,
    name: &str,
    observer: ObserverSite,
    window: Period<ModifiedJulianDate>,
) {
    println!("─────────────────────────────────────────────────");
    println!("  {}", name);
    println!("─────────────────────────────────────────────────");

    // Above horizon
    let up_periods = body.above_threshold(observer, window, Degrees::new(0.0));
    let total_up: f64 = up_periods.iter().map(|p| p.duration_days() * 24.0).sum();
    println!(
        "  Above horizon: {:.1} hours in {} periods",
        total_up,
        up_periods.len()
    );

    // High altitude (> 45°)
    let high_periods = body.above_threshold(observer, window, Degrees::new(45.0));
    let total_high: f64 = high_periods.iter().map(|p| p.duration_days() * 24.0).sum();
    println!(
        "  Above 45°:     {:.1} hours in {} periods",
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
    let total_low: f64 = low_periods.iter().map(|p| p.duration_days() * 24.0).sum();
    println!(
        "  Low (0-10°):   {:.1} hours in {} periods",
        total_low,
        low_periods.len()
    );

    // Peak altitude during window
    let mid_jd =
        siderust::astro::JulianDate::new((window.start.value() + window.end.value()) / 2.0);
    let alt = body.altitude_at(&observer, mid_jd).to::<Degree>();
    println!("  Altitude at window midpoint: {:.2}°", alt.value());
    println!();
}

fn main() {
    println!("╔════════════════════════════════════════════════════════╗");
    println!("║    Generic Altitude Analysis: Sun, Moon, Star         ║");
    println!("╚════════════════════════════════════════════════════════╝\n");

    let observatory = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0));
    println!(
        "Observatory: Greenwich ({:.2}°N, {:.2}°E)",
        observatory.lat.value(),
        observatory.lon.value()
    );

    let window = Period::new(
        ModifiedJulianDate::new(60000.0),
        ModifiedJulianDate::new(60001.0),
    );
    println!(
        "Window: 24 hours (MJD {:.1} to {:.1})\n",
        window.start.value(),
        window.end.value()
    );

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
