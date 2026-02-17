// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Star Observability Planner
//!
//! Demonstrates using the `AltitudePeriodsProvider` trait to plan optimal
//! observing windows for multiple stars at a given observatory.
//!
//! Run with: `cargo run --example star_observability`

use siderust::bodies::catalog::{ALTAIR, BETELGEUSE, POLARIS, RIGEL, SIRIUS, VEGA};
use siderust::bodies::solar_system::Sun;
use siderust::bodies::Star;
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

fn main() {
    println!("╔════════════════════════════════════════════════════════╗");
    println!("║         Star Observability Planner                    ║");
    println!("╚════════════════════════════════════════════════════════╝\n");

    // Observatory: Greenwich
    let observatory = ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0));
    println!("Observatory: Greenwich Royal Observatory");
    println!("  Location: {}, {}\n", observatory.lat, observatory.lon);

    // Tonight: one night starting at MJD 60000.5 (local midnight approximation)
    let start = ModifiedJulianDate::new(60000.5);
    let end = ModifiedJulianDate::new(60001.5);
    let night = Period::new(start, end);

    println!("Observation window: {} to {}", start, end);
    println!("  (approximately tonight's darkness)\n");

    // Find astronomical night (Sun below -18°)
    let dark_periods = Sun.below_threshold(observatory, night, Degrees::new(-18.0));

    if dark_periods.is_empty() {
        println!("⚠ No astronomical darkness available tonight!");
        println!("  (Twilight all night or polar summer)\n");
        return;
    }

    let total_dark_hours = dark_periods.iter().fold(Hours::new(0.0), |acc, p| {
        acc + p.duration_days().to::<Hour>()
    });
    println!("✓ Astronomical night duration: {}\n", total_dark_hours);

    // Target stars
    let sirius = &SIRIUS;
    let vega = &VEGA;
    let altair = &ALTAIR;
    let betelgeuse = &BETELGEUSE;
    let rigel = &RIGEL;
    let polaris = &POLARIS;

    let targets: Vec<(&str, &Star)> = vec![
        ("Sirius", sirius),
        ("Vega", vega),
        ("Altair", altair),
        ("Betelgeuse", betelgeuse),
        ("Rigel", rigel),
        ("Polaris", polaris),
    ];

    println!("═══════════════════════════════════════════════════════");
    println!("  Target Visibility During Astronomical Night");
    println!("═══════════════════════════════════════════════════════\n");

    // Minimum altitude for good observation: 30° (to reduce atmospheric extinction)
    let min_altitude = Degrees::new(30.0);

    for &(name, star) in &targets {
        // Find when star is above minimum altitude
        let visible_periods = star.above_threshold(observatory, night, min_altitude);

        // Filter to only dark periods (intersection would be better, but this demonstrates the API)
        let observable_hours = visible_periods.iter().fold(Hours::new(0.0), |acc, p| {
            acc + p.duration_days().to::<Hour>()
        });

        print!("{:12} ", name);

        if visible_periods.is_empty() {
            println!("❌ Not observable (never above {})", min_altitude);
        } else if observable_hours < Hours::new(1.0) {
            println!(
                "⚠  Limited window: {} above {}",
                observable_hours.to::<Minute>(),
                min_altitude
            );
        } else {
            println!(
                "✓  Observable for {} above {}",
                observable_hours, min_altitude
            );

            // Show peak altitude during the night
            if let Some(period) = visible_periods.first() {
                let mid_mjd = siderust::time::ModifiedJulianDate::new(
                    (period.start.value() + period.end.value()) / 2.0,
                );
                let peak_alt = star.altitude_at(&observatory, mid_mjd).to::<Degree>();
                println!("             Peak altitude: {}", peak_alt);
            }
        }
    }

    println!("\n═══════════════════════════════════════════════════════");
    println!("  Observing Strategy Recommendations");
    println!("═══════════════════════════════════════════════════════\n");

    // Find the best 2-hour window for multiple targets
    println!("Best observing window: middle of astronomical night");

    if let Some(dark) = dark_periods.first() {
        let mid_mjd =
            siderust::time::ModifiedJulianDate::new((dark.start.value() + dark.end.value()) / 2.0);

        println!("  Around {}\n", mid_mjd);
        println!("Altitudes at this time:");

        for &(name, star) in &targets {
            let alt = star.altitude_at(&observatory, mid_mjd).to::<Degree>();
            let status = if alt > Degrees::new(30.0) {
                "✓ Good"
            } else if alt > Degrees::new(0.0) {
                "⚠ Low"
            } else {
                "❌ Below horizon"
            };
            println!("  {:12} {:>8}  {}", name, alt, status);
        }
    }

    println!("\n═══════════════════════════════════════════════════════");
    println!("Notes:");
    println!("  • Altitudes > 30° recommended for best image quality");
    println!("  • Consider weather, seeing, and Moon phase");
    println!("  • Circumpolar stars (e.g., Polaris) are always visible");
    println!("═══════════════════════════════════════════════════════");
}
