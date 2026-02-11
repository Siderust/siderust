// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Periods — Trait-Based API Example
//!
//! Demonstrates the unified `AltitudePeriodsProvider` trait for finding time
//! intervals when celestial bodies are within specific altitude ranges.
//!
//! Run with: `cargo run --example altitude_periods_trait`

use siderust::bodies::catalog::{POLARIS, SIRIUS, VEGA};
use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::{AltitudePeriodsProvider, AltitudeQuery};
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::spherical::direction;
use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::time::{ModifiedJulianDate, Period, MJD};

use qtty::*;

fn main() {
    println!("=== Altitude Periods API Examples ===\n");

    // Observer: Roque de los Muchachos Observatory (La Palma)
    let observer = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
    println!("Observer: Roque de los Muchachos Observatory");
    println!(
        "  Location: {:.3}°N, {:.3}°E, {} m\n",
        observer.lat.value(),
        observer.lon.value(),
        observer.height.value()
    );

    // Time window: one week starting from MJD 60000
    let start = ModifiedJulianDate::new(60000.0);
    let end = ModifiedJulianDate::new(60007.0);
    let window = Period::new(start, end);
    println!(
        "Time window: MJD {:.1} to {:.1} (7 days)\n",
        start.value(),
        end.value()
    );

    // -----------------------------------------------------------------------
    // Example 1: Find astronomical night periods (Sun below -18°)
    // -----------------------------------------------------------------------
    println!("--- Example 1: Astronomical Nights ---");
    let astro_nights = Sun.below_threshold(observer, window, Degrees::new(-18.0));

    println!("Found {} astronomical night periods:", astro_nights.len());
    for (i, period) in astro_nights.iter().enumerate().take(3) {
        let hours = period.duration_days() * 24.0;
        println!(
            "  Night {}: {:.2} hours (JD {:.4} to {:.4})",
            i + 1,
            hours,
            period.start.value(),
            period.end.value()
        );
    }
    println!();

    // -----------------------------------------------------------------------
    // Example 2: Find when Sirius is above 30° altitude
    // -----------------------------------------------------------------------
    println!("--- Example 2: Sirius High Above Horizon ---");
    let sirius_high = SIRIUS.above_threshold(observer, window, Degrees::new(30.0));

    println!("Sirius above 30° altitude:");
    println!("  Found {} periods", sirius_high.len());
    if let Some(first) = sirius_high.first() {
        let hours = first.duration_days() * 24.0;
        println!("  First period: {:.2} hours", hours);
    }
    println!();

    // -----------------------------------------------------------------------
    // Example 3: Use direction::ICRS for a custom target
    // -----------------------------------------------------------------------
    println!("--- Example 3: Custom ICRS Direction (Betelgeuse) ---");
    // Betelgeuse: RA ≈ 88.79°, Dec ≈ +7.41°
    let betelgeuse = direction::ICRS::new(Degrees::new(88.79), Degrees::new(7.41));

    let betelgeuse_visible = betelgeuse.above_threshold(
        observer,
        window,
        Degrees::new(0.0), // above horizon
    );

    println!("Betelgeuse above horizon:");
    println!("  Found {} periods in 7 days", betelgeuse_visible.len());
    let total_hours: f64 = betelgeuse_visible
        .iter()
        .map(|p| p.duration_days().value() * 24.0)
        .sum();
    println!("  Total visible time: {:.1} hours", total_hours);
    println!();

    // -----------------------------------------------------------------------
    // Example 4: Moon altitude range query (low moon: 0° to 20°)
    // -----------------------------------------------------------------------
    println!("--- Example 4: Low Moon Periods (0° to 20°) ---");
    let query = AltitudeQuery {
        observer,
        window,
        min_altitude: Degrees::new(0.0),
        max_altitude: Degrees::new(20.0),
    };
    let low_moon = Moon.altitude_periods(&query);

    println!("Moon between 0° and 20° altitude:");
    println!("  Found {} periods", low_moon.len());
    for (i, period) in low_moon.iter().enumerate().take(2) {
        let minutes = period.duration_days() * 1440.0;
        println!("  Period {}: {:.1} minutes", i + 1, minutes);
    }
    println!();

    // -----------------------------------------------------------------------
    // Example 5: Circumpolar star check (Polaris)
    // -----------------------------------------------------------------------
    println!("--- Example 5: Circumpolar Star (Polaris) ---");
    let polaris_up = POLARIS.above_threshold(observer, window, Degrees::new(0.0));

    println!("Polaris above horizon at {:.1}°N:", observer.lat.value());
    if polaris_up.len() == 1 && (polaris_up[0].duration_days() - Days::new(7.0)).abs() < 0.1 {
        println!("  ✓ Circumpolar (continuously visible for entire week)");
    } else {
        println!(
            "  Found {} periods (total {:.2} days)",
            polaris_up.len(),
            polaris_up
                .iter()
                .map(|p| p.duration_days().value())
                .sum::<f64>()
        );
    }
    println!();

    // -----------------------------------------------------------------------
    // Example 6: Twilight band (nautical to astronomical)
    // -----------------------------------------------------------------------
    println!("--- Example 6: Twilight Band (−18° to −12°) ---");
    let twilight_query = AltitudeQuery {
        observer,
        window: Period::new(start, ModifiedJulianDate::new(60002.0)), // 2 days
        min_altitude: Degrees::new(-18.0),
        max_altitude: Degrees::new(-12.0),
    };
    let twilight = Sun.altitude_periods(&twilight_query);

    println!("Sun in nautical-to-astronomical twilight band:");
    println!("  Found {} twilight periods in 2 days", twilight.len());
    for (i, period) in twilight.iter().enumerate() {
        let minutes = period.duration_days() * 1440.0;
        println!("  Period {}: {:.1} minutes", i + 1, minutes);
    }
    println!();

    // -----------------------------------------------------------------------
    // Example 7: Point altitude evaluation
    // -----------------------------------------------------------------------
    println!("--- Example 7: Single-Point Altitude Queries ---");

    let sun_alt = Sun.altitude_at(&observer, start).to::<Degree>();
    let moon_alt = Moon.altitude_at(&observer, start).to::<Degree>();
    let vega_alt = VEGA.altitude_at(&observer, start).to::<Degree>();

    println!("Altitudes at MJD {:.1}:", start.value());
    println!("  Sun:  {:.2}°", sun_alt.value());
    println!("  Moon: {:.2}°", moon_alt.value());
    println!("  Vega: {:.2}°", vega_alt.value());
    println!();

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    println!("=== Summary ===");
    println!("The AltitudePeriodsProvider trait provides a unified interface for:");
    println!("  • Sun, Moon, and any Star from the catalog");
    println!("  • Lightweight direction::ICRS for custom RA/Dec coordinates");
    println!("  • Consistent API: above_threshold, below_threshold, altitude_periods");
    println!("  • Single-point queries: altitude_at(observer, mjd)");
}

fn _print_period_details(periods: &[Period<MJD>], label: &str) {
    println!("{}: {} periods", label, periods.len());
    for (i, p) in periods.iter().enumerate() {
        println!(
            "  {}: MJD {:.4} → {:.4} ({:.2} hours)",
            i + 1,
            p.start.value(),
            p.end.value(),
            p.duration_days() * 24.0
        );
    }
}
