// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Night Quality Scoring
//!
//! Practical example: Score each night in a month based on Moon interference
//! and astronomical darkness duration. Uses the trait API to compute both
//! Sun and Moon altitude periods.
//!
//! Run with: `cargo run --example night_quality_scoring`

use siderust::bodies::solar_system::{Moon, Sun};
use siderust::calculus::altitude::AltitudePeriodsProvider;
use siderust::coordinates::centers::ObserverSite;
use siderust::time::{ModifiedJulianDate, Period};

use qtty::*;

#[derive(Clone)]
struct NightScore {
    date: ModifiedJulianDate,
    dark_hours: f64,
    moon_up_hours: f64,
    score: f64,
}

fn score_night(night_start: ModifiedJulianDate, observer: ObserverSite) -> NightScore {
    let night_end = ModifiedJulianDate::new(night_start.value() + 1.0);
    let window = Period::new(night_start, night_end);

    // Astronomical darkness (Sun below -18°)
    let dark_periods = Sun.below_threshold(observer, window, Degrees::new(-18.0));
    let dark_hours: f64 = dark_periods
        .iter()
        .map(|p| p.duration_days().to::<Hour>())
        .sum();

    // Moon above horizon during dark time
    let moon_up = Moon.above_threshold(observer, window, Degrees::new(0.0));
    let moon_hours: f64 = moon_up
        .iter()
        .map(|p| p.duration_days().to::<Hour>())
        .sum();

    // Simple scoring: darkness good, Moon bad
    // Score = dark_hours * (1 - moon_interference_factor)
    let moon_interference = (moon_hours / 24.0).min(1.0);
    let score = dark_hours * (1.0 - 0.7 * moon_interference);

    NightScore {
        date: night_start,
        dark_hours,
        moon_up_hours: moon_hours,
        score,
    }
}

fn main() {
    println!("╔════════════════════════════════════════════════════════╗");
    println!("║         Monthly Observing Conditions Report           ║");
    println!("╚════════════════════════════════════════════════════════╝\n");

    // Mauna Kea Observatory, Hawaii
    let observatory = ObserverSite::new(
        Degrees::new(-155.472),
        Degrees::new(19.826),
        Meters::new(4207.0),
    );
    println!("Observatory: Mauna Kea, Hawaii");
    println!(
        "  Location: {}, {}, {} elevation\n",
        observatory.lat, observatory.lon, observatory.height
    );

    // Score 30 consecutive nights
    let start_mjd = 60000.0;
    let mut scores: Vec<NightScore> = Vec::new();

    println!("Analyzing 30 nights starting MJD {:.0}...\n", start_mjd);

    for day in 0..30 {
        let night_mjd = ModifiedJulianDate::new(start_mjd + day as f64);
        scores.push(score_night(night_mjd, observatory));
    }

    // Display results
    println!("═══════════════════════════════════════════════════════════════");
    println!(" Night │ MJD      │ Dark   │ Moon   │ Score │ Quality");
    println!("       │          │ (hrs)  │ (hrs)  │       │");
    println!("───────┼──────────┼────────┼────────┼───────┼─────────────");

    for (i, score) in scores.iter().enumerate() {
        let quality = if score.score > 8.0 {
            "★★★★★ Excellent"
        } else if score.score > 6.0 {
            "★★★★  Very Good"
        } else if score.score > 4.0 {
            "★★★   Good"
        } else if score.score > 2.0 {
            "★★    Fair"
        } else {
            "★     Poor"
        };

        println!(
            " {:5} │ {:>8} │ {:6.2} │ {:6.2} │ {:5.2} │ {}",
            i + 1,
            score.date,
            score.dark_hours,
            score.moon_up_hours,
            score.score,
            quality
        );
    }

    println!("═══════════════════════════════════════════════════════════════");

    // Find best nights
    let mut sorted_scores = scores.clone();
    sorted_scores.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    println!("\n🌟 Top 5 Best Nights for Observing:\n");
    for (i, score) in sorted_scores.iter().take(5).enumerate() {
        let night_num = ((score.date.value() - start_mjd) as usize) + 1;
        println!(
            "  {}. Night {} ({}): {:.2} hours dark, {:.2} hours Moon",
            i + 1,
            night_num,
            score.date,
            score.dark_hours,
            score.moon_up_hours
        );
    }

    println!("\n═══════════════════════════════════════════════════════════════");
    println!("Scoring formula:");
    println!("  Score = dark_hours × (1 - 0.7 × moon_interference)");
    println!("  where moon_interference = min(moon_up_hours / 24, 1.0)");
    println!();
    println!("This example demonstrates using the trait API to build");
    println!("practical observing planning tools with minimal code.");
    println!("═══════════════════════════════════════════════════════════════");
}
