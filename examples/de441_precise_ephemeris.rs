// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE441 Precise Ephemeris Example
//!
//! This example demonstrates how to use JPL DE441 ephemeris for high-precision
//! calculations of solar system body positions and velocities. DE441 provides
//! significantly better accuracy than analytical theories like VSOP87/ELP2000.
//!
//! ## Features Demonstrated
//! - Computing Sun, Earth, and Moon positions using DE441
//! - Comparing DE441 with VSOP87 to show precision improvement
//! - Working with different coordinate centers (barycentric, heliocentric, geocentric)
//! - Computing velocities for aberration and Doppler corrections
//! - Distance calculations and relative positions
//!
//! ## Requirements
//! This example requires the `de441` feature to be enabled:
//! ```bash
//! cargo run --example de441_precise_ephemeris --features de441
//! ```

#[cfg(not(feature = "de441"))]
fn main() {
    eprintln!("ERROR: This example requires the 'de441' feature.");
    eprintln!("Run with: cargo run --example de441_precise_ephemeris --features de441");
    std::process::exit(1);
}

#[cfg(feature = "de441")]
fn main() {
    use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
    use qtty::*;
    use siderust::calculus::ephemeris::{De441Ephemeris, Ephemeris, Vsop87Ephemeris};
    use siderust::time::JulianDate;

    /// Build a [`JulianDate`] from a calendar date (year, month, day at 12:00 TT).
    fn jd_from_ymd(year: i32, month: u32, day: u32) -> JulianDate {
        let naive = NaiveDate::from_ymd_opt(year, month, day)
            .expect("invalid date")
            .and_time(NaiveTime::from_hms_opt(12, 0, 0).unwrap());
        JulianDate::from_utc(Utc.from_utc_datetime(&naive))
    }

    println!("═══════════════════════════════════════════════════════════");
    println!("  DE441 Precise Ephemeris Example");
    println!("═══════════════════════════════════════════════════════════\n");

    // =========================================================================
    // 1. Time Setup
    // =========================================================================
    println!("1. TIME SETUP");
    println!("─────────────");

    let j2000 = JulianDate::J2000;
    let current = JulianDate::from_utc(chrono::Utc::now());
    let test_date = jd_from_ymd(2026, 2, 15);

    println!("  J2000.0 epoch: JD {:.6}", j2000);
    println!("  Current time:  JD {:.6}", current);
    println!("  Test date:     JD {:.6} (2026-02-15)\n", test_date);

    // =========================================================================
    // 2. Sun Position (Barycentric)
    // =========================================================================
    println!("2. SUN POSITION (Barycentric Ecliptic)");
    println!("───────────────────────────────────────");

    let sun_de441 = De441Ephemeris::sun_barycentric(test_date);
    println!("  Using DE441:");
    println!("    X = {:15.10} AU", sun_de441.get_position().x());
    println!("    Y = {:15.10} AU", sun_de441.get_position().y());
    println!("    Z = {:15.10} AU", sun_de441.get_position().z());
    println!(
        "    Distance from SSB: {:.10} AU\n",
        sun_de441.get_position().distance()
    );

    // =========================================================================
    // 3. Earth Position (Multiple Centers)
    // =========================================================================
    println!("3. EARTH POSITION (Multiple Reference Frames)");
    println!("──────────────────────────────────────────────");

    // Earth barycentric
    let earth_bary_de441 = De441Ephemeris::earth_barycentric(test_date);
    println!("  Earth Barycentric (DE441):");
    println!("    X = {:15.10} AU", earth_bary_de441.get_position().x());
    println!("    Y = {:15.10} AU", earth_bary_de441.get_position().y());
    println!("    Z = {:15.10} AU", earth_bary_de441.get_position().z());
    println!(
        "    Distance from SSB: {:.10} AU",
        earth_bary_de441.get_position().distance()
    );

    // Earth heliocentric
    let earth_helio_de441 = De441Ephemeris::earth_heliocentric(test_date);
    println!("\n  Earth Heliocentric (DE441):");
    println!("    X = {:15.10} AU", earth_helio_de441.get_position().x());
    println!("    Y = {:15.10} AU", earth_helio_de441.get_position().y());
    println!("    Z = {:15.10} AU", earth_helio_de441.get_position().z());
    println!(
        "    Distance from Sun: {:.10} AU (should be ~1.0)\n",
        earth_helio_de441.get_position().distance()
    );

    // =========================================================================
    // 4. Moon Position (Geocentric)
    // =========================================================================
    println!("4. MOON POSITION (Geocentric Ecliptic)");
    println!("───────────────────────────────────────");

    let moon_de441 = De441Ephemeris::moon_geocentric(test_date);
    println!("  Using DE441:");
    println!("    X = {:15.6} km", moon_de441.x());
    println!("    Y = {:15.6} km", moon_de441.y());
    println!("    Z = {:15.6} km", moon_de441.z());
    println!("    Distance from Earth: {:.6} km", moon_de441.distance());

    // Convert to Earth radii for perspective
    const EARTH_RADIUS_KM: f64 = 6371.0;
    println!(
        "    Distance: {:.2} Earth radii\n",
        moon_de441.distance().value() / EARTH_RADIUS_KM
    );

    // =========================================================================
    // 5. Earth Velocity (for aberration corrections)
    // =========================================================================
    println!("5. EARTH VELOCITY (Barycentric Ecliptic)");
    println!("─────────────────────────────────────────");

    let earth_vel = De441Ephemeris::earth_barycentric_velocity(test_date);
    println!("  Using DE441:");
    println!("    Vx = {:15.10} AU/day", earth_vel.x());
    println!("    Vy = {:15.10} AU/day", earth_vel.y());
    println!("    Vz = {:15.10} AU/day", earth_vel.z());

    // Compute speed in km/s
    const AU_TO_KM: f64 = 149_597_870.7;
    const DAY_TO_SECONDS: f64 = 86_400.0;
    let speed_au_day = (earth_vel.x().value().powi(2)
        + earth_vel.y().value().powi(2)
        + earth_vel.z().value().powi(2))
    .sqrt();
    let speed_km_s = speed_au_day * AU_TO_KM / DAY_TO_SECONDS;
    println!(
        "    Speed: {:.6} AU/day = {:.3} km/s\n",
        speed_au_day, speed_km_s
    );

    // =========================================================================
    // 6. Precision Comparison: DE441 vs VSOP87
    // =========================================================================
    println!("6. PRECISION COMPARISON: DE441 vs VSOP87");
    println!("────────────────────────────────────────");

    println!("  Comparing Earth heliocentric position at J2000:");

    let earth_vsop = Vsop87Ephemeris::earth_heliocentric(j2000);
    let earth_de441_j2000 = De441Ephemeris::earth_heliocentric(j2000);

    println!("\n  VSOP87:");
    println!("    X = {:15.10} AU", earth_vsop.get_position().x());
    println!("    Y = {:15.10} AU", earth_vsop.get_position().y());
    println!("    Z = {:15.10} AU", earth_vsop.get_position().z());

    println!("\n  DE441:");
    println!("    X = {:15.10} AU", earth_de441_j2000.get_position().x());
    println!("    Y = {:15.10} AU", earth_de441_j2000.get_position().y());
    println!("    Z = {:15.10} AU", earth_de441_j2000.get_position().z());

    // Calculate differences
    let dx = (earth_vsop.get_position().x() - earth_de441_j2000.get_position().x()).value();
    let dy = (earth_vsop.get_position().y() - earth_de441_j2000.get_position().y()).value();
    let dz = (earth_vsop.get_position().z() - earth_de441_j2000.get_position().z()).value();
    let diff_au = (dx * dx + dy * dy + dz * dz).sqrt();
    let diff_km = diff_au * AU_TO_KM;

    println!("\n  Position Difference:");
    println!("    ΔX = {:15.10} AU ({:12.3} km)", dx, dx * AU_TO_KM);
    println!("    ΔY = {:15.10} AU ({:12.3} km)", dy, dy * AU_TO_KM);
    println!("    ΔZ = {:15.10} AU ({:12.3} km)", dz, dz * AU_TO_KM);
    println!("    Total: {:.6} AU = {:.3} km", diff_au, diff_km);
    println!("    (DE441 is significantly more accurate)\n");

    // =========================================================================
    // 7. Moon Precision Comparison
    // =========================================================================
    println!("7. MOON PRECISION COMPARISON: DE441 vs ELP2000");
    println!("───────────────────────────────────────────────");

    let moon_vsop = Vsop87Ephemeris::moon_geocentric(j2000);
    let moon_de441_j2000 = De441Ephemeris::moon_geocentric(j2000);

    println!("  Comparing Moon geocentric position at J2000:");

    println!("\n  ELP2000-82B:");
    println!("    X = {:15.6} km", moon_vsop.x());
    println!("    Y = {:15.6} km", moon_vsop.y());
    println!("    Z = {:15.6} km", moon_vsop.z());
    println!("    Distance: {:.6} km", moon_vsop.distance());

    println!("\n  DE441:");
    println!("    X = {:15.6} km", moon_de441_j2000.x());
    println!("    Y = {:15.6} km", moon_de441_j2000.y());
    println!("    Z = {:15.6} km", moon_de441_j2000.z());
    println!("    Distance: {:.6} km", moon_de441_j2000.distance());

    let moon_dx = (moon_vsop.x() - moon_de441_j2000.x()).value();
    let moon_dy = (moon_vsop.y() - moon_de441_j2000.y()).value();
    let moon_dz = (moon_vsop.z() - moon_de441_j2000.z()).value();
    let moon_diff_km = (moon_dx * moon_dx + moon_dy * moon_dy + moon_dz * moon_dz).sqrt();

    println!("\n  Position Difference:");
    println!("    ΔX = {:15.6} km", moon_dx);
    println!("    ΔY = {:15.6} km", moon_dy);
    println!("    ΔZ = {:15.6} km", moon_dz);
    println!("    Total: {:.6} km", moon_diff_km);
    println!("    (DE441 provides meter-level precision)\n");

    // =========================================================================
    // 8. Time Series: Earth-Sun Distance Over a Year
    // =========================================================================
    println!("8. TIME SERIES: Earth-Sun Distance Over a Year");
    println!("───────────────────────────────────────────────");

    println!("  Computing Earth-Sun distance for 12 months (2026):\n");
    println!("   Month    JD            Distance (AU)   Variation (%)");
    println!("  ────────────────────────────────────────────────────────");

    let mut min_dist = f64::MAX;
    let mut max_dist = f64::MIN;

    for month in 1..=12 {
        let jd = jd_from_ymd(2026, month, 15);
        let earth = De441Ephemeris::earth_heliocentric(jd);
        let dist = earth.get_position().distance().value();

        min_dist = min_dist.min(dist);
        max_dist = max_dist.max(dist);

        let variation = (dist - 1.0) / 1.0 * 100.0;
        println!(
            "   {:2}       {:.2}    {:.10}   {:+6.3}%",
            month, jd, dist, variation
        );
    }

    println!("\n  Perihelion: {:.10} AU", min_dist);
    println!("  Aphelion:   {:.10} AU", max_dist);
    println!(
        "  Eccentricity effect: {:.6} AU ({:.3}%)\n",
        max_dist - min_dist,
        (max_dist - min_dist) / ((max_dist + min_dist) / 2.0) * 100.0
    );

    // =========================================================================
    // Summary
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════");
    println!("  Summary");
    println!("═══════════════════════════════════════════════════════════");
    println!();
    println!("  DE441 provides:");
    println!("  • Kilometer-level precision for planetary positions");
    println!("  • Meter-level precision for lunar positions");
    println!("  • Coverage: DE441 part-2 kernel (modern/future epochs)");
    println!("  • Based on JPL's latest planetary ephemeris");
    println!();
    println!("  Use DE441 for:");
    println!("  • Precise aberration and light-time corrections");
    println!("  • High-accuracy satellite tracking");
    println!("  • Scientific observations requiring meter-level precision");
    println!("  • Modern epoch ephemeris (near current date)");
    println!();
    println!("  Use VSOP87/ELP2000 when:");
    println!("  • Sub-kilometer precision is sufficient");
    println!("  • Historical dates outside DE441 coverage");
    println!("  • Smaller binary size is important");
    println!();
}
