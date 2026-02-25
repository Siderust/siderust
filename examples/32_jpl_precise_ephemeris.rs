// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! JPL Precise Ephemeris Example
//!
//! This example demonstrates how to use JPL DE440/DE441 ephemeris backends for
//! high-precision calculations of solar system body positions and velocities,
//! and compares them with the always-available VSOP87/ELP2000 analytical series.
//!
//! ## Features Demonstrated
//! - Computing Sun, Earth, and Moon positions using JPL ephemerides
//! - Comparing DE4xx precision with VSOP87/ELP2000
//! - Working with different coordinate centers (barycentric, heliocentric, geocentric)
//! - Computing velocities for aberration and Doppler corrections
//! - Side-by-side backend comparison when both DE440 and DE441 are available
//!
//! ## Requirements
//! This example requires at least one JPL feature to be enabled:
//! ```bash
//! cargo run --example jpl_precise_ephemeris --features de440
//! cargo run --example jpl_precise_ephemeris --features de441
//! cargo run --example jpl_precise_ephemeris --features de440,de441
//! ```

// Gate on at least one JPL backend
#[cfg(not(any(feature = "de440", feature = "de441")))]
fn main() {
    eprintln!("ERROR: This example requires the 'de440' and/or 'de441' feature.");
    eprintln!("Run with: cargo run --example jpl_precise_ephemeris --features de440");
    eprintln!("     or:  cargo run --example jpl_precise_ephemeris --features de441");
    eprintln!("     or:  cargo run --example jpl_precise_ephemeris --features de440,de441");
    std::process::exit(1);
}

#[cfg(any(feature = "de440", feature = "de441"))]
fn main() {
    use chrono::{NaiveDate, NaiveTime, TimeZone, Utc};
    use qtty::*;
    use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
    use siderust::time::JulianDate;

    #[cfg(feature = "de440")]
    use siderust::calculus::ephemeris::De440Ephemeris;
    #[cfg(feature = "de441")]
    use siderust::calculus::ephemeris::De441Ephemeris;

    fn stub_enabled_for(prefix: &str) -> bool {
        let Ok(raw) = std::env::var("SIDERUST_JPL_STUB") else {
            return false;
        };
        let raw = raw.trim();
        if raw.is_empty() {
            return false;
        }

        let lower = raw.to_ascii_lowercase();
        if lower == "all" || lower == "1" || lower == "true" || lower == "yes" || lower == "on" {
            return true;
        }

        lower
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .any(|tok| tok == prefix)
    }

    /// Build a [`JulianDate`] from a calendar date (year, month, day at 12:00 TT).
    fn jd_from_ymd(year: i32, month: u32, day: u32) -> JulianDate {
        let naive = NaiveDate::from_ymd_opt(year, month, day)
            .expect("invalid date")
            .and_time(NaiveTime::from_hms_opt(12, 0, 0).unwrap());
        JulianDate::from_utc(Utc.from_utc_datetime(&naive))
    }

    println!("═══════════════════════════════════════════════════════════");
    println!("  JPL Precise Ephemeris — Backend Comparison");
    println!("═══════════════════════════════════════════════════════════\n");

    // =========================================================================
    // 1. Available Backends
    // =========================================================================
    println!("1. AVAILABLE BACKENDS");
    println!("─────────────────────");

    let stub_de440 = stub_enabled_for("de440");
    let stub_de441_env = stub_enabled_for("de441");
    let stub_de441_cfg = cfg!(siderust_mock_de441);
    let stub_de441 = stub_de441_env || stub_de441_cfg;

    println!("  VSOP87/ELP2000 : always available (analytical series)");
    #[cfg(feature = "de440")]
    println!(
        "  DE440          : enabled (JPL, 1550–2650 CE){}",
        if stub_de440 {
            " — STUBBED via SIDERUST_JPL_STUB (runtime calls skipped)"
        } else {
            ""
        }
    );
    #[cfg(not(feature = "de440"))]
    println!("  DE440          : not enabled");
    #[cfg(feature = "de441")]
    println!(
        "  DE441          : enabled (JPL, extended coverage){}",
        if stub_de441 {
            " — STUBBED/MOCKED (runtime calls compare against VSOP87/ELP2000)"
        } else {
            ""
        }
    );
    #[cfg(not(feature = "de441"))]
    println!("  DE441          : not enabled");
    println!();

    // =========================================================================
    // 2. Time Setup
    // =========================================================================
    println!("2. TIME SETUP");
    println!("─────────────");

    let j2000 = JulianDate::J2000;
    let test_date = jd_from_ymd(2026, 2, 15);
    println!("  J2000.0 epoch: JD {:.6}", j2000);
    println!("  Test date:     JD {:.6} (2026-02-15)\n", test_date);

    // =========================================================================
    // 3. Earth Heliocentric Position Comparison
    // =========================================================================
    println!("3. EARTH HELIOCENTRIC POSITION (at J2000)");
    println!("──────────────────────────────────────────");

    let earth_vsop = Vsop87Ephemeris::earth_heliocentric(j2000);
    println!("  VSOP87:");
    println!("    X = {:15.10} AU", earth_vsop.x());
    println!("    Y = {:15.10} AU", earth_vsop.y());
    println!("    Z = {:15.10} AU", earth_vsop.z());

    #[cfg(feature = "de440")]
    {
        if stub_de440 {
            println!("\n  DE440:");
            println!("    (stubbed) Skipping DE440 runtime calls to avoid panics.");
        } else {
            let earth_de440 = De440Ephemeris::earth_heliocentric(j2000);
            println!("\n  DE440:");
            println!("    X = {:15.10} AU", earth_de440.x());
            println!("    Y = {:15.10} AU", earth_de440.y());
            println!("    Z = {:15.10} AU", earth_de440.z());

            let diff_km = earth_vsop.distance_to(&earth_de440).to::<Kilometer>();
            println!("    Δ(VSOP87−DE440) = {}", diff_km);
        }
    }

    #[cfg(feature = "de441")]
    {
        if stub_de441 {
            println!("\n  DE441:");
            println!("    (stubbed/mocked) Skipping DE441 runtime calls.");
        } else {
            let earth_de441 = De441Ephemeris::earth_heliocentric(j2000);
            println!("\n  DE441:");
            println!("    X = {:15.10} AU", earth_de441.x());
            println!("    Y = {:15.10} AU", earth_de441.y());
            println!("    Z = {:15.10} AU", earth_de441.z());

            let diff_km = earth_vsop.distance_to(&earth_de441).to::<Kilometer>();
            println!("    Δ(VSOP87−DE441) = {}", diff_km);
        }
    }

    // When both are available, compare DE440 vs DE441
    #[cfg(all(feature = "de440", feature = "de441"))]
    {
        if stub_de440 || stub_de441 {
            println!("\n  Δ(DE440−DE441) = (skipped: stubbed/mocked backend enabled)");
        } else {
            let e440 = De440Ephemeris::earth_heliocentric(j2000);
            let e441 = De441Ephemeris::earth_heliocentric(j2000);
            let diff_km = e440.distance_to(&e441).to::<Kilometer>();
            println!("\n  Δ(DE440−DE441) = {}", diff_km);
        }
    }
    println!();

    // =========================================================================
    // 4. Moon Geocentric Position Comparison
    // =========================================================================
    println!("4. MOON GEOCENTRIC POSITION (at J2000)");
    println!("───────────────────────────────────────");

    let moon_vsop = Vsop87Ephemeris::moon_geocentric(j2000);
    println!("  ELP2000-82B:");
    println!("    X = {:15.6} km", moon_vsop.x());
    println!("    Y = {:15.6} km", moon_vsop.y());
    println!("    Z = {:15.6} km", moon_vsop.z());
    println!("    Distance: {:.6} km", moon_vsop.distance());

    #[cfg(feature = "de440")]
    {
        if stub_de440 {
            println!("\n  DE440:");
            println!("    (stubbed) Skipping DE440 runtime calls to avoid panics.");
        } else {
            let moon_de440 = De440Ephemeris::moon_geocentric(j2000);
            println!("\n  DE440:");
            println!("    X = {:15.6} km", moon_de440.x());
            println!("    Y = {:15.6} km", moon_de440.y());
            println!("    Z = {:15.6} km", moon_de440.z());
            println!("    Distance: {:.6} km", moon_de440.distance());

            let diff_km = moon_vsop.distance_to(&moon_de440);
            println!("    Δ(ELP2000−DE440) = {}", diff_km);
        }
    }

    #[cfg(feature = "de441")]
    {
        if stub_de441 {
            println!("\n  DE441:");
            println!("    (stubbed/mocked) Skipping DE441 runtime calls.");
        } else {
            let moon_de441 = De441Ephemeris::moon_geocentric(j2000);
            println!("\n  DE441:");
            println!("    X = {:15.6} km", moon_de441.x());
            println!("    Y = {:15.6} km", moon_de441.y());
            println!("    Z = {:15.6} km", moon_de441.z());
            println!("    Distance: {:.6} km", moon_de441.distance());

            let diff_km = moon_vsop.distance_to(&moon_de441);
            println!("    Δ(ELP2000−DE441) = {}", diff_km);
        }
    }
    println!();

    // =========================================================================
    // 5. Earth Velocity (for aberration corrections)
    // =========================================================================
    println!("5. EARTH VELOCITY (Barycentric, at test date)");
    println!("──────────────────────────────────────────────");

    type KmPerSecond = qtty::Per<Kilometer, Second>;

    let vel_vsop = Vsop87Ephemeris::earth_barycentric_velocity(test_date);
    let speed_vsop = vel_vsop.magnitude();
    println!(
        "  VSOP87:  {} = {}",
        speed_vsop,
        speed_vsop.to::<KmPerSecond>()
    );

    #[cfg(feature = "de440")]
    {
        if stub_de440 {
            println!("  DE440:   (stubbed) Skipping DE440 runtime calls to avoid panics.");
        } else {
            let vel = De440Ephemeris::earth_barycentric_velocity(test_date);
            let speed = vel.magnitude();
            println!("  DE440:   {} = {}", speed, speed.to::<KmPerSecond>());
        }
    }

    #[cfg(feature = "de441")]
    {
        if stub_de441 {
            println!("  DE441:   (stubbed/mocked) Skipping DE441 runtime calls.");
        } else {
            let vel = De441Ephemeris::earth_barycentric_velocity(test_date);
            let speed = vel.magnitude();
            println!("  DE441:   {} = {}", speed, speed.to::<KmPerSecond>());
        }
    }
    println!();

    // =========================================================================
    // 6. Time Series: Earth-Sun Distance Over a Year
    // =========================================================================
    println!("6. EARTH–SUN DISTANCE OVER 2026 (primary JPL backend)");
    println!("──────────────────────────────────────────────────────");
    println!("   Month    JD            Distance (AU)   Variation (%)");
    println!("  ────────────────────────────────────────────────────────");

    let mut min_dist: Option<AstronomicalUnits> = None;
    let mut max_dist: Option<AstronomicalUnits> = None;

    for month in 1..=12 {
        let jd = jd_from_ymd(2026, month, 15);

        // Use the highest-fidelity backend that is actually available at runtime.
        #[cfg(feature = "de441")]
        let dist = if !stub_de441 {
            De441Ephemeris::earth_heliocentric(jd).distance()
        } else {
            Vsop87Ephemeris::earth_heliocentric(jd).distance()
        };
        #[cfg(all(feature = "de440", not(feature = "de441")))]
        let dist = if !stub_de440 {
            De440Ephemeris::earth_heliocentric(jd).distance()
        } else {
            Vsop87Ephemeris::earth_heliocentric(jd).distance()
        };

        if min_dist.map_or(true, |min| dist < min) {
            min_dist = Some(dist);
        }
        if max_dist.map_or(true, |max| dist > max) {
            max_dist = Some(dist);
        }

        let variation =
            ((dist - AstronomicalUnits::new(1.0)) / AstronomicalUnits::new(1.0)) * 100.0;
        println!(
            "   {:2}       {:.2}    {}   {:+6.3}%",
            month, jd, dist, variation
        );
    }

    let min_dist = min_dist.expect("time series has at least one sample");
    let max_dist = max_dist.expect("time series has at least one sample");
    let range = max_dist - min_dist;
    let mean = (max_dist + min_dist) / 2.0;

    println!("\n  Perihelion: {}", min_dist);
    println!("  Aphelion:   {}", max_dist);
    println!(
        "  Eccentricity effect: {} ({:.3}%)\n",
        range,
        (range / mean) * 100.0
    );

    // =========================================================================
    // Summary
    // =========================================================================
    println!("═══════════════════════════════════════════════════════════");
    println!("  Summary");
    println!("═══════════════════════════════════════════════════════════");
    println!();
    println!("  DE4xx backends provide higher-fidelity solar-system ephemerides");
    println!("  than the analytical VSOP87/ELP2000 series (at the cost of larger datasets).");
    println!();
    println!("  Backend selection:");
    println!("  • DE440: 1550–2650 CE (modern focus)");
    println!("  • DE441: Extended coverage (deep past/future)");
    println!("  • VSOP87/ELP2000: Always available analytical series");
    println!();
    println!("  DefaultEphemeris auto-selects: DE441 > DE440 > VSOP87");
    println!();
}
