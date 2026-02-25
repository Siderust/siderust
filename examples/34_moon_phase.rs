// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Moon Phase Example
//!
//! Demonstrates:
//!
//! 1. **Point-in-time phase** — illuminated fraction, phase angle, waxing/waning,
//!    and eight-label classification at a given UTC instant and observer location.
//!
//! 2. **Phase event finder** — locating the next New Moon, First Quarter,
//!    Full Moon, and Last Quarter in a 35-day window.
//!
//! 3. **Illumination period finder** — finding time windows when the Moon's
//!    illuminated fraction is within a specified range (analogous to the
//!    `altitude_ranges` / `above_threshold` API for solar/stellar bodies).
//!
//! ## Usage
//!
//! ```bash
//! cargo run --example moon_phase
//! cargo run --example moon_phase -- 2026-03-01 28.76 -17.89 2396
//! ```
//!
//! ## Arguments (all optional, positional)
//!
//! | # | Description             | Default                       |
//! |---|-------------------------|-------------------------------|
//! | 1 | Start date `YYYY-MM-DD` | today UTC                     |
//! | 2 | Observer latitude (°N)  | 28.762 (Roque de los Muchachos)|
//! | 3 | Observer longitude (°E) | −17.892                       |
//! | 4 | Observer height (m)     | 2396                          |

use chrono::{NaiveDate, NaiveDateTime, NaiveTime, TimeZone, Utc};
use qtty::*;
use siderust::calculus::ephemeris::Vsop87Ephemeris;
use siderust::calculus::lunar::phase::{
    find_phase_events, illumination_range, moon_phase_geocentric, moon_phase_topocentric,
    PhaseSearchOpts,
};
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use siderust::time::{JulianDate, ModifiedJulianDate, Period};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // -----------------------------------------------------------------------
    // Parse arguments
    // -----------------------------------------------------------------------
    let today = Utc::now().date_naive();
    let start_date = args
        .get(1)
        .map(|s| NaiveDate::parse_from_str(s, "%Y-%m-%d").expect("Invalid date (use YYYY-MM-DD)"))
        .unwrap_or(today);

    let lat: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(28.762);
    let lon: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(-17.892);
    let h_m: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(2396.0);

    let site = Geodetic::<ECEF>::new(
        Degrees::new(lon),
        Degrees::new(lat),
        Quantity::<Meter>::new(h_m),
    );

    // Convert start date (midnight UTC) → JulianDate / MJD
    let midnight = NaiveDateTime::new(start_date, NaiveTime::from_hms_opt(0, 0, 0).unwrap());
    let jd_start = JulianDate::from_utc(Utc.from_utc_datetime(&midnight));
    let mjd_start = ModifiedJulianDate::from(jd_start);
    let window = Period::new(mjd_start, mjd_start + Days::new(35.0));

    // -----------------------------------------------------------------------
    // 1. Point-in-time phase (geocentric + topocentric) at midnight of start
    // -----------------------------------------------------------------------
    let geo = moon_phase_geocentric::<Vsop87Ephemeris>(jd_start);
    let topo = moon_phase_topocentric::<Vsop87Ephemeris>(jd_start, site);

    println!("═══════════════════════════════════════════════════════════");
    println!(" Moon Phase — {start_date} 00:00 UTC");
    println!("═══════════════════════════════════════════════════════════");
    println!(" Observer : lat {lat:.3}°, lon {lon:.3}°, height {h_m:.0} m");
    println!();
    println!(" Geocentric:");
    println!("   Phase label        : {}", geo.label());
    println!("   Illuminated        : {:.1} %", geo.illuminated_percent());
    println!("   Phase angle        : {}", geo.phase_angle.to::<Degree>());
    println!("   Elongation         : {}", geo.elongation.to::<Degree>());
    println!(
        "   Waxing             : {}",
        if geo.waxing { "yes" } else { "no" }
    );
    println!();
    println!(" Topocentric (parallax corrected):");
    println!("   Phase label        : {}", topo.label());
    println!(
        "   Illuminated        : {:.1} %",
        topo.illuminated_percent()
    );
    println!("   Elongation         : {}", topo.elongation.to::<Degree>());
    println!(
        "   Illumination diff  : {:+.3} %",
        (topo.illuminated_fraction - geo.illuminated_fraction) * 100.0
    );

    // -----------------------------------------------------------------------
    // 2. Principal phase events in the next 35 days
    // -----------------------------------------------------------------------
    let events = find_phase_events::<Vsop87Ephemeris>(window, PhaseSearchOpts::default());

    println!();
    println!("───────────────────────────────────────────────────────────");
    println!(" Principal phase events — next 35 days");
    println!("───────────────────────────────────────────────────────────");
    if events.is_empty() {
        println!("  (none found in window)");
    } else {
        for ev in &events {
            if let Some(utc) = ev.mjd.to_utc() {
                println!(
                    "  {:15}  {}",
                    ev.kind.to_string(),
                    utc.format("%Y-%m-%d %H:%M UTC")
                );
            } else {
                println!("  {:15}  MJD {}", ev.kind.to_string(), ev.mjd);
            }
        }
    }

    // -----------------------------------------------------------------------
    // 3. Illumination range periods — crescent and gibbous windows
    // -----------------------------------------------------------------------
    let opts = PhaseSearchOpts::default();

    // "Dark sky" window: Moon less than 25% illuminated
    let dark_windows = illumination_range::<Vsop87Ephemeris>(window, 0.0, 0.25, opts);

    // "Bright Moon" window: more than 75% illuminated
    let bright_windows = illumination_range::<Vsop87Ephemeris>(window, 0.75, 1.0, opts);

    // "Crescent" window: 5–40% illuminated
    let crescent_windows = illumination_range::<Vsop87Ephemeris>(window, 0.05, 0.40, opts);

    println!();
    println!("───────────────────────────────────────────────────────────");
    println!(" Illumination range periods (geocentric) — next 35 days");
    println!("───────────────────────────────────────────────────────────");

    print_periods("Dark sky (0–25%)", &dark_windows);
    print_periods("Crescent  (5–40%)", &crescent_windows);
    print_periods("Bright Moon (75–100%)", &bright_windows);
}

// ---------------------------------------------------------------------------
// Formatting helper
// ---------------------------------------------------------------------------

fn print_periods(label: &str, periods: &[siderust::time::Period<siderust::time::MJD>]) {
    println!();
    println!("  {label}:");
    if periods.is_empty() {
        println!("    (none in window)");
        return;
    }
    let total_days: Days = periods
        .iter()
        .map(|p| p.end - p.start)
        .fold(Days::new(0.0), |acc, d| acc + d);
    for p in periods {
        let dur = p.end - p.start;
        match (p.start.to_utc(), p.end.to_utc()) {
            (Some(s), Some(e)) => println!(
                "    {} → {}  ({})",
                s.format("%m-%d %H:%M"),
                e.format("%m-%d %H:%M"),
                dur
            ),
            _ => println!("    MJD {} → {}  ({})", p.start, p.end, dur),
        }
    }
    println!("    Total : {}", total_days);
}
