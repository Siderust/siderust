// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Trackable Trait Demo
//!
//! Demonstrates the `Trackable` trait — a unified interface for any
//! astronomical object that can produce coordinates at a given time.
//!
//! Shows how planets, stars, fixed directions, and precomputed coordinate
//! samples can all be used interchangeably through trait-based generics.

use siderust::bodies::{catalog, solar_system};
use siderust::coordinates::spherical::direction;
use siderust::targets::Trackable;
use siderust::time::JulianDate;

fn main() {
    let jd = JulianDate::J2000;
    println!("=== Siderust Trackable Trait Demo ===");
    println!("Epoch: J2000.0 (JD {})\n", jd);

    // ─── 1. Planet (VSOP87 ephemeris) ───────────────────────────────────
    let earth_pos = solar_system::Earth.track(jd);
    println!(
        "Earth (barycentric ecliptic):\n  x = {}, y = {}, z = {}\n  distance from SSB = {}\n",
        earth_pos.position.x(),
        earth_pos.position.y(),
        earth_pos.position.z(),
        earth_pos.position.distance(),
    );

    // ─── 2. Moon (ELP2000) ──────────────────────────────────────────────
    let moon_pos = solar_system::Moon.track(jd);
    println!(
        "Moon (geocentric ecliptic):\n  distance from Earth center = {}\n",
        moon_pos.distance(),
    );

    // ─── 3. Star (catalog direction) ────────────────────────────────────
    let sirius = &catalog::SIRIUS;
    let dir = sirius.track(jd);
    println!(
        "Sirius (ICRS direction):\n  RA = {}, Dec = {}\n",
        dir.ra(),
        dir.dec(),
    );

    // ─── 4. Fixed ICRS direction (time-invariant) ───────────────────────
    let zenith_dir = direction::ICRS::new(qtty::Degrees::new(180.0), qtty::Degrees::new(45.0));
    let same = zenith_dir.track(jd);
    println!(
        "Fixed ICRS direction:\n  RA = {}, Dec = {} (unchanged at any epoch)\n",
        same.ra(),
        same.dec(),
    );

    // ─── 5. Generic function over Trackable ─────────────────────────────
    println!("--- Generic dispatch demo ---");
    print_tracking_info("Mars", &solar_system::Mars, jd);
    print_tracking_info("Jupiter", &solar_system::Jupiter, jd);
    print_tracking_info("Saturn", &solar_system::Saturn, jd);
}

/// Example of a generic function constrained by `Trackable`.
fn print_tracking_info<T: Trackable>(name: &str, obj: &T, jd: JulianDate)
where
    T::Coords: std::fmt::Debug,
{
    let coords = obj.track(jd);
    println!("  {name}: {coords:?}");
}
