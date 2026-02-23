// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Targets + Proper Motion Example
//!
//! Demonstrates:
//! - `CoordinateWithPM<T>` as a (position, epoch, proper-motion) bundle
//! - `ProperMotion` construction for Gaia/Hipparcos-style `µα⋆`
//! - Applying proper motion from J2000.0 to a future epoch

use qtty::*;
use siderust::astro::proper_motion::{set_proper_motion_since_j2000, ProperMotion};
use siderust::coordinates::spherical::position;
use siderust::targets::CoordinateWithPM;
use siderust::time::JulianDate;

fn main() {
    println!("=== Siderust Targets + Proper Motion Example ===\n");

    // Example target: a star-like equatorial mean-J2000 position.
    // (Numbers here are illustrative; use catalog values for real work.)
    let mean_j2000 = position::EquatorialMeanJ2000::<LightYear>::new(
        88.792939 * DEG, // RA
        7.407064 * DEG,  // Dec
        LightYears::new(548.0),
    );

    // Gaia/Hipparcos convention: µα⋆ = µα cos(δ), µδ.
    type MasPerYear = qtty::Per<qtty::MilliArcsecond, qtty::Year>;
    type MasPerYearQ = qtty::Quantity<MasPerYear>;

    let pm = ProperMotion::from_mu_alpha_star::<MasPerYear>(
        MasPerYearQ::new(27.54), // µα⋆ (mas/year)
        MasPerYearQ::new(10.86), // µδ  (mas/year)
    );

    let t0 = JulianDate::J2000;
    let mut target = CoordinateWithPM::new(mean_j2000, t0, pm.clone());

    println!("At J2000.0:");
    println!("  RA  = {:.6} deg", target.position.ra());
    println!("  Dec = {:.6} deg", target.position.dec());
    println!("  d   = {}", target.position.distance);
    println!();

    // Propagate 25 Julian years into the future.
    let t1 = t0 + 25.0 * JulianDate::JULIAN_YEAR;
    let moved = set_proper_motion_since_j2000(target.position, pm, t1)
        .expect("proper motion propagation failed");

    target.update(moved, t1);

    println!("After 25 Julian years:");
    println!("  JD  = {:.6}", target.time);
    println!("  RA  = {:.6} deg", target.position.ra());
    println!("  Dec = {:.6} deg", target.position.dec());
}
