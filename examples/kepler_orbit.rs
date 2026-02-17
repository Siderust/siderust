// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Kepler Orbit Example
//!
//! Demonstrates:
//! - Solving Kepler's equation (`solve_keplers_equation`)
//! - Propagating a simple Keplerian orbit with `Orbit::kepler_position`

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::calculus::kepler_equations::solve_keplers_equation;
use siderust::time::JulianDate;

fn main() {
    println!("=== Siderust Kepler Orbit Example ===\n");

    // ---------------------------------------------------------------------
    // 1) Kepler solver (E - e sin E = M)
    // ---------------------------------------------------------------------
    let m = 1.047 * RAD; // mean anomaly
    let e = 0.0167; // eccentricity (Earth-like)
    let e_anom = solve_keplers_equation(m, e);

    println!("Kepler's equation:");
    println!("  M  = {}", m);
    println!("  e  = {:.6}", e);
    println!("  E  = {}", e_anom);
    println!("  residual = {:.3e}\n", (e_anom - Radians::new(e * e_anom.sin()) - m).value());

    // ---------------------------------------------------------------------
    // 2) Simple orbit propagation in heliocentric ecliptic coordinates
    // ---------------------------------------------------------------------
    let earth_like = Orbit::new(
        1.0 * AU,                // a
        0.0167,                  // e
        Degrees::new(0.0),       // i
        Degrees::new(0.0),       // Ω
        Degrees::new(102.94719), // ω (illustrative)
        Degrees::new(100.46435), // M0 at epoch (illustrative)
        JulianDate::J2000,
    );

    let jd0 = JulianDate::J2000;
    let jd1 = jd0 + Days::new(100.0);

    let p0 = earth_like.kepler_position(jd0);
    let p1 = earth_like.kepler_position(jd1);

    println!("Orbit propagation (heliocentric ecliptic Cartesian):");
    println!("  at JD {:.1}: (x,y,z) = ({:+.6}, {:+.6}, {:+.6}) AU", jd0, p0.x(), p0.y(), p0.z());
    println!("  at JD {:.1}: (x,y,z) = ({:+.6}, {:+.6}, {:+.6}) AU", jd1, p1.x(), p1.y(), p1.z());
    println!("  |Δr| = {}", (p1 - p0).magnitude().to::<AstronomicalUnit>());
}
