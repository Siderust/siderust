// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Mars and Saturn geocentric physical orientation.
//!
//! Run with: `cargo run --example 20_mars_saturn_orientation`
#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::{Mars, Saturn};
use siderust::ephemeris::Vsop87Ephemeris;
use siderust::event::planetary::{mars_orientation, saturn_ring_orientation};
use siderust::qtty::{Degree, Degrees};
use siderust::time::JulianDate;

fn main() {
    // Meeus Example 42.a / 45.a epochs (TT ≈ JDE in the book examples).
    let mars_jd = JulianDate::new(2_448_935.500_683);
    let saturn_jd = JulianDate::new(2_448_972.500_68);

    let mars = Mars::orientation(mars_jd);
    println!("Mars orientation @ JD {}", mars_jd);
    if let Some(p) = mars.north_pole_position_angle {
        println!("  north pole PA:  {:.2}", p.to::<Degree>());
    }
    if let Some(b) = mars.bright_limb_position_angle {
        println!("  bright limb PA: {:.2}", b.to::<Degree>());
    }

    let saturn = saturn_ring_orientation::<Vsop87Ephemeris>(saturn_jd);
    println!("\nSaturn rings @ JD {}", saturn_jd);
    if let Some(p) = saturn.north_pole_position_angle {
        println!("  ring-axis PA:   {:.3}", p.to::<Degree>());
    }
    println!(
        "  opening Earth:  {:+.3}",
        saturn.ring_opening_earth.to::<Degree>()
    );
    println!(
        "  opening Sun:    {:+.3}",
        saturn.ring_opening_sun.to::<Degree>()
    );

    // Expert backend path matches the default wrapper for VSOP87.
    let mars_e = mars_orientation::<Vsop87Ephemeris>(mars_jd);
    assert_eq!(mars, mars_e);
    let _ = Saturn::ring_orientation(siderust::J2000);
    let _ = Degrees::new(0.0);
}
