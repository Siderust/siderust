// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Mars and Saturn physical orientation (generic body / center / frame API).
//!
//! Run with: `cargo run --example 20_mars_saturn_orientation`
#![allow(clippy::print_stdout)]

use siderust::bodies::solar_system::{Jupiter, Mars};
use siderust::coordinates::centers::{Geocentric, Geodetic, Topocentric};
use siderust::coordinates::frames::{EquatorialTrueOfDate, ECEF, ICRS};
use siderust::ephemeris::Vsop87Ephemeris;
use siderust::event::planetary::{body_orientation, saturn_ring_orientation, saturn_rings};
use siderust::qtty::{Degree, Degrees, Meters};
use siderust::time::JulianDate;

fn main() {
    // Meeus Example 42.a / 45.a epochs (TT ≈ JDE in the book examples).
    let mars_jd = JulianDate::new(2_448_935.500_683);
    let saturn_jd = JulianDate::new(2_448_972.500_68);

    let mars = Mars::orientation(mars_jd);
    println!("Mars (geocentric, true-of-date) @ JD {}", mars_jd);
    if let Some(p) = mars.north_pole_position_angle {
        println!("  north pole PA:  {:.2}", p.to::<Degree>());
    }
    if let Some(b) = mars.bright_limb_position_angle {
        println!("  bright limb PA: {:.2}", b.to::<Degree>());
    }

    let site = Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.48), Meters::new(0.0));
    let mars_topo = Mars::orientation_at::<Topocentric, EquatorialTrueOfDate>((mars_jd, site));
    if let (Some(g), Some(t)) = (
        mars.north_pole_position_angle,
        mars_topo.north_pole_position_angle,
    ) {
        println!(
            "  topo − geo pole PA: {:.4}",
            (t - g).to::<Degree>().wrap_signed()
        );
    }

    let saturn = saturn_ring_orientation::<Vsop87Ephemeris>(saturn_jd);
    println!(
        "\nSaturn rings (geocentric, true-of-date) @ JD {}",
        saturn_jd
    );
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

    let saturn_icrs = saturn_rings::<Vsop87Ephemeris, Geocentric, ICRS>(saturn_jd);
    println!(
        "  opening Earth (ICRS path): {:+.3}  [frame-invariant]",
        saturn_icrs.ring_opening_earth.to::<Degree>()
    );

    let jupiter =
        body_orientation::<Jupiter, Vsop87Ephemeris, Geocentric, EquatorialTrueOfDate>(mars_jd);
    println!(
        "\nJupiter pole PA (generic body API): {:.2}",
        jupiter.north_pole_position_angle.unwrap().to::<Degree>()
    );
}
