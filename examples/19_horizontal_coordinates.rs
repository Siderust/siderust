// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Horizontal Coordinates Example
//!
//! Demonstrates equatorial ↔ horizontal transforms for spherical and cartesian
//! directions and positions, with round-trip checks.
//!
//! Positions use **mean-of-date + GMST** via [`.transform(jd)`]. Directions use
//! **true-of-date + GAST** via [`.to_horizontal`] / [`.to_equatorial`] (observer
//! site passed explicitly). See `doc/frames_and_centers.md` for frame details.
//!
//! Run with: `cargo run --example 19_horizontal_coordinates`
#![allow(clippy::print_stdout)]

use siderust::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Topocentric;
use siderust::coordinates::frames::{EquatorialMeanOfDate, Horizontal};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::{DirectionAstroExt, SphericalDirectionAstroExt, Transform};
use siderust::qtty::{Degrees, Kilometer, KM};
use siderust::time::JulianDate;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
    let jd = JulianDate::new(2_459_015.5);
    let ra = Degrees::new(101.287);
    let dec = Degrees::new(-16.716);

    println!("Observer: {}", ROQUE_DE_LOS_MUCHACHOS.name);
    println!("Epoch (TT): {jd:.6}");

    // Sample topocentric position: mean-of-date ↔ horizontal, using GMST.
    let equatorial = spherical::direction::EquatorialMeanOfDate::new(ra, dec)
        .position_with_params::<Topocentric, Kilometer>(site, 384_400.0 * KM);
    let horizontal: spherical::Position<Topocentric, Horizontal, Kilometer> =
        equatorial.transform(jd);
    let back: spherical::Position<Topocentric, EquatorialMeanOfDate, Kilometer> =
        horizontal.transform(jd);
    println!("\nSpherical position (GMST): {horizontal:.4}");
    println!(
        "Round-trip position error: {}",
        equatorial
            .to_cartesian()
            .try_distance_to(&back.to_cartesian())?
    );

    let equatorial_cart = equatorial.to_cartesian();
    let horizontal_cart: Position<Topocentric, Horizontal, Kilometer> =
        equatorial_cart.transform(jd);
    let back_cart: Position<Topocentric, EquatorialMeanOfDate, Kilometer> =
        horizontal_cart.transform(jd);
    println!("\nCartesian position (GMST): {horizontal_cart:.4}");
    println!(
        "Round-trip position error: {}",
        equatorial_cart.try_distance_to(&back_cart)?
    );

    // Direction transforms use true-of-date coordinates and GAST.
    let direction = spherical::direction::EquatorialTrueOfDate::new(ra, dec);
    let horizontal_direction = direction.to_horizontal(&jd, &site);
    let back_direction = horizontal_direction.to_equatorial(&jd, &site);
    println!("\nSpherical direction (GAST): {horizontal_direction:.4}");
    println!(
        "Round-trip angular error: {}",
        direction.angular_separation(&back_direction)
    );

    let cart_direction = direction.to_cartesian();
    let horizontal_cart_direction = cart_direction.to_horizontal(&jd, &site);
    let back_cart_direction = horizontal_cart_direction.to_equatorial(&jd, &site);
    println!("\nCartesian direction (GAST): {horizontal_cart_direction:.4}");
    println!(
        "Round-trip angular error: {}",
        cart_direction.angle_to(&back_cart_direction)
    );
    Ok(())
}
