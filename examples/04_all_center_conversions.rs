// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example: all currently supported center conversions.
//!
//! This demonstrates all center-shift pairs implemented in `providers.rs`:
//! - Barycentric <-> Heliocentric
//! - Barycentric <-> Geocentric
//! - Heliocentric <-> Geocentric
//! - Identity shifts for each center

use qtty::*;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::{Barycentric, Geocentric, Heliocentric, ReferenceCenter};
use siderust::coordinates::frames::EclipticMeanJ2000;
use siderust::coordinates::transform::PositionAstroExt;
use siderust::time::JulianDate;

type F = EclipticMeanJ2000;
type U = AstronomicalUnit;

fn position_error<C: ReferenceCenter<Params = ()>>(
    a: &Position<C, F, U>,
    b: &Position<C, F, U>,
) -> Quantity<U> {
    (a - b).magnitude()
}

fn main() {
    let jd = JulianDate::new(2_460_000.5);

    // A sample physical point in barycentric ecliptic coordinates.
    // We derive equivalent representations in other centers to compare
    // true center-shift roundtrips for the same physical location.

    println!("Center conversion demo at JD(TT) = {:.1}", jd);

    let p_bary = Position::<Barycentric, F, U>::new(0.40, -0.10, 1.20);
    let identity = p_bary.to_center::<Barycentric>(&jd);
    let p_helio = p_bary.to_center::<Heliocentric>(&jd);
    let p_geo = p_bary.to_center::<Geocentric>(&jd);

    println!("{}", p_bary);
    println!(
        "{:.3e}.\tReverse error: {:.3e}",
        identity,
        position_error(&p_bary, &identity.to_center::<Barycentric>(&jd))
    );
    println!(
        "{:.3e}.\tReverse error: {:.3e}",
        p_helio,
        position_error(&p_bary, &p_helio.to_center::<Barycentric>(&jd))
    );
    println!(
        "{:.3e}.\tReverse error: {:.3e}",
        p_geo,
        position_error(&p_bary, &p_geo.to_center::<Barycentric>(&jd))
    );
}
