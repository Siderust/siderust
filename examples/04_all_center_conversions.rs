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
use siderust::coordinates::transform::{CenterShiftProvider, PositionAstroExt};
use siderust::time::JulianDate;

type F = EclipticMeanJ2000;
type U = AstronomicalUnit;

fn show_center_conversion<C1, C2>(jd: &JulianDate, src: &Position<C1, F, U>)
where
    C1: ReferenceCenter<Params = ()>,
    C2: ReferenceCenter<Params = ()>,
    (): CenterShiftProvider<C1, C2, F>,
    (): CenterShiftProvider<C2, C1, F>,
{
    let out: Position<C2, F, U> = src.to_center(jd);
    let back: Position<C1, F, U> = out.to_center(jd);
    let err = (*src - back).magnitude();

    println!(
        "{:<12} -> {:<12} out=({:+.9})  roundtrip={:.3e}",
        C1::center_name(),
        C2::center_name(),
        out,
        err
    );
}

fn main() {
    let jd = JulianDate::new(2_460_000.5);

    // A sample physical point in barycentric ecliptic coordinates.
    // We derive equivalent representations in other centers to compare
    // true center-shift roundtrips for the same physical location.

    println!("Center conversion demo at JD(TT) = {:.1}", jd.value());

    let p_bary = Position::<Barycentric, F, U>::new(0.40, -0.10, 1.20);
    let p_helio: Position<Heliocentric, F, U> = p_bary.to_center(&jd);
    let p_geo: Position<Geocentric, F, U> = p_bary.to_center(&jd);

    // Barycentric source
    show_center_conversion::<Barycentric, Barycentric>(&jd, &p_bary);
    show_center_conversion::<Barycentric, Heliocentric>(&jd, &p_bary);
    show_center_conversion::<Barycentric, Geocentric>(&jd, &p_bary);

    // Heliocentric source
    show_center_conversion::<Heliocentric, Heliocentric>(&jd, &p_helio);
    show_center_conversion::<Heliocentric, Barycentric>(&jd, &p_helio);
    show_center_conversion::<Heliocentric, Geocentric>(&jd, &p_helio);

    // Geocentric source
    show_center_conversion::<Geocentric, Geocentric>(&jd, &p_geo);
    show_center_conversion::<Geocentric, Barycentric>(&jd, &p_geo);
    show_center_conversion::<Geocentric, Heliocentric>(&jd, &p_geo);
}
