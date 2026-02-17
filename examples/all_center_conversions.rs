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

fn position_error<C: ReferenceCenter<Params = ()>>(
    a: &Position<C, F, U>,
    b: &Position<C, F, U>,
) -> f64 {
    let dx = (a.x() - b.x()).value();
    let dy = (a.y() - b.y()).value();
    let dz = (a.z() - b.z()).value();
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn show_center_conversion<C1, C2>(jd: &JulianDate, src: &Position<C1, F, U>)
where
    C1: ReferenceCenter<Params = ()>,
    C2: ReferenceCenter<Params = ()>,
    (): CenterShiftProvider<C1, C2, F>,
    (): CenterShiftProvider<C2, C1, F>,
{
    let out: Position<C2, F, U> = src.to_center(jd);
    let back: Position<C1, F, U> = out.to_center(jd);
    let err = position_error(src, &back);

    println!(
        "{:<12} -> {:<12} out=({:+.9}, {:+.9}, {:+.9})  roundtrip={:.3e}",
        C1::center_name(),
        C2::center_name(),
        out.x().value(),
        out.y().value(),
        out.z().value(),
        err
    );
}

fn main() {
    let jd = JulianDate::new(2_460_000.5);

    // A sample physical point in barycentric ecliptic coordinates.
    // We derive equivalent representations in other centers to compare
    // true center-shift roundtrips for the same physical location.
    let p_bary = Position::<Barycentric, F, U>::new(0.40, -0.10, 1.20);
    let p_helio: Position<Heliocentric, F, U> = p_bary.to_center(&jd);
    let p_geo: Position<Geocentric, F, U> = p_bary.to_center(&jd);

    println!("Center conversion demo at JD(TT) = {:.1}", jd.value());

    macro_rules! c {
        ($src:ident, $from:ty => $to:ty) => {
            show_center_conversion::<$from, $to>(&jd, &$src);
        };
    }

    // Barycentric source
    c!(p_bary, Barycentric => Barycentric);
    c!(p_bary, Barycentric => Heliocentric);
    c!(p_bary, Barycentric => Geocentric);

    // Heliocentric source
    c!(p_helio, Heliocentric => Heliocentric);
    c!(p_helio, Heliocentric => Barycentric);
    c!(p_helio, Heliocentric => Geocentric);

    // Geocentric source
    c!(p_geo, Geocentric => Geocentric);
    c!(p_geo, Geocentric => Barycentric);
    c!(p_geo, Geocentric => Heliocentric);
}
