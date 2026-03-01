// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example: all currently supported frame conversions.
//!
//! This demonstrates every direct frame-rotation pair implemented in
//! `providers.rs`, plus identity rotations for each frame.

use qtty::*;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::Barycentric;
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate,
    ReferenceFrame, ICRF, ICRS,
};
use siderust::coordinates::transform::{FrameRotationProvider, PositionAstroExt};
use siderust::time::JulianDate;

type C = Barycentric;
type U = AstronomicalUnit;

fn show_frame_conversion<F1, F2>(jd: &JulianDate, src: &Position<C, F1, U>)
where
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    (): FrameRotationProvider<F1, F2>,
    (): FrameRotationProvider<F2, F1>,
{
    let out: Position<C, F2, U> = src.to_frame(jd);
    let back: Position<C, F1, U> = out.to_frame(jd);
    let err = (*src - back).magnitude();

    println!(
        "{:<24} -> {:<24} out=({:+.9})  roundtrip={:.3e}",
        F1::frame_name(),
        F2::frame_name(),
        out,
        err
    );
}

fn main() {
    let jd = JulianDate::new(2_460_000.5);
    println!("Frame conversion demo at JD(TT) = {:.1}", jd.value());

    let p_icrs = Position::<C, ICRS, U>::new(0.30, -0.70, 0.64);
    let p_icrf: Position<C, ICRF, U> = p_icrs.to_frame(&jd);
    let p_ecl: Position<C, EclipticMeanJ2000, U> = p_icrs.to_frame(&jd);
    let p_eq_j2000: Position<C, EquatorialMeanJ2000, U> = p_icrs.to_frame(&jd);
    let p_eq_mod: Position<C, EquatorialMeanOfDate, U> = p_icrs.to_frame(&jd);
    let p_eq_tod: Position<C, EquatorialTrueOfDate, U> = p_icrs.to_frame(&jd);

    // Identity
    show_frame_conversion::<ICRS, ICRS>(&jd, &p_icrs);
    show_frame_conversion::<ICRF, ICRF>(&jd, &p_icrf);
    show_frame_conversion::<EclipticMeanJ2000, EclipticMeanJ2000>(&jd, &p_ecl);
    show_frame_conversion::<EquatorialMeanJ2000, EquatorialMeanJ2000>(&jd, &p_eq_j2000);
    show_frame_conversion::<EquatorialMeanOfDate, EquatorialMeanOfDate>(&jd, &p_eq_mod);
    show_frame_conversion::<EquatorialTrueOfDate, EquatorialTrueOfDate>(&jd, &p_eq_tod);

    // All direct non-identity provider pairs
    show_frame_conversion::<ICRS, EclipticMeanJ2000>(&jd, &p_icrs);
    show_frame_conversion::<EclipticMeanJ2000, ICRS>(&jd, &p_ecl);
    show_frame_conversion::<ICRS, EquatorialMeanJ2000>(&jd, &p_icrs);
    show_frame_conversion::<EquatorialMeanJ2000, ICRS>(&jd, &p_eq_j2000);
    show_frame_conversion::<EquatorialMeanJ2000, EclipticMeanJ2000>(&jd, &p_eq_j2000);
    show_frame_conversion::<EclipticMeanJ2000, EquatorialMeanJ2000>(&jd, &p_ecl);
    show_frame_conversion::<EquatorialMeanJ2000, EquatorialMeanOfDate>(&jd, &p_eq_j2000);
    show_frame_conversion::<EquatorialMeanOfDate, EquatorialMeanJ2000>(&jd, &p_eq_mod);
    show_frame_conversion::<EquatorialMeanOfDate, EquatorialTrueOfDate>(&jd, &p_eq_mod);
    show_frame_conversion::<EquatorialTrueOfDate, EquatorialMeanOfDate>(&jd, &p_eq_tod);
    show_frame_conversion::<EquatorialMeanJ2000, EquatorialTrueOfDate>(&jd, &p_eq_j2000);
    show_frame_conversion::<EquatorialTrueOfDate, EquatorialMeanJ2000>(&jd, &p_eq_tod);
    show_frame_conversion::<ICRS, EquatorialMeanOfDate>(&jd, &p_icrs);
    show_frame_conversion::<EquatorialMeanOfDate, ICRS>(&jd, &p_eq_mod);
    show_frame_conversion::<ICRS, EquatorialTrueOfDate>(&jd, &p_icrs);
    show_frame_conversion::<EquatorialTrueOfDate, ICRS>(&jd, &p_eq_tod);
    show_frame_conversion::<ICRF, ICRS>(&jd, &p_icrf);
    show_frame_conversion::<ICRS, ICRF>(&jd, &p_icrs);
    show_frame_conversion::<ICRF, EquatorialMeanJ2000>(&jd, &p_icrf);
    show_frame_conversion::<EquatorialMeanJ2000, ICRF>(&jd, &p_eq_j2000);
    show_frame_conversion::<ICRF, EclipticMeanJ2000>(&jd, &p_icrf);
    show_frame_conversion::<EclipticMeanJ2000, ICRF>(&jd, &p_ecl);
    show_frame_conversion::<ICRF, EquatorialMeanOfDate>(&jd, &p_icrf);
    show_frame_conversion::<EquatorialMeanOfDate, ICRF>(&jd, &p_eq_mod);
    show_frame_conversion::<ICRF, EquatorialTrueOfDate>(&jd, &p_icrf);
    show_frame_conversion::<EquatorialTrueOfDate, ICRF>(&jd, &p_eq_tod);
}
