// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example: all currently supported frame conversions.
//!
//! This demonstrates:
//! - Every direct frame-rotation pair implemented in `providers.rs`.
//! - Specialized time-dependent frame conversions:
//!   - EquatorialMeanOfDate <-> EclipticTrueOfDate
//!   - EquatorialTrueOfDate <-> Horizontal

use qtty::*;
use siderust::coordinates::cartesian::Direction;
use siderust::coordinates::centers::ObserverSite;
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EclipticTrueOfDate, EquatorialMeanJ2000, EquatorialMeanOfDate,
    EquatorialTrueOfDate, Horizontal, ReferenceFrame, ICRF, ICRS,
};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::ecliptic_of_date::FromEclipticTrueOfDate;
use siderust::coordinates::transform::FromHorizontal;
use siderust::coordinates::transform::{DirectionAstroExt, FrameRotationProvider};
use siderust::time::JulianDate;

fn direction_error<F: ReferenceFrame>(a: &Direction<F>, b: &Direction<F>) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    let dz = a.z() - b.z();
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn show_provider_conversion<F1, F2>(jd: &JulianDate, src: &Direction<F1>)
where
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    (): FrameRotationProvider<F1, F2>,
    (): FrameRotationProvider<F2, F1>,
{
    let out: Direction<F2> = src.to_frame(jd);
    let back: Direction<F1> = out.to_frame(jd);
    let err = direction_error(src, &back);

    println!(
        "{:<24} -> {:<24} out=({:+.9}, {:+.9}, {:+.9})  roundtrip={:.3e}",
        F1::frame_name(),
        F2::frame_name(),
        out.x(),
        out.y(),
        out.z(),
        err
    );
}

fn show_provider_conversions(jd: &JulianDate) {
    println!("=== Provider-Based Frame Conversions ===");

    let d_icrs = Direction::<ICRS>::normalize(0.3, -0.7, 0.64);
    let d_icrf: Direction<ICRF> = d_icrs.to_frame(jd);
    let d_ecl: Direction<EclipticMeanJ2000> = d_icrs.to_frame(jd);
    let d_eq_j2000: Direction<EquatorialMeanJ2000> = d_icrs.to_frame(jd);
    let d_eq_mod: Direction<EquatorialMeanOfDate> = d_icrs.to_frame(jd);
    let d_eq_tod: Direction<EquatorialTrueOfDate> = d_icrs.to_frame(jd);

    macro_rules! c {
        ($src:ident, $from:ty => $to:ty) => {
            show_provider_conversion::<$from, $to>(jd, &$src);
        };
    }

    // Identity
    c!(d_icrs, ICRS => ICRS);
    c!(d_icrf, ICRF => ICRF);
    c!(d_ecl, EclipticMeanJ2000 => EclipticMeanJ2000);
    c!(d_eq_j2000, EquatorialMeanJ2000 => EquatorialMeanJ2000);
    c!(d_eq_mod, EquatorialMeanOfDate => EquatorialMeanOfDate);
    c!(d_eq_tod, EquatorialTrueOfDate => EquatorialTrueOfDate);

    // All direct non-identity provider pairs
    c!(d_icrs, ICRS => EclipticMeanJ2000);
    c!(d_ecl, EclipticMeanJ2000 => ICRS);
    c!(d_icrs, ICRS => EquatorialMeanJ2000);
    c!(d_eq_j2000, EquatorialMeanJ2000 => ICRS);
    c!(d_eq_j2000, EquatorialMeanJ2000 => EclipticMeanJ2000);
    c!(d_ecl, EclipticMeanJ2000 => EquatorialMeanJ2000);
    c!(d_eq_j2000, EquatorialMeanJ2000 => EquatorialMeanOfDate);
    c!(d_eq_mod, EquatorialMeanOfDate => EquatorialMeanJ2000);
    c!(d_eq_mod, EquatorialMeanOfDate => EquatorialTrueOfDate);
    c!(d_eq_tod, EquatorialTrueOfDate => EquatorialMeanOfDate);
    c!(d_eq_j2000, EquatorialMeanJ2000 => EquatorialTrueOfDate);
    c!(d_eq_tod, EquatorialTrueOfDate => EquatorialMeanJ2000);
    c!(d_icrs, ICRS => EquatorialMeanOfDate);
    c!(d_eq_mod, EquatorialMeanOfDate => ICRS);
    c!(d_icrs, ICRS => EquatorialTrueOfDate);
    c!(d_eq_tod, EquatorialTrueOfDate => ICRS);
    c!(d_icrf, ICRF => ICRS);
    c!(d_icrs, ICRS => ICRF);
    c!(d_icrf, ICRF => EquatorialMeanJ2000);
    c!(d_eq_j2000, EquatorialMeanJ2000 => ICRF);
    c!(d_icrf, ICRF => EclipticMeanJ2000);
    c!(d_ecl, EclipticMeanJ2000 => ICRF);
    c!(d_icrf, ICRF => EquatorialMeanOfDate);
    c!(d_eq_mod, EquatorialMeanOfDate => ICRF);
    c!(d_icrf, ICRF => EquatorialTrueOfDate);
    c!(d_eq_tod, EquatorialTrueOfDate => ICRF);
}

fn show_specialized_conversions(jd_tt: &JulianDate) {
    println!("\n=== Specialized Frame Conversions ===");

    // EquatorialMeanOfDate <-> EclipticTrueOfDate
    let eq_mod =
        spherical::Direction::<EquatorialMeanOfDate>::new(120.0 * DEG, 23.0 * DEG).to_cartesian();
    let ecl_of_date: Direction<EclipticTrueOfDate> = eq_mod.to_ecliptic_of_date(jd_tt);
    let eq_mod_back: Direction<EquatorialMeanOfDate> =
        ecl_of_date.to_equatorial_mean_of_date(jd_tt);
    println!(
        "{:24} -> {:24} roundtrip={:.3e}",
        EquatorialMeanOfDate::frame_name(),
        EclipticTrueOfDate::frame_name(),
        direction_error(&eq_mod, &eq_mod_back)
    );

    // EquatorialTrueOfDate <-> Horizontal (observer-dependent)
    let site = ObserverSite::new(2.3522 * DEG, 48.8566 * DEG, 35.0 * M);
    let jd_ut1 = siderust::astro::earth_rotation::jd_ut1_from_tt(*jd_tt);
    let eq_tod =
        spherical::Direction::<EquatorialTrueOfDate>::new(210.0 * DEG, -10.0 * DEG).to_cartesian();
    let horiz: Direction<Horizontal> = eq_tod.to_horizontal(jd_tt, &site);
    let eq_tod_back: Direction<EquatorialTrueOfDate> = horiz.to_equatorial(&jd_ut1, jd_tt, &site);
    println!(
        "{:24} -> {:24} roundtrip={:.3e}",
        EquatorialTrueOfDate::frame_name(),
        Horizontal::frame_name(),
        direction_error(&eq_tod, &eq_tod_back)
    );
}

fn main() {
    let jd_tt = JulianDate::new(2_460_000.5);

    println!("Frame conversion demo at JD(TT) = {:.1}", jd_tt.value());
    show_provider_conversions(&jd_tt);
    show_specialized_conversions(&jd_tt);
}
