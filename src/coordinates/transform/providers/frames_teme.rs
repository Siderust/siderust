// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-rotation provider for **TEME** (True Equator Mean Equinox).
//!
//! TEME is the output frame of SGP4/SDP4 propagators used with TLE data.
//! It shares the true equator of date but uses a simplified equinox
//! definition (the "equation of the equinoxes" shortcut).
//!
//! The rotation chain is:
//!
//! ```text
//! TEME → EquatorialTrueOfDate → (precession/nutation) → ICRS
//! ```
//!
//! The first step is Rz(Δψ cos εA), which removes the approximate
//! equinox offset.

use super::*;

impl FrameRotationProvider<TEME, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let nut = nutation::nutation_iau2000b(jd);
        Rotation3::rz(nut.dpsi * nut.mean_obliquity.cos())
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, TEME> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, TEME, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<TEME, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        compose_rotation::<TEME, EquatorialTrueOfDate, ICRS, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRS, TEME> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<ICRS, TEME, Eph, Eop, Nut>(jd, ctx)
    }
}
