// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame-rotation providers for **inertial and quasi-inertial** frames.
//!
//! This module covers rotations between ICRS (the hub), the J2000 ecliptic
//! and equatorial frames, the mean/true equators of date, and the
//! identity-equivalent frames ICRF and GCRS.
//!
//! Key rotation chains:
//!
//! - **Frame bias**: ICRS ↔ EquatorialMeanJ2000 (constant ~80 µas matrix).
//! - **Obliquity**: EquatorialMeanJ2000 ↔ EclipticMeanJ2000 (Rx by −ε₀).
//! - **Precession**: EquatorialMeanJ2000 → EquatorialMeanOfDate (IAU 2006).
//! - **Nutation**: EquatorialMeanOfDate → EquatorialTrueOfDate (IAU 2000B).
//! - **EME2000 ≡ EquatorialMeanJ2000**: identity alias.
//! - **ICRF ≡ ICRS**: identity alias.
//! - **GCRS ≈ ICRS**: treated as identity (no aberration modelling).

use super::*;

/// ICRS → EclipticMeanJ2000 rotation (J2000 mean ecliptic).
impl FrameRotationProvider<ICRS, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::rx(-j2000_obliquity()) * FRAME_BIAS_ICRS_TO_J2000
    }
}

impl FrameRotationProvider<EclipticMeanJ2000, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EclipticMeanJ2000, ICRS, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRS, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FRAME_BIAS_ICRS_TO_J2000
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FRAME_BIAS_ICRS_TO_J2000.inverse()
    }
}

impl FrameRotationProvider<EME2000, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanJ2000, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, ICRS>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<ICRS, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, EclipticMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EclipticMeanJ2000, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EclipticMeanJ2000, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, EquatorialMeanOfDate>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanOfDate, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanOfDate, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, EquatorialTrueOfDate>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, ICRF>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<ICRF, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EME2000, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, GCRSFrame>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<GCRSFrame, EME2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<GCRSFrame, EME2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::rx(-j2000_obliquity())
    }
}

impl FrameRotationProvider<EclipticMeanJ2000, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::rx(j2000_obliquity())
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        precession::precession_matrix_iau2006(jd)
    }
}

impl FrameRotationProvider<EquatorialMeanOfDate, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        precession::precession_matrix_iau2006(jd).inverse()
    }
}

impl FrameRotationProvider<EquatorialMeanOfDate, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let nut = nutation::nutation_iau2000b(jd);
        Rotation3::fused_rx_rz_rx(nut.mean_obliquity + nut.deps, nut.dpsi, -nut.mean_obliquity)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, EquatorialMeanOfDate, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        compose_rotation::<
            EquatorialMeanJ2000,
            EquatorialMeanOfDate,
            EquatorialTrueOfDate,
            Eph,
            Eop,
            Nut,
        >(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, EquatorialMeanJ2000, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRS, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        compose_rotation::<ICRS, EquatorialMeanJ2000, EquatorialMeanOfDate, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanOfDate, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanOfDate, ICRS, Eph, Eop, Nut>(jd, ctx)
    }
}

/// ICRS → EquatorialTrueOfDate rotation.
///
/// Combines IAU 2000B nutation with IERS celestial-pole corrections
/// (dX, dY) from the EOP provider into a full precession–nutation matrix.
///
/// The first-order correction is:
///   dψ_eop = dX / sin(εA),  dε_eop = dY
impl FrameRotationProvider<ICRS, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        let nut = nutation::nutation_iau2000b(jd);
        let mut dpsi = nut.dpsi;
        let mut deps = nut.deps;
        let eop = ctx.eop_at(jd);
        let dx_rad = qtty::Radians::from(eop.dx);
        let dy_rad = qtty::Radians::from(eop.dy);

        let zero = qtty::Radians::new(0.0);
        if dx_rad != zero || dy_rad != zero {
            let sin_eps = nut.mean_obliquity.sin();
            if sin_eps.abs() > 1e-15 {
                dpsi += dx_rad / sin_eps;
            }
            deps += dy_rad;
        }

        precession::precession_nutation_matrix(jd, dpsi, deps)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, ICRS, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

impl FrameRotationProvider<ICRS, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<ICRS, ICRF, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanJ2000, ICRF, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EclipticMeanJ2000, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EclipticMeanJ2000, ICRF, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, EquatorialMeanOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanOfDate>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanOfDate, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanOfDate, ICRF, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<ICRF, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, ICRF> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialTrueOfDate, ICRF, Eph, Eop, Nut>(jd, ctx)
    }
}

impl FrameRotationProvider<GCRSFrame, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

impl FrameRotationProvider<ICRS, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        Rotation3::IDENTITY
    }
}

impl FrameRotationProvider<GCRSFrame, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialMeanJ2000, ICRS>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<GCRSFrame, EquatorialTrueOfDate> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EquatorialTrueOfDate, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EquatorialTrueOfDate, ICRS>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<GCRSFrame, EclipticMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<ICRS, EclipticMeanJ2000>>::rotation(jd, ctx)
    }
}

impl FrameRotationProvider<EclipticMeanJ2000, GCRSFrame> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        <() as FrameRotationProvider<EclipticMeanJ2000, ICRS>>::rotation(jd, ctx)
    }
}
