// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Center-shift providers for the three **standard** reference centers:
//! [`Barycentric`], [`Heliocentric`] and [`Geocentric`].
//!
//! These form the hub of the center-shift graph. Every other center
//! (planetary, selenocentric, …) composes through [`Barycentric`], so
//! only the three pairwise directions need explicit implementations here.
//!
//! # Data source
//!
//! The Sun and Earth barycentric positions come from the
//! [`Ephemeris`] trait (VSOP87 by default, DE441 optionally).
//! Both return `Position<Barycentric, EclipticMeanJ2000, AU>`.

use super::*;

// ---------------------------------------------------------------------------
// Heliocentric ↔ Barycentric
// ---------------------------------------------------------------------------

/// Heliocentric → Barycentric shift.
///
/// Returns the Sun's barycentric position converted from the ecliptic to
/// the target frame `F`.
impl<F> CenterShiftProvider<Heliocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(Eph::sun_barycentric(jd), jd, ctx)
    }
}

impl<F> CenterShiftProvider<Barycentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        inverse_shift::<Barycentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

// ---------------------------------------------------------------------------
// Geocentric ↔ Barycentric
// ---------------------------------------------------------------------------

/// Geocentric → Barycentric shift.
///
/// Returns the Earth's barycentric position converted from the ecliptic to
/// the target frame `F`.
impl<F> CenterShiftProvider<Geocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(Eph::earth_barycentric(jd), jd, ctx)
    }
}

impl<F> CenterShiftProvider<Barycentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        inverse_shift::<Barycentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

// ---------------------------------------------------------------------------
// Heliocentric ↔ Geocentric (composed through Barycentric)
// ---------------------------------------------------------------------------

impl<F> CenterShiftProvider<Heliocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        compose_shift::<Heliocentric, Barycentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Geocentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        inverse_shift::<Geocentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}
