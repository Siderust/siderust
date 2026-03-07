// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Center-shift providers for **planetary**, **plutocentric** and
//! **selenocentric** reference centers.
//!
//! # VSOP87 planets (Mercury – Neptune)
//!
//! Each planet center type (e.g. [`Mercurycentric`]) exposes `vsop87e(jd)`
//! which returns a typed `Position<Barycentric, EclipticMeanJ2000, AU>`.
//! The [`impl_planet_center_shift_vsop`] macro generates the direct
//! `Center → Barycentric` shift plus all reverse and composed impls via
//! [`impl_reverse_center_shifts`].
//!
//! # Pluto
//!
//! Pluto has no VSOP87 series. Its heliocentric position is obtained from
//! a Keplerian orbit, then combined with the Sun's barycentric position.
//!
//! # Moon (Selenocentric)
//!
//! The Moon's geocentric position is provided by the ephemeris in
//! kilometres and must be converted to AU before combining with the
//! Earth's barycentric position. The conversion uses
//! [`qtty`]'s typed unit conversion (`to_unit::<AstronomicalUnit>()`),
//! avoiding any hardcoded magic constants.

use super::*;

// ---------------------------------------------------------------------------
// VSOP87 planets (Mercury – Neptune)
// ---------------------------------------------------------------------------

/// Generates all six `CenterShiftProvider` impls for a VSOP87 planet center.
///
/// The direct `Center → Barycentric` impl calls `vsop87e(jd)` which returns
/// a typed `Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>` and
/// passes it directly to [`rotate_shift_from_ecliptic`]. The five
/// reverse/composed impls are delegated to [`impl_reverse_center_shifts`].
macro_rules! impl_planet_center_shift_vsop {
    ($center:ty) => {
        impl<F> CenterShiftProvider<$center, Barycentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> AuShift {
                rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(<$center>::vsop87e(jd), jd, ctx)
            }
        }

        impl_reverse_center_shifts!($center);
    };
}

impl_planet_center_shift_vsop!(Mercurycentric);
impl_planet_center_shift_vsop!(Venuscentric);
impl_planet_center_shift_vsop!(Marscentric);
impl_planet_center_shift_vsop!(Jovicentric);
impl_planet_center_shift_vsop!(Saturnocentric);
impl_planet_center_shift_vsop!(Uranocentric);
impl_planet_center_shift_vsop!(Neptunocentric);

// ---------------------------------------------------------------------------
// Pluto (Keplerian orbit + Sun barycentric)
// ---------------------------------------------------------------------------

/// Pluto → Barycentric shift.
///
/// Pluto lacks a VSOP87 series. Instead its heliocentric ecliptic position
/// is computed from a Keplerian orbit and offset by the Sun's barycentric
/// position to yield the full barycentric vector.
impl<F> CenterShiftProvider<Plutocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        use crate::bodies::solar_system;

        let helio_pos = solar_system::PLUTO.orbit.kepler_position(jd);
        let sun_bary = Eph::sun_barycentric(jd);
        // Combine the heliocentric Keplerian position with the Sun's
        // barycentric offset using AU quantity arithmetic — no raw f64 needed.
        let bary_pos = Position::<Barycentric, EclipticMeanJ2000, qtty::AstronomicalUnit>::new(
            helio_pos.x() + sun_bary.x(),
            helio_pos.y() + sun_bary.y(),
            helio_pos.z() + sun_bary.z(),
        );

        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(bary_pos, jd, ctx)
    }
}

impl_reverse_center_shifts!(Plutocentric);

// ---------------------------------------------------------------------------
// Moon (Selenocentric)
// ---------------------------------------------------------------------------

/// Selenocentric → Barycentric shift.
///
/// The Moon's geocentric position (in km) is converted to AU using
/// `qtty`'s typed unit conversion, then added to the Earth's barycentric
/// position (already in AU).
impl<F> CenterShiftProvider<Selenocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        let moon_geo_au = Eph::moon_geocentric(jd).to_unit::<qtty::AstronomicalUnit>();
        let earth_bary = Eph::earth_barycentric(jd);
        // Combine geocentric Moon (now in AU) with Earth's barycentric offset
        // using AU quantity arithmetic — no raw f64 or magic constants needed.
        let seleno_bary = Position::<Barycentric, EclipticMeanJ2000, qtty::AstronomicalUnit>::new(
            moon_geo_au.x() + earth_bary.x(),
            moon_geo_au.y() + earth_bary.y(),
            moon_geo_au.z() + earth_bary.z(),
        );

        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(seleno_bary, jd, ctx)
    }
}

/// Geocentric → Selenocentric override (direct, avoids double composition).
///
/// Negates the Moon's geocentric position in AU directly, which is more
/// efficient than routing through Barycentric.
impl<F> CenterShiftProvider<Geocentric, Selenocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        let moon_geo_au = Eph::moon_geocentric(jd).to_unit::<qtty::AstronomicalUnit>();
        // The Geo→Seleno shift is the negated Moon geocentric position:
        // to move from Earth-centred to Moon-centred, subtract moon's location.
        let geo_to_seleno = Position::<Geocentric, EclipticMeanJ2000, qtty::AstronomicalUnit>::new(
            -moon_geo_au.x(),
            -moon_geo_au.y(),
            -moon_geo_au.z(),
        );

        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(geo_to_seleno, jd, ctx)
    }
}

impl<F> CenterShiftProvider<Barycentric, Selenocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        inverse_shift::<Barycentric, Selenocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Heliocentric, Selenocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        compose_shift::<Heliocentric, Barycentric, Selenocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Selenocentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        inverse_shift::<Selenocentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

/// Selenocentric → Geocentric override (inverse of the direct Geo→Seleno).
impl<F> CenterShiftProvider<Selenocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> AuShift {
        inverse_shift::<Selenocentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}
