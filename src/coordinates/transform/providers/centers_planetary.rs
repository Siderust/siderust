// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Center-shift providers for **planetary**, **plutocentric** and
//! **selenocentric** reference centers.
//!
//! # Planet centers (Mercury – Neptune)
//!
//! Translations request the selected [`Ephemeris::try_planet_barycentric`]
//! capability. Missing capabilities fail explicitly; there is no fallback to
//! another backend. Reverse and composed shifts use the same position.
//!
//! # Pluto
//!
//! Pluto has no VSOP87 series. Its heliocentric position is obtained from
//! the selected backend. The default backend retains its coarse Keplerian
//! orbit plus Sun model; other backends must opt in explicitly.
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

/// Generate direct planetary-center shifts and their reverse/composed paths.
macro_rules! impl_planet_center_shift {
    ($center:ty, $planet:ident) => {
        impl<F> CenterShiftProvider<$center, Barycentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop>,
            ) -> AuShift {
                let position =
                    Eph::try_planet_barycentric(crate::ephemeris::MajorPlanet::$planet, jd)
                        .expect("selected ephemeris cannot supply the requested planet center");
                rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(position, jd, ctx)
            }
        }

        impl_reverse_center_shifts!($center);
    };
}

impl_planet_center_shift!(Mercurycentric, Mercury);
impl_planet_center_shift!(Venuscentric, Venus);
impl_planet_center_shift!(Marscentric, Mars);
impl_planet_center_shift!(Jovicentric, Jupiter);
impl_planet_center_shift!(Saturnocentric, Saturn);
impl_planet_center_shift!(Uranocentric, Uranus);
impl_planet_center_shift!(Neptunocentric, Neptune);

// ---------------------------------------------------------------------------
// Pluto (Keplerian orbit + Sun barycentric)
// ---------------------------------------------------------------------------

/// Pluto → Barycentric shift.
///
/// Uses the selected backend's explicit Pluto capability. The default backend
/// retains its catalogue Keplerian approximation, not a high precision DE state.
impl<F> CenterShiftProvider<Plutocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        let bary_pos = Eph::pluto_barycentric(jd);
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
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        let moon_geo_au = Eph::moon_geocentric(jd).to_unit::<crate::qtty::AstronomicalUnit>();
        let earth_bary = Eph::earth_barycentric(jd);
        // Combine geocentric Moon (now in AU) with Earth's barycentric offset
        // using AU quantity arithmetic, no raw f64 or magic constants needed.
        let seleno_bary =
            Position::<Barycentric, EclipticMeanJ2000, crate::qtty::AstronomicalUnit>::new(
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
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        let moon_geo_au = Eph::moon_geocentric(jd).to_unit::<crate::qtty::AstronomicalUnit>();
        // The Geo→Seleno shift is the negated Moon geocentric position:
        // to move from Earth-centred to Moon-centred, subtract moon's location.
        let geo_to_seleno =
            Position::<Geocentric, EclipticMeanJ2000, crate::qtty::AstronomicalUnit>::new(
                -moon_geo_au.x(),
                -moon_geo_au.y(),
                -moon_geo_au.z(),
            );

        rotate_shift_from_ecliptic::<_, F, Eph, Eop, Nut>(geo_to_seleno, jd, ctx)
    }
}

impl_reverse_center_shifts_without_geocentric!(Selenocentric);

/// Selenocentric → Geocentric override (inverse of the direct Geo→Seleno).
impl<F> CenterShiftProvider<Selenocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut: NutationModel>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop>,
    ) -> AuShift {
        inverse_shift::<Selenocentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}
