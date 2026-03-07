// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

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
            ) -> [f64; 3] {
                let bary_pos = <$center>::vsop87e(jd);
                rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
                    [
                        bary_pos.x().value(),
                        bary_pos.y().value(),
                        bary_pos.z().value(),
                    ],
                    jd,
                    ctx,
                )
            }
        }

        impl<F> CenterShiftProvider<Barycentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                inverse_shift::<Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<Heliocentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                compose_shift::<Heliocentric, Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<$center, Heliocentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                inverse_shift::<$center, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<Geocentric, $center, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                compose_shift::<Geocentric, Barycentric, $center, F, Eph, Eop, Nut>(jd, ctx)
            }
        }

        impl<F> CenterShiftProvider<$center, Geocentric, F> for ()
        where
            F: affn::ReferenceFrame,
            (): FrameRotationProvider<EclipticMeanJ2000, F>,
        {
            #[inline]
            fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> [f64; 3] {
                inverse_shift::<$center, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
            }
        }
    };
}

impl_planet_center_shift_vsop!(Mercurycentric);
impl_planet_center_shift_vsop!(Venuscentric);
impl_planet_center_shift_vsop!(Marscentric);
impl_planet_center_shift_vsop!(Jovicentric);
impl_planet_center_shift_vsop!(Saturnocentric);
impl_planet_center_shift_vsop!(Uranocentric);
impl_planet_center_shift_vsop!(Neptunocentric);

impl<F> CenterShiftProvider<Plutocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        use crate::bodies::solar_system;

        let helio_pos = solar_system::PLUTO.orbit.kepler_position(jd);
        let sun_bary = Eph::sun_barycentric(jd);

        rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
            [
                helio_pos.x().value() + sun_bary.x().value(),
                helio_pos.y().value() + sun_bary.y().value(),
                helio_pos.z().value() + sun_bary.z().value(),
            ],
            jd,
            ctx,
        )
    }
}

impl<F> CenterShiftProvider<Barycentric, Plutocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Barycentric, Plutocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Heliocentric, Plutocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        compose_shift::<Heliocentric, Barycentric, Plutocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Plutocentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Plutocentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Geocentric, Plutocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        compose_shift::<Geocentric, Barycentric, Plutocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Plutocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Plutocentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Selenocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let moon_geo = Eph::moon_geocentric(jd);
        let earth_bary = Eph::earth_barycentric(jd);
        let km_per_au = 149_597_870.7;

        rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
            [
                moon_geo.x().value() / km_per_au + earth_bary.x().value(),
                moon_geo.y().value() / km_per_au + earth_bary.y().value(),
                moon_geo.z().value() / km_per_au + earth_bary.z().value(),
            ],
            jd,
            ctx,
        )
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
    ) -> [f64; 3] {
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
    ) -> [f64; 3] {
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
    ) -> [f64; 3] {
        inverse_shift::<Selenocentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Geocentric, Selenocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let moon_geo = Eph::moon_geocentric(jd);
        let km_per_au = 149_597_870.7;

        rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
            [
                -moon_geo.x().value() / km_per_au,
                -moon_geo.y().value() / km_per_au,
                -moon_geo.z().value() / km_per_au,
            ],
            jd,
            ctx,
        )
    }
}

impl<F> CenterShiftProvider<Selenocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Selenocentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}
