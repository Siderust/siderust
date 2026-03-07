// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

impl<F> CenterShiftProvider<Heliocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let sun_bary = Eph::sun_barycentric(jd);
        rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
            [
                sun_bary.x().value(),
                sun_bary.y().value(),
                sun_bary.z().value(),
            ],
            jd,
            ctx,
        )
    }
}

impl<F> CenterShiftProvider<Barycentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Barycentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Geocentric, Barycentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        let earth_bary = Eph::earth_barycentric(jd);
        rotate_shift_from_ecliptic::<F, Eph, Eop, Nut>(
            [
                earth_bary.x().value(),
                earth_bary.y().value(),
                earth_bary.z().value(),
            ],
            jd,
            ctx,
        )
    }
}

impl<F> CenterShiftProvider<Barycentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Barycentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Heliocentric, Geocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        compose_shift::<Heliocentric, Barycentric, Geocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}

impl<F> CenterShiftProvider<Geocentric, Heliocentric, F> for ()
where
    F: affn::ReferenceFrame,
    (): FrameRotationProvider<EclipticMeanJ2000, F>,
{
    #[inline]
    fn shift<Eph: Ephemeris, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> [f64; 3] {
        inverse_shift::<Geocentric, Heliocentric, F, Eph, Eop, Nut>(jd, ctx)
    }
}
