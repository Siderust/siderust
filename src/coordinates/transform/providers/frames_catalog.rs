// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

const GALACTIC_TO_ICRS: Rotation3 = Rotation3::from_matrix([
    [
        -0.054_875_560_416_215_4,
        -0.873_437_090_234_885_1,
        -0.483_835_015_548_713_2,
    ],
    [
        0.494_109_427_875_583_7,
        -0.444_829_629_960_011_2,
        0.746_982_244_580_286_6,
    ],
    [
        -0.867_666_149_019_004_7,
        -0.198_076_373_431_201_5,
        0.455_983_776_175_066_9,
    ],
]);

impl FrameRotationProvider<Galactic, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        GALACTIC_TO_ICRS
    }
}

impl FrameRotationProvider<ICRS, Galactic> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        GALACTIC_TO_ICRS.inverse()
    }
}

impl_via_icrs_bidirectional!(Galactic, EquatorialMeanJ2000);
impl_via_icrs_bidirectional!(Galactic, EclipticMeanJ2000);

const FK4_TO_FK5: Rotation3 = Rotation3::from_matrix([
    [
        0.999_925_679_495_687_7,
        -0.011_181_483_220_466_2,
        -0.004_859_003_815_359_2,
    ],
    [
        0.011_181_483_239_171_7,
        0.999_937_484_893_313_5,
        -0.000_027_162_594_714_2,
    ],
    [
        0.004_859_003_772_314_3,
        -0.000_027_170_293_744_0,
        0.999_988_194_602_374_2,
    ],
]);

impl FrameRotationProvider<FK4B1950, ICRS> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FRAME_BIAS_ICRS_TO_J2000.inverse() * FK4_TO_FK5
    }
}

impl FrameRotationProvider<ICRS, FK4B1950> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FK4_TO_FK5.inverse() * FRAME_BIAS_ICRS_TO_J2000
    }
}

impl FrameRotationProvider<FK4B1950, EquatorialMeanJ2000> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        _jd: JulianDate,
        _ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        FK4_TO_FK5
    }
}

impl FrameRotationProvider<EquatorialMeanJ2000, FK4B1950> for () {
    #[inline]
    fn rotation<Eph, Eop: EopProvider, Nut>(
        jd: JulianDate,
        ctx: &AstroContext<Eph, Eop, Nut>,
    ) -> Rotation3 {
        inverse_rotation::<EquatorialMeanJ2000, FK4B1950, Eph, Eop, Nut>(jd, ctx)
    }
}
