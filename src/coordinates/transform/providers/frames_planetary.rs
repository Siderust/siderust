// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

macro_rules! impl_body_fixed_rotation {
    ($frame:ty, $body_ty:ty) => {
        impl FrameRotationProvider<$frame, ICRS> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut>(
                jd: JulianDate,
                _ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                iau_body_fixed_to_icrs(&<$body_ty as HasIauRotation>::ROTATION, jd)
            }
        }

        impl FrameRotationProvider<ICRS, $frame> for () {
            #[inline]
            fn rotation<Eph, Eop: EopProvider, Nut>(
                jd: JulianDate,
                ctx: &AstroContext<Eph, Eop, Nut>,
            ) -> Rotation3 {
                inverse_rotation::<ICRS, $frame, Eph, Eop, Nut>(jd, ctx)
            }
        }
    };
}

impl_body_fixed_rotation!(MercuryFixed, crate::bodies::solar_system::Mercury);
impl_body_fixed_rotation!(VenusFixed, crate::bodies::solar_system::Venus);
impl_body_fixed_rotation!(MarsFixed, crate::bodies::solar_system::Mars);
impl_body_fixed_rotation!(MoonPrincipalAxes, crate::bodies::solar_system::Moon);
impl_body_fixed_rotation!(JupiterSystemIII, crate::bodies::solar_system::Jupiter);
impl_body_fixed_rotation!(SaturnFixed, crate::bodies::solar_system::Saturn);
impl_body_fixed_rotation!(UranusFixed, crate::bodies::solar_system::Uranus);
impl_body_fixed_rotation!(NeptuneFixed, crate::bodies::solar_system::Neptune);
impl_body_fixed_rotation!(PlutoFixed, crate::bodies::solar_system::Pluto);

impl_via_icrs_bidirectional!(EclipticMeanJ2000, MercuryFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, VenusFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, MarsFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, MoonPrincipalAxes);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, JupiterSystemIII);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, SaturnFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, UranusFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, NeptuneFixed);
impl_via_icrs_bidirectional!(EclipticMeanJ2000, PlutoFixed);

impl_via_icrs_bidirectional!(EME2000, MercuryFixed);
impl_via_icrs_bidirectional!(EME2000, VenusFixed);
impl_via_icrs_bidirectional!(EME2000, MarsFixed);
impl_via_icrs_bidirectional!(EME2000, MoonPrincipalAxes);
impl_via_icrs_bidirectional!(EME2000, JupiterSystemIII);
impl_via_icrs_bidirectional!(EME2000, SaturnFixed);
impl_via_icrs_bidirectional!(EME2000, UranusFixed);
impl_via_icrs_bidirectional!(EME2000, NeptuneFixed);
impl_via_icrs_bidirectional!(EME2000, PlutoFixed);
