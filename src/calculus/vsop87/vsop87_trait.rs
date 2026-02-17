// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::bodies::solar_system::*;
use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::targets::Target;
use crate::time::JulianDate;
use qtty::AstronomicalUnit;

pub trait VSOP87 {
    fn vsop87a(
        &self,
        jd: JulianDate,
    ) -> Target<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>>;
    fn vsop87e(
        &self,
        jd: JulianDate,
    ) -> Target<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>>;
}

macro_rules! impl_vsop87_for_planet {
    ($planet:ident) => {
        impl VSOP87 for $planet {
            fn vsop87a(
                &self,
                jd: JulianDate,
            ) -> Target<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>> {
                $planet::vsop87a(jd)
            }

            fn vsop87e(
                &self,
                jd: JulianDate,
            ) -> Target<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>> {
                $planet::vsop87e(jd)
            }
        }
    };
}

impl_vsop87_for_planet!(Mercury);
impl_vsop87_for_planet!(Venus);
impl_vsop87_for_planet!(Earth);
impl_vsop87_for_planet!(Mars);
impl_vsop87_for_planet!(Jupiter);
impl_vsop87_for_planet!(Saturn);
impl_vsop87_for_planet!(Uranus);
impl_vsop87_for_planet!(Neptune);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::macros::assert_cartesian_eq;
    use crate::time::JulianDate;

    const PRECISION: f64 = 1.0e-12;

    macro_rules! test_dispatch {
        ($planet:ident) => {{
            let body = $planet;
            let dyn_ref: &dyn VSOP87 = &body;
            let jd = JulianDate::J2000;

            let via_trait_a = dyn_ref.vsop87a(jd).get_position().clone();
            let via_trait_e = dyn_ref.vsop87e(jd).get_position().clone();
            let inherent_a = $planet::vsop87a(jd).get_position().clone();
            let inherent_e = $planet::vsop87e(jd).get_position().clone();

            assert_cartesian_eq!(via_trait_a, inherent_a, PRECISION);
            assert_cartesian_eq!(via_trait_e, inherent_e, PRECISION);
        }};
    }

    #[test]
    fn trait_dispatch_for_all_planets() {
        test_dispatch!(Mercury);
        test_dispatch!(Venus);
        test_dispatch!(Earth);
        test_dispatch!(Mars);
        test_dispatch!(Jupiter);
        test_dispatch!(Saturn);
        test_dispatch!(Uranus);
        test_dispatch!(Neptune);
    }
}
