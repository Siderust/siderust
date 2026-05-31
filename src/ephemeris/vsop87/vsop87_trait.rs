// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # VSOP87 Trait and Planet Dispatch
//!
//! ## Scientific scope
//!
//! This module exposes the `VSOP87` object-safe trait, which provides both the
//! heliocentric (VSOP87A) and barycentric (VSOP87E) ecliptic rectangular
//! coordinates for every supported planet.  The coordinates are in the mean
//! ecliptic / equinox of J2000.0 and expressed in astronomical units.
//!
//! ## Technical scope
//!
//! - [`VSOP87`] — trait with two methods:
//!   - `vsop87a(jd: JulianDate) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>`
//!   - `vsop87e(jd: JulianDate) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>`
//!
//! Blanket `impl_vsop87_for_planet!` macro instantiates the trait for
//! `Mercury`, `Venus`, `Earth`, `Mars`, `Jupiter`, `Saturn`, `Uranus`,
//! and `Neptune`, delegating to each planet's inherent `vsop87a` / `vsop87e`
//! methods generated in `vsop87a.rs` and `vsop87e.rs`.
//!
//! The `jd` argument is a `JulianDate` in Terrestrial Time (TT); the
//! implementation converts to TDB via `JulianDate::tt_to_tdb` before computing
//! the series argument T (Julian millennia).
//!
//! ## References
//!
//! - Bretagnon, P., & Francou, G. (1988). "Planetary theories in rectangular
//!   and spherical variables: VSOP87 solutions".
//!   *Astronomy and Astrophysics* 202, 309–315.
use crate::bodies::solar_system::*;
use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::qtty::AstronomicalUnit;
use crate::time::JulianDate;

/// Object-safe trait for planets that support VSOP87 ephemeris computation.
pub trait VSOP87 {
    /// Heliocentric ecliptic rectangular coordinates (VSOP87A solution).
    fn vsop87a(
        &self,
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;
    /// Barycentric ecliptic rectangular coordinates (VSOP87E solution).
    fn vsop87e(&self, jd: JulianDate)
        -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;
}

macro_rules! impl_vsop87_for_planet {
    ($planet:ident) => {
        impl VSOP87 for $planet {
            fn vsop87a(
                &self,
                jd: JulianDate,
            ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
                $planet::vsop87a(jd)
            }

            fn vsop87e(
                &self,
                jd: JulianDate,
            ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
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

    const PRECISION: f64 = 1.0e-12;

    macro_rules! test_dispatch {
        ($planet:ident) => {{
            let body = $planet;
            let dyn_ref: &dyn VSOP87 = &body;
            let jd = crate::J2000;

            let via_trait_a = dyn_ref.vsop87a(jd);
            let via_trait_e = dyn_ref.vsop87e(jd);
            let inherent_a = $planet::vsop87a(jd);
            let inherent_e = $planet::vsop87e(jd);

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
