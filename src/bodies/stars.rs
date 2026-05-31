// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Stars
//!
//! ## Scientific scope
//!
//! Stars are self-luminous plasma spheroids held together by their own gravity
//! and powered by nuclear fusion. From an observational standpoint, a catalog
//! star is characterised by its angular position (right ascension, declination),
//! distance, proper motion, radial velocity, and photometric/physical properties
//! (mass, radius, luminosity). This module targets the Hipparcos/Gaia regime:
//! equatorial coordinates in the geocentric mean equator and equinox of J2000.0
//! (ICRS-aligned), distances in light-years, and proper motions in degrees/year.
//!
//! When both trigonometric parallax and radial velocity are available, the full
//! 6D space-motion propagation (Hipparcos-style, honouring secular acceleration)
//! is applied by [`Star::position_at`]; otherwise a linear transverse-only
//! approximation is used.
//!
//! ## Technical scope
//!
//! - [`Star`] — typed catalog star with name, distance ([`LightYears`]),
//!   mass ([`SolarMasses`]), radius ([`SolarRadiuses`]),
//!   luminosity ([`SolarLuminosities`]), ICRS coordinate with optional
//!   proper motion, parallax ([`MilliArcseconds`]), and radial velocity.
//! - [`Star::position_at`] — propagates position to an arbitrary epoch
//!   ([`JulianDate`]) using proper motion and, when available, the full
//!   Hipparcos-style 6D space-motion model.
//! - [`Star::has_full_space_motion`] — `true` when parallax, radial
//!   velocity, and proper motion are all available.
//!
//! Pre-built star constants live in [`crate::bodies::catalog`].
//!
//! ## References
//!
//! - Hipparcos Catalogue, ESA SP-1200 (1997).
//!   <https://www.cosmos.esa.int/web/hipparcos/catalogues>
//! - van Leeuwen, F. (2007). *Hipparcos, the New Reduction of the Raw Data*.
//!   Astrophysics and Space Science Library 350. doi:10.1007/978-1-4020-6342-8
//! - Luri, X., et al. (2018). "Gaia Data Release 2: using Gaia parallaxes".
//!   *A&A* 616, A9. doi:10.1051/0004-6361/201832964

use crate::astro::proper_motion::{
    propagate_space_motion_since_j2000, ProperMotionError, StarSpaceMotion,
};
use crate::coordinates::spherical::direction;
use crate::coordinates::spherical::position;
use crate::coordinates::{centers::Geocentric, frames::EquatorialMeanJ2000, spherical::Position};
use crate::qtty::length::nominal::SolarRadiuses;
use crate::qtty::velocity::Velocity;
use crate::qtty::*;
use crate::targets::CoordinateWithPM;
use crate::time::JulianDate;

use std::borrow::Cow;

/// Represents a **Star** characterized by its distance, mass, radius, luminosity and position in the sky.
#[derive(Clone, Debug)]
pub struct Star<'a> {
    /// Human-readable name or designation (e.g. `"Alpha Centauri A"`).
    pub name: Cow<'a, str>,
    /// Distance from Earth.
    pub distance: LightYears,
    /// Total mass.
    pub mass: SolarMasses,
    /// Mean photospheric radius.
    pub radius: SolarRadiuses,
    /// Total bolometric luminosity.
    pub luminosity: SolarLuminosities,
    /// J2000 ICRS coordinate, with optional proper motion.
    pub coordinate: CoordinateWithPM<Position<Geocentric, EquatorialMeanJ2000, LightYear>>,
    /// Annual trigonometric parallax.
    ///
    /// When set, callers should apply the parallax ellipse correction against
    /// the observer's position in the Solar System. Required for sub-arcsecond
    /// accuracy on stars closer than ~200 pc.
    pub parallax: Option<MilliArcseconds>,
    /// Radial (line-of-sight) velocity.
    ///
    /// Used to apply the secular acceleration to proper motion: as the star
    /// approaches or recedes its angular PM changes even if the transverse
    /// velocity is constant.
    pub radial_velocity: Option<Velocity<Kilometer, unit::Second>>,
}

impl<'a> Star<'a> {
    /// Compile‐time constructor: only works with `'static` string literals.
    pub const fn new_const(
        name: &'static str,
        distance: LightYears,
        mass: SolarMasses,
        radius: SolarRadiuses,
        luminosity: SolarLuminosities,
        coordinate: CoordinateWithPM<Position<Geocentric, EquatorialMeanJ2000, LightYear>>,
    ) -> Star<'static> {
        Star {
            name: Cow::Borrowed(name),
            distance,
            mass,
            radius,
            luminosity,
            coordinate,
            parallax: None,
            radial_velocity: None,
        }
    }

    /// Runtime constructor: accepts any string-like thing.
    pub fn new<N>(
        name: N,
        distance: LightYears,
        mass: SolarMasses,
        radius: SolarRadiuses,
        luminosity: SolarLuminosities,
        coordinate: CoordinateWithPM<Position<Geocentric, EquatorialMeanJ2000, LightYear>>,
    ) -> Star<'a>
    where
        N: Into<Cow<'a, str>>,
    {
        Star {
            name: name.into(),
            distance,
            mass,
            radius,
            luminosity,
            coordinate,
            parallax: None,
            radial_velocity: None,
        }
    }

    /// Builder: attach an annual trigonometric parallax.
    pub const fn with_parallax(mut self, parallax: MilliArcseconds) -> Self {
        self.parallax = Some(parallax);
        self
    }

    /// Builder: attach a radial velocity (positive = receding).
    pub const fn with_radial_velocity(mut self, rv: Velocity<Kilometer, unit::Second>) -> Self {
        self.radial_velocity = Some(rv);
        self
    }

    /// Compute the BCRS J2000 catalog position at the requested epoch.
    ///
    /// When the catalog provides parallax **and** radial velocity, this
    /// performs a full 6D space-motion propagation using
    /// [`propagate_space_motion_since_j2000`], honouring the parallax-aware
    /// secular acceleration of the proper motion.
    ///
    /// When either parallax or radial velocity is missing, it falls back to
    /// the transverse-only linear propagation in
    /// [`set_proper_motion_since_j2000`](crate::astro::proper_motion::set_proper_motion_since_j2000),
    /// which is appropriate for the vast majority of catalog stars over
    /// human time horizons but can drift several mas/century for fast
    /// nearby objects.
    ///
    /// # Errors
    /// Propagates [`ProperMotionError`] from the underlying routine, which
    /// can occur when:
    /// - `µα⋆` is used at/near the celestial poles, or
    /// - the parallax is non-positive (a noisy fit, not a physical value).
    pub fn position_at(
        &self,
        jd: JulianDate,
    ) -> Result<position::EquatorialMeanJ2000<AstronomicalUnit>, ProperMotionError> {
        // Project the catalog position to AU. Any radial coordinate works:
        // when full space motion is applied, distance is reset from parallax;
        // when proper-motion-only is used, distance is preserved on the AU
        // axis to keep the typed signature happy.
        let pos_ly = self.coordinate.get_position();
        let r_au = pos_ly.distance.to::<AstronomicalUnit>();
        let pos_au = position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            pos_ly.ra(),
            pos_ly.dec(),
            r_au.value(),
        );

        // Full Hipparcos-style propagation when we have both parallax and
        // radial velocity. Otherwise fall back to transverse-only PM.
        match (
            self.parallax,
            self.radial_velocity,
            &self.coordinate.proper_motion,
        ) {
            (Some(parallax), Some(rv), Some(pm)) => {
                // ProperMotion stores deg/yr, but the full-motion API needs
                // mas/yr in the catalog μα⋆ convention. Convert via μα → μα⋆.
                let mas_per_yr_factor = 3_600_000.0_f64; // deg → mas
                let pm_dec_mas = pm.pm_dec.value() * mas_per_yr_factor;
                let pm_ra_mas = match pm.ra_convention {
                    crate::astro::proper_motion::RaProperMotionConvention::MuAlphaStar => {
                        pm.pm_ra.value() * mas_per_yr_factor
                    }
                    crate::astro::proper_motion::RaProperMotionConvention::MuAlpha => {
                        let cos_dec = pos_au.dec().cos();
                        pm.pm_ra.value() * mas_per_yr_factor * cos_dec
                    }
                };
                let motion = StarSpaceMotion {
                    pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate::new(pm_ra_mas),
                    pm_dec: crate::qtty::angular_rate::AngularRate::new(pm_dec_mas),
                    parallax,
                    radial_velocity: rv,
                };
                propagate_space_motion_since_j2000(pos_au, motion, jd)
            }
            (_, _, Some(pm)) => {
                crate::astro::proper_motion::set_proper_motion_since_j2000(pos_au, pm.clone(), jd)
            }
            (_, _, None) => Ok(pos_au),
        }
    }

    /// Returns whether the catalog entry carries enough information for a
    /// full 6D space-motion propagation (parallax + radial velocity + PM).
    pub fn has_full_space_motion(&self) -> bool {
        self.parallax.is_some()
            && self.radial_velocity.is_some()
            && self.coordinate.proper_motion.is_some()
    }
}

impl From<&Star<'_>> for direction::ICRS {
    /// Extracts the J2000 RA/Dec from a [`Star`]'s coordinate position.
    ///
    /// The position's *azimuth* is RA and *polar* is Dec in the
    /// `EquatorialMeanJ2000` frame convention used throughout the crate.
    fn from(star: &Star<'_>) -> Self {
        let pos = star.coordinate.get_position();
        Self::new(pos.azimuth, pos.polar)
    }
}

impl crate::targets::Trackable for Star<'_> {
    type Coords = direction::ICRS;

    #[inline]
    fn track(&self, _jd: JulianDate) -> Self::Coords {
        direction::ICRS::from(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::proper_motion::{ProperMotion, RaProperMotionConvention};
    use crate::qtty::angular::Degrees;
    use crate::targets::CoordinateWithPM;
    use crate::targets::Trackable;
    use crate::J2000;

    type DegreesPerYear =
        crate::qtty::Quantity<crate::qtty::Per<crate::qtty::unit::Degree, crate::qtty::unit::Year>>;

    fn make_pos() -> position::EquatorialMeanJ2000<LightYear> {
        // new_unchecked(polar=dec, azimuth=ra, distance)
        position::EquatorialMeanJ2000::<LightYear>::new_unchecked(
            Degrees::new(5.7),
            Degrees::new(45.0),
            LightYears::new(4.37),
        )
    }

    fn simple_coord() -> CoordinateWithPM<position::EquatorialMeanJ2000<LightYear>> {
        CoordinateWithPM::new_static(make_pos(), J2000)
    }

    fn simple_star() -> Star<'static> {
        Star::new(
            "TestStar",
            LightYears::new(4.37),
            SolarMasses::new(1.1),
            SolarRadiuses::new(1.2),
            SolarLuminosities::new(1.5),
            simple_coord(),
        )
    }

    #[test]
    fn star_name_is_set_correctly() {
        let s = simple_star();
        assert_eq!(s.name, "TestStar");
    }

    #[test]
    fn with_parallax_builder_attaches_value() {
        let s = simple_star().with_parallax(MilliArcseconds::new(742.0));
        assert!(s.parallax.is_some());
        assert!((s.parallax.unwrap().value() - 742.0).abs() < 1e-9);
    }

    #[test]
    fn with_radial_velocity_builder_attaches_value() {
        let rv = Velocity::<Kilometer, unit::Second>::new(-22.0);
        let s = simple_star().with_radial_velocity(rv);
        assert!(s.radial_velocity.is_some());
    }

    #[test]
    fn has_full_space_motion_false_without_parallax_rv_pm() {
        let s = simple_star();
        assert!(!s.has_full_space_motion());
    }

    #[test]
    fn has_full_space_motion_true_with_all_fields() {
        let pm = ProperMotion {
            pm_ra: DegreesPerYear::new(1e-6),
            pm_dec: DegreesPerYear::new(1e-6),
            ra_convention: RaProperMotionConvention::MuAlphaStar,
        };
        let coord = CoordinateWithPM::new(make_pos(), J2000, pm);
        let s = Star::new(
            "FullMotionStar",
            LightYears::new(4.37),
            SolarMasses::new(1.0),
            SolarRadiuses::new(1.0),
            SolarLuminosities::new(1.0),
            coord,
        )
        .with_parallax(MilliArcseconds::new(742.0))
        .with_radial_velocity(Velocity::<Kilometer, unit::Second>::new(-22.0));
        assert!(s.has_full_space_motion());
    }

    #[test]
    fn from_star_ref_extracts_direction() {
        let s = simple_star();
        let dir = direction::ICRS::from(&s);
        let pos = s.coordinate.get_position();
        assert!((dir.azimuth.value() - pos.azimuth.value()).abs() < 1e-9);
        assert!((dir.polar.value() - pos.polar.value()).abs() < 1e-9);
    }

    #[test]
    fn track_matches_from_ref() {
        let s = simple_star();
        let tracked = s.track(J2000);
        let from_ref = direction::ICRS::from(&s);
        assert_eq!(tracked.azimuth, from_ref.azimuth);
        assert_eq!(tracked.polar, from_ref.polar);
    }

    #[test]
    fn position_at_no_pm_returns_ok() {
        let s = simple_star();
        assert!(s.position_at(J2000).is_ok());
    }

    #[test]
    fn position_at_with_pm_only_returns_ok() {
        let pm = ProperMotion {
            pm_ra: DegreesPerYear::new(1e-6),
            pm_dec: DegreesPerYear::new(1e-6),
            ra_convention: RaProperMotionConvention::MuAlphaStar,
        };
        let coord = CoordinateWithPM::new(make_pos(), J2000, pm);
        let s = Star::new(
            "PmOnlyStar",
            LightYears::new(4.37),
            SolarMasses::new(1.0),
            SolarRadiuses::new(1.0),
            SolarLuminosities::new(1.0),
            coord,
        );
        assert!(s.position_at(J2000).is_ok());
    }
}
