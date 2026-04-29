// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Star type and related functionality.
//!
//! Represents stars with name, distance, mass, radius, luminosity, and target.
//! - `name`: Name of the star (borrowed or owned).
//! - `distance`: LengthUnit from Earth in light-years (`LightYear`).
//! - `mass`: Stellar mass in solar masses (`SolarMasses`).
//! - `radius`: Stellar radius in solar radii (`SolarRadiuses`).
//! - `luminosity`: Stellar luminosity in solar luminosities (`SolarLuminosity`).
//! - `coordinate`: [`CoordinateWithPM`] pointing to a Spherical coordinates (see [`Position`]), using degrees and Julian Day.
//! - `parallax`: Optional annual parallax in milliarcseconds.
//! - `radial_velocity`: Optional radial velocity in km/s.

use crate::astro::proper_motion::{
    propagate_space_motion_since_j2000, ProperMotionError, StarSpaceMotion,
};
use crate::coordinates::spherical::direction;
use crate::coordinates::spherical::position;
use crate::coordinates::{centers::Geocentric, frames::EquatorialMeanJ2000, spherical::Position};
use crate::targets::CoordinateWithPM;
use crate::time::JulianDate;
use crate::qtty::length::nominal::SolarRadiuses;
use crate::qtty::velocity::Velocity;
use crate::qtty::*;

use std::borrow::Cow;

/// Represents a **Star** characterized by its distance, mass, radius, luminosity and position in the sky.
#[derive(Clone, Debug)]
pub struct Star<'a> {
    pub name: Cow<'a, str>,
    pub distance: LightYears,
    pub mass: SolarMasses,
    pub radius: SolarRadiuses,
    pub luminosity: SolarLuminosities,
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
    pub radial_velocity: Option<Velocity<Kilometer, Second>>,
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
    pub const fn with_radial_velocity(mut self, rv: Velocity<Kilometer, Second>) -> Self {
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
        match (self.parallax, self.radial_velocity, &self.coordinate.proper_motion) {
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
                        let cos_dec = pos_au.dec().value().to_radians().cos();
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
            (_, _, Some(pm)) => crate::astro::proper_motion::set_proper_motion_since_j2000(
                pos_au,
                pm.clone(),
                jd,
            ),
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
