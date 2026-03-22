// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Additional orbit models used when a plain elliptic Keplerian orbit is not
//! the right semantic container.
//!
//! The reusable conic geometry itself lives in `affn::conic`. This module keeps
//! the astronomy-specific state that layers epoch and anomaly semantics on top
//! of that geometry for siderust propagation code.

use crate::time::JulianDate;
use affn::conic::{
    ConicOrientation, ConicValidationError, OrientedPeriapsisConic, OrientedSemiMajorAxisConic,
};
use qtty::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use affn::conic::ConicKind;

/// Validation and propagation errors for conic-based orbit models.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ConicError {
    /// Eccentricity must be finite and non-negative.
    InvalidEccentricity,
    /// Semi-major axis must be finite and strictly positive.
    InvalidSemiMajorAxis,
    /// Periapsis distance must be finite and strictly positive.
    InvalidPeriapsisDistance,
    /// Parabolic orbits are intentionally not supported yet.
    ParabolicUnsupported,
}

impl std::fmt::Display for ConicError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidEccentricity => write!(f, "invalid eccentricity"),
            Self::InvalidSemiMajorAxis => write!(f, "invalid semi-major axis"),
            Self::InvalidPeriapsisDistance => write!(f, "invalid periapsis distance"),
            Self::ParabolicUnsupported => write!(f, "parabolic orbits are not supported"),
        }
    }
}

impl std::error::Error for ConicError {}

pub(crate) fn map_validation_error(error: ConicValidationError) -> ConicError {
    match error {
        ConicValidationError::InvalidEccentricity => ConicError::InvalidEccentricity,
        ConicValidationError::InvalidSemiMajorAxis => ConicError::InvalidSemiMajorAxis,
        ConicValidationError::InvalidPeriapsisDistance => ConicError::InvalidPeriapsisDistance,
    }
}

/// Unified conic elements expressed using periapsis distance.
///
/// The conic geometry itself is represented by
/// [`affn::conic::OrientedPeriapsisConic`]. `siderust::ConicOrbit` adds the
/// mean-anomaly-at-epoch and epoch fields needed for astronomical propagation.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ConicOrbit {
    /// Periapsis distance.
    pub periapsis_distance: AstronomicalUnits,
    /// Orbital eccentricity.
    pub eccentricity: f64,
    /// Inclination.
    pub inclination: Degrees,
    /// Longitude of the ascending node.
    pub longitude_of_ascending_node: Degrees,
    /// Argument of periapsis.
    pub argument_of_periapsis: Degrees,
    /// Mean anomaly at `epoch`.
    pub mean_anomaly_at_epoch: Degrees,
    /// Reference epoch.
    pub epoch: JulianDate,
}

impl ConicOrbit {
    /// Creates a new conic-element set.
    pub const fn new(
        periapsis_distance: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Self {
        Self {
            periapsis_distance,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
            mean_anomaly_at_epoch,
            epoch,
        }
    }

    /// Classifies the orbit from its eccentricity.
    pub fn kind(&self) -> Result<ConicKind, ConicError> {
        self.geometry().kind().map_err(map_validation_error)
    }

    pub(crate) fn geometry(&self) -> OrientedPeriapsisConic<AstronomicalUnit> {
        OrientedPeriapsisConic::new(
            self.periapsis_distance,
            self.eccentricity,
            ConicOrientation::new(
                self.inclination,
                self.longitude_of_ascending_node,
                self.argument_of_periapsis,
            ),
        )
    }

    pub(crate) fn validated_geometry(
        &self,
    ) -> Result<OrientedPeriapsisConic<AstronomicalUnit>, ConicError> {
        let geometry = self.geometry();
        geometry.validate().map_err(map_validation_error)?;
        Ok(geometry)
    }
}

/// Mean-motion-driven elliptic elements.
///
/// This is intentionally distinct from [`crate::astro::orbit::Orbit`]. Here the
/// stored mean daily motion is authoritative and the epoch corresponds to zero
/// mean anomaly, which matches a number of catalog and legacy data sets. The
/// reusable geometry fields are modeled by
/// [`affn::conic::OrientedSemiMajorAxisConic`].
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MeanMotionOrbit {
    /// Semi-major axis.
    pub semi_major_axis: AstronomicalUnits,
    /// Orbital eccentricity.
    pub eccentricity: f64,
    /// Inclination.
    pub inclination: Degrees,
    /// Longitude of the ascending node.
    pub longitude_of_ascending_node: Degrees,
    /// Argument of periapsis.
    pub argument_of_periapsis: Degrees,
    /// Mean motion in degrees per day.
    pub mean_motion_deg_per_day: f64,
    /// Epoch at which the mean anomaly is defined to be zero.
    pub epoch: JulianDate,
}

impl MeanMotionOrbit {
    /// Creates a new mean-motion orbit.
    pub const fn new(
        semi_major_axis: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_motion_deg_per_day: f64,
        epoch: JulianDate,
    ) -> Self {
        Self {
            semi_major_axis,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
            mean_motion_deg_per_day,
            epoch,
        }
    }

    pub(crate) fn geometry(&self) -> OrientedSemiMajorAxisConic<AstronomicalUnit> {
        OrientedSemiMajorAxisConic::new(
            self.semi_major_axis,
            self.eccentricity,
            ConicOrientation::new(
                self.inclination,
                self.longitude_of_ascending_node,
                self.argument_of_periapsis,
            ),
        )
    }

    pub(crate) fn validated_geometry(
        &self,
    ) -> Result<OrientedSemiMajorAxisConic<AstronomicalUnit>, ConicError> {
        let geometry = self.geometry();
        geometry.validate().map_err(map_validation_error)?;
        if !matches!(
            geometry.kind().map_err(map_validation_error)?,
            ConicKind::Elliptic
        ) {
            return Err(ConicError::InvalidEccentricity);
        }
        Ok(geometry)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_elliptic() {
        let orbit = ConicOrbit::new(
            1.0 * AU,
            0.5,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        assert_eq!(orbit.kind().unwrap(), ConicKind::Elliptic);
    }

    #[test]
    fn classify_hyperbolic() {
        let orbit = ConicOrbit::new(
            1.0 * AU,
            1.5,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        assert_eq!(orbit.kind().unwrap(), ConicKind::Hyperbolic);
    }

    #[test]
    fn negative_eccentricity_is_invalid() {
        assert_eq!(
            ConicOrbit::new(
                1.0 * AU,
                -0.1,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            )
            .kind(),
            Err(ConicError::InvalidEccentricity)
        );
    }
}
