// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Conic Orbit Models
//!
//! ## Scientific scope
//!
//! This module provides orbit model types for the two cases where the
//! plain elliptic six-element [`KeplerianOrbit`](crate::astro::orbit::KeplerianOrbit)
//! is not the right semantic container:
//!
//! - **[`ConicOrbit`]** — a unified conic (elliptic or hyperbolic) expressed
//!   using **periapsis distance** as the shape parameter.  Used for comets,
//!   near-parabolic intruders, and any body whose eccentricity is not
//!   definitively below 1.
//!
//! - **[`MeanMotionOrbit`]** — an elliptic orbit where the **mean daily
//!   motion** is the authoritative element and the epoch corresponds to zero
//!   mean anomaly.  Used for minor planets from the MPC catalog where mean
//!   motion is the directly fitted quantity.
//!
//! The reusable conic geometry (shape + orientation abstractions) lives in
//! `affn::conic`.  This module layers epoch and anomaly semantics on top.
//!
//! ## Technical scope
//!
//! - [`ConicOrbit`] — validated [`OrientedConic<PeriapsisParam<AU>, Ecliptic>`]
//!   plus mean anomaly at epoch and reference epoch.
//! - [`MeanMotionOrbit`] — validated elliptic oriented conic plus typed mean
//!   motion [`AngularRate<Degree, Day>`] and epoch.  The mean motion field
//!   `pub mean_motion: AngularRate<Degree, Day>` replaces the old untyped
//!   `mean_motion_deg_per_day: f64`; callers needing a raw `f64` call
//!   `.mean_motion.value()`.
//! - [`ConicError`] — unified validation error enum for both types.
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann-Bell.
//! - Minor Planet Center. *MPC Orbit (MPCORB) Database*.
//!   <https://minorplanetcenter.net/iau/MPCORB.html>

use crate::qtty::angular_rate::AngularRate;
use crate::qtty::*;
use crate::time::JulianDate;
use affn::conic::{
    ClassifiedSemiMajorAxisParam, ConicOrientation, ConicValidationError, Elliptic, OrientedConic,
    PeriapsisParam, SemiMajorAxisParam, TypedSemiMajorAxisParam,
};
use affn::frames::EclipticMeanJ2000;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use affn::conic::ConicKind;

/// Validation and propagation errors for conic-based orbit models.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ConicError {
    /// Eccentricity must be finite and non-negative.
    InvalidEccentricity,
    /// Semi-major axis must be finite and strictly positive.
    InvalidSemiMajorAxis,
    /// Periapsis distance must be finite and strictly positive.
    InvalidPeriapsisDistance,
    /// Semi-major axis is undefined for parabolic conics (`e == 1`).
    ParabolicSemiMajorAxis,
    /// Orientation angles must be finite.
    InvalidOrientation,
    /// Parabolic orbits are intentionally not supported yet.
    ParabolicUnsupported,
    /// Hyperbolic eccentricity (`e ≥ 1`) is not valid for this orbit type.
    ///
    /// [`KeplerianOrbit`](crate::astro::orbit::KeplerianOrbit) and
    /// [`MeanMotionOrbit`] only support elliptic orbits (`0 ≤ e < 1`).
    /// Use [`ConicOrbit`] for hyperbolic trajectories.
    HyperbolicNotSupported,
    /// Epoch must be finite (not NaN or infinity).
    InvalidEpoch,
    /// Mean anomaly at epoch must be finite.
    InvalidMeanAnomaly,
    /// Mean motion must be finite and positive.
    InvalidMeanMotion,
    /// The hyperbolic anomaly solver failed to converge.
    HyperbolicSolverFailed,
    /// A parameter value is outside its valid range.
    OutOfRange {
        /// Name of the out-of-range field.
        field: &'static str,
        /// The offending value.
        value: f64,
    },
}

impl std::fmt::Display for ConicError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidEccentricity => write!(f, "invalid eccentricity"),
            Self::InvalidSemiMajorAxis => write!(f, "invalid semi-major axis"),
            Self::InvalidPeriapsisDistance => write!(f, "invalid periapsis distance"),
            Self::ParabolicSemiMajorAxis => {
                write!(
                    f,
                    "semi-major axis is undefined for parabolic conics (e == 1)"
                )
            }
            Self::InvalidOrientation => write!(f, "orientation angles must be finite"),
            Self::ParabolicUnsupported => write!(f, "parabolic orbits are not supported"),
            Self::HyperbolicNotSupported => {
                write!(
                    f,
                    "hyperbolic eccentricity (e >= 1) is not supported by this orbit type; \
                     use ConicOrbit instead"
                )
            }
            Self::InvalidEpoch => write!(f, "epoch must be finite"),
            Self::InvalidMeanAnomaly => write!(f, "mean anomaly at epoch must be finite"),
            Self::InvalidMeanMotion => write!(f, "mean motion must be finite and positive"),
            Self::HyperbolicSolverFailed => {
                write!(f, "hyperbolic anomaly solver failed to converge")
            }
            Self::OutOfRange { field, value } => {
                write!(f, "parameter '{field}' has out-of-range value {value}")
            }
        }
    }
}

impl std::error::Error for ConicError {}

pub(crate) fn map_validation_error(error: ConicValidationError) -> ConicError {
    match error {
        ConicValidationError::InvalidEccentricity => ConicError::InvalidEccentricity,
        ConicValidationError::InvalidSemiMajorAxis => ConicError::InvalidSemiMajorAxis,
        ConicValidationError::InvalidPeriapsisDistance => ConicError::InvalidPeriapsisDistance,
        ConicValidationError::ParabolicSemiMajorAxis => ConicError::ParabolicSemiMajorAxis,
        ConicValidationError::InvalidOrientation => ConicError::InvalidOrientation,
        ConicValidationError::OutOfRange { .. } => ConicError::InvalidOrientation,
    }
}

// =============================================================================
// ConicOrbit
// =============================================================================

/// Unified conic elements expressed using periapsis distance.
///
/// The conic geometry is stored as a validated
/// [`OrientedConic<PeriapsisParam<AstronomicalUnit>, EclipticMeanJ2000>`].
/// After construction the geometry is always valid; `kind()` is infallible.
///
/// Supports elliptic and hyperbolic orbits. Parabolic eccentricity (`e == 1`)
/// is rejected at construction time by [`try_new`](Self::try_new).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ConicOrbit {
    geometry: OrientedConic<PeriapsisParam<AstronomicalUnit>, EclipticMeanJ2000>,
    /// Mean anomaly at `epoch`.
    pub mean_anomaly_at_epoch: Degrees,
    /// Reference epoch.
    pub epoch: JulianDate,
}

impl ConicOrbit {
    /// Creates a new validated conic-element set.
    ///
    /// Returns an error if any parameter is invalid (non-positive periapsis
    /// distance, non-finite eccentricity, non-finite orientation angles,
    /// non-finite mean anomaly or epoch, or parabolic eccentricity).
    pub fn try_new(
        periapsis_distance: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Result<Self, ConicError> {
        let shape = PeriapsisParam::try_new(periapsis_distance, eccentricity)
            .map_err(map_validation_error)?;
        if eccentricity == 1.0 {
            return Err(ConicError::ParabolicUnsupported);
        }
        let orientation = ConicOrientation::try_new(
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
        )
        .map_err(map_validation_error)?;
        if !mean_anomaly_at_epoch.is_finite() {
            return Err(ConicError::InvalidMeanAnomaly);
        }
        if !epoch.raw().value().is_finite() {
            return Err(ConicError::InvalidEpoch);
        }
        Ok(Self {
            geometry: OrientedConic::new(shape, orientation),
            mean_anomaly_at_epoch,
            epoch,
        })
    }

    /// Creates a new conic-element set **without validation**.
    ///
    /// Intended for compile-time body constants with known-correct values.
    pub const fn new_unchecked(
        periapsis_distance: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Self {
        Self {
            geometry: OrientedConic::new(
                PeriapsisParam::new_unchecked(periapsis_distance, eccentricity),
                ConicOrientation::new(
                    inclination,
                    longitude_of_ascending_node,
                    argument_of_periapsis,
                ),
            ),
            mean_anomaly_at_epoch,
            epoch,
        }
    }

    /// The validated oriented geometry.
    #[inline]
    pub fn geometry(&self) -> &OrientedConic<PeriapsisParam<AstronomicalUnit>, EclipticMeanJ2000> {
        &self.geometry
    }

    /// Classifies the orbit from its eccentricity. Infallible.
    pub fn kind(&self) -> ConicKind {
        self.geometry.kind()
    }
}

// =============================================================================
// MeanMotionOrbit
// =============================================================================

/// Mean-motion-driven elliptic elements.
///
/// This is intentionally distinct from [`crate::astro::orbit::KeplerianOrbit`].
/// Here the stored mean daily motion is authoritative and the epoch corresponds
/// to zero mean anomaly.
///
/// Stores `OrientedConic<TypedSemiMajorAxisParam<AstronomicalUnit, Elliptic>,
/// EclipticMeanJ2000>` — elliptic geometry is enforced at construction time.
///
/// The mean motion is typed as [`AngularRate<Degree, Day>`] so callers cannot
/// accidentally use it as radians-per-day or confuse it with other rates.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MeanMotionOrbit {
    geometry: OrientedConic<TypedSemiMajorAxisParam<AstronomicalUnit, Elliptic>, EclipticMeanJ2000>,
    /// Mean motion (degrees per day), typed.
    pub mean_motion: AngularRate<Degree, Day>,
    /// Epoch at which the mean anomaly is defined to be zero.
    pub epoch: JulianDate,
}

impl MeanMotionOrbit {
    /// Creates a new validated elliptic mean-motion orbit.
    ///
    /// Returns an error if any parameter is invalid, if the eccentricity is
    /// not elliptic (`e >= 1`), or if the mean motion or epoch is non-finite.
    pub fn try_new(
        semi_major_axis: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_motion: AngularRate<Degree, Day>,
        epoch: JulianDate,
    ) -> Result<Self, ConicError> {
        let sma = SemiMajorAxisParam::try_new(semi_major_axis, eccentricity)
            .map_err(map_validation_error)?;
        let typed = match sma.classify() {
            ClassifiedSemiMajorAxisParam::Elliptic(t) => t,
            ClassifiedSemiMajorAxisParam::Hyperbolic(_) => {
                return Err(ConicError::HyperbolicNotSupported);
            }
        };
        let orientation = ConicOrientation::try_new(
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
        )
        .map_err(map_validation_error)?;
        if !mean_motion.value().is_finite() || mean_motion.value() <= 0.0 {
            return Err(ConicError::InvalidMeanMotion);
        }
        if !epoch.raw().value().is_finite() {
            return Err(ConicError::InvalidEpoch);
        }
        Ok(Self {
            geometry: OrientedConic::new(typed, orientation),
            mean_motion,
            epoch,
        })
    }

    /// Creates a new mean-motion orbit **without validation**.
    ///
    /// Intended for compile-time body constants with known-correct elliptic values.
    pub const fn new_unchecked(
        semi_major_axis: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_motion: AngularRate<Degree, Day>,
        epoch: JulianDate,
    ) -> Self {
        Self {
            geometry: OrientedConic::new(
                TypedSemiMajorAxisParam::new_unchecked(SemiMajorAxisParam::new_unchecked(
                    semi_major_axis,
                    eccentricity,
                )),
                ConicOrientation::new(
                    inclination,
                    longitude_of_ascending_node,
                    argument_of_periapsis,
                ),
            ),
            mean_motion,
            epoch,
        }
    }

    /// The validated oriented geometry.
    #[inline]
    pub fn geometry(
        &self,
    ) -> &OrientedConic<TypedSemiMajorAxisParam<AstronomicalUnit, Elliptic>, EclipticMeanJ2000>
    {
        &self.geometry
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_elliptic() {
        let orbit = ConicOrbit::try_new(
            1.0 * AU,
            0.5,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        )
        .unwrap();
        assert_eq!(orbit.kind(), ConicKind::Elliptic);
    }

    #[test]
    fn classify_hyperbolic() {
        let orbit = ConicOrbit::try_new(
            1.0 * AU,
            1.5,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        )
        .unwrap();
        assert_eq!(orbit.kind(), ConicKind::Hyperbolic);
    }

    #[test]
    fn negative_eccentricity_is_invalid() {
        assert_eq!(
            ConicOrbit::try_new(
                1.0 * AU,
                -0.1,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            ),
            Err(ConicError::InvalidEccentricity)
        );
    }

    #[test]
    fn mean_motion_rejects_hyperbolic() {
        assert_eq!(
            MeanMotionOrbit::try_new(
                1.0 * AU,
                1.1,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(1.0),
                JulianDate::J2000,
            ),
            Err(ConicError::HyperbolicNotSupported)
        );
    }

    #[test]
    fn conic_rejects_nan_epoch() {
        assert_eq!(
            ConicOrbit::try_new(
                1.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::from_raw_unchecked(qtty::Day::new(f64::NAN)),
            ),
            Err(ConicError::InvalidEpoch)
        );
    }

    #[test]
    fn conic_rejects_inf_mean_anomaly() {
        assert_eq!(
            ConicOrbit::try_new(
                1.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(f64::INFINITY),
                JulianDate::J2000,
            ),
            Err(ConicError::InvalidMeanAnomaly)
        );
    }

    #[test]
    fn conic_rejects_parabolic() {
        assert_eq!(
            ConicOrbit::try_new(
                1.0 * AU,
                1.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            ),
            Err(ConicError::ParabolicUnsupported)
        );
    }

    #[test]
    fn mean_motion_rejects_nan_mean_motion() {
        assert_eq!(
            MeanMotionOrbit::try_new(
                1.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(f64::NAN),
                JulianDate::J2000,
            ),
            Err(ConicError::InvalidMeanMotion)
        );
    }

    #[test]
    fn mean_motion_rejects_negative_mean_motion() {
        assert_eq!(
            MeanMotionOrbit::try_new(
                1.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(-1.0),
                JulianDate::J2000,
            ),
            Err(ConicError::InvalidMeanMotion)
        );
    }

    #[test]
    fn mean_motion_rejects_nan_epoch() {
        assert_eq!(
            MeanMotionOrbit::try_new(
                1.0 * AU,
                0.5,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                AngularRate::<Degree, Day>::new(1.0),
                JulianDate::from_raw_unchecked(qtty::Day::new(f64::NAN)),
            ),
            Err(ConicError::InvalidEpoch)
        );
    }
}
