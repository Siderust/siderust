// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Keplerian Orbit Model
//!
//! This module defines the `KeplerianOrbit` struct, which encapsulates the **six classical Keplerian orbital elements**
//! used to describe the motion of a celestial object around a central body, such as a planet around the Sun.
//!
//! These elements are:
//!
//! 1. **Semi-major axis (`a`)**
//!    - Defines the size of the orbit.
//!    - It is half the longest diameter of the ellipse.
//!    - Expressed in astronomical units (AstronomicalUnits).
//!
//! 2. **Eccentricity (`e`)**
//!    - Defines the shape of the orbit.
//!    - `KeplerianOrbit` only supports **elliptic** motion (`0 ≤ e < 1`).
//!    - For parabolic or hyperbolic trajectories, use
//!      [`ConicOrbit`](crate::astro::conic::ConicOrbit) instead.
//!    - Values:
//!       - `e = 0`: circular
//!       - `0 < e < 1`: elliptical
//!
//! 3. **Inclination (`i`)**
//!    - The angle between the orbital plane and a reference plane (typically the ecliptic).
//!    - Describes the tilt of the orbit relative to the reference frame.
//!    - Expressed in degrees.
//!
//! 4. **Longitude of the ascending node (`Ω`)**
//!    - Angle from a fixed reference direction (e.g., the vernal equinox) to the ascending node —
//!      the point where the orbit crosses the reference plane going north.
//!    - Expressed in degrees.
//!
//! 5. **Argument of periapsis (`ω`)**
//!    - The angle from the ascending node to the periapsis (the point of closest approach).
//!    - Measured in the direction of motion.
//!    - Expressed in degrees.
//!
//! 6. **Mean anomaly at epoch (`M₀`)**
//!    - Represents the position of the object along its orbit at a specific reference time (epoch).
//!    - It evolves linearly over time and is used to compute the true anomaly.
//!    - Expressed in degrees.
//!
//! ## Epoch
//!
//! The `epoch` is the reference point in time (given in Julian Day) at which the `mean_anomaly_at_epoch` applies.
//! From this point, the object's position can be propagated using Kepler's equation.
//!
//! ## Coordinate Calculation
//!
//! The `KeplerianOrbit::kepler_position(jd)` method (implemented in
//! [`calculus::kepler_equations`](crate::calculus::kepler_equations)) returns the
//! **heliocentric ecliptic Cartesian coordinates** of the orbiting body at a given
//! Julian Day (`jd`), based on the orbital elements and epoch.
//!
//! ## Units
//!
//! This module assumes that:
//! - **Angles** are expressed in degrees (`Degrees`).
//! - **Distances** use astronomical units (`AstronomicalUnits`).
//! - **Time** is expressed as Julian Days (`JulianDate`).
//!
//! ## Usage Example
//!
//! This example computes Earth's position on a given Julian date.
//!
//! ```rust
//! use siderust::astro::orbit::KeplerianOrbit;
//! use siderust::time::JulianDate;
//! use qtty::*;
//!
//! let earth_orbit = KeplerianOrbit::new(
//!     1.0*AU,                    // a
//!     0.0167,                    // e
//!     Degrees::new(0.00005),     // i
//!     Degrees::new(-11.26064),   // Ω
//!     Degrees::new(102.94719),   // ω
//!     Degrees::new(100.46435),   // M₀
//!     JulianDate::J2000,         // epoch (J2000)
//! );
//!
//! let coords = earth_orbit.kepler_position(JulianDate::new(2459200.5));
//! ```

use crate::time::JulianDate;
use affn::conic::{ConicOrientation, ConicShape, SemiMajorAxisParam};
use affn::frames::EclipticMeanJ2000;
use qtty::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Keplerian orbital elements in semi-major-axis + mean-anomaly-at-epoch form.
///
/// This is the canonical orbital element set for elliptic Solar System bodies
/// (`0 ≤ e < 1`). For parabolic or hyperbolic trajectories, use
/// [`ConicOrbit`](crate::astro::conic::ConicOrbit) instead.
///
/// The shape and orientation are composed from `affn` conic types, tagged
/// to the `EclipticMeanJ2000` frame matching JPL/NASA catalog conventions.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct KeplerianOrbit<U: LengthUnit = AstronomicalUnit> {
    /// Shape: semi-major axis and eccentricity.
    pub shape: SemiMajorAxisParam<U>,
    /// 3-D orientation in the ecliptic mean J2000 frame.
    pub orientation: ConicOrientation<EclipticMeanJ2000>,
    /// Mean anomaly at `epoch`.
    pub mean_anomaly_at_epoch: Degrees,
    /// Reference epoch.
    pub epoch: JulianDate,
}

impl<U: LengthUnit> std::fmt::Display for KeplerianOrbit<U>
where
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "a={}, e={:.6}, i={}, \u{03a9}={}, \u{03c9}={}, M\u{2080}={}, epoch={}",
            self.shape.semi_major_axis,
            self.shape.eccentricity,
            self.orientation.inclination,
            self.orientation.longitude_of_ascending_node,
            self.orientation.argument_of_periapsis,
            self.mean_anomaly_at_epoch,
            self.epoch,
        )
    }
}

impl<U: LengthUnit> KeplerianOrbit<U> {
    /// Creates a new set of Keplerian orbital elements.
    ///
    /// Note: this constructor does not validate. Call [`validate()`](Self::validate)
    /// before propagation to ensure the orbit is elliptic.
    pub const fn new(
        semi_major_axis: Quantity<U>,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Self {
        Self {
            shape: SemiMajorAxisParam::new(semi_major_axis, eccentricity),
            orientation: ConicOrientation::new(
                inclination,
                longitude_of_ascending_node,
                argument_of_periapsis,
            ),
            mean_anomaly_at_epoch,
            epoch,
        }
    }

    /// Validates that this orbit is elliptic (`0 ≤ e < 1`) with a positive
    /// semi-major axis and finite orientation angles.
    pub fn validate(&self) -> Result<(), crate::astro::conic::ConicError> {
        use crate::astro::conic::{map_validation_error, ConicError};
        use affn::conic::ConicKind;

        self.shape.validate().map_err(map_validation_error)?;
        self.orientation.validate().map_err(map_validation_error)?;

        if !matches!(
            self.shape.kind().map_err(map_validation_error)?,
            ConicKind::Elliptic
        ) {
            return Err(ConicError::InvalidEccentricity);
        }
        Ok(())
    }
}
