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
//!    - Angle from a fixed reference direction (e.g., the vernal equinox) to the ascending node,
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
//! use siderust::qtty::*;
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

use crate::astro::conic::ConicError;
use crate::astro::units::GaussianYears;
use crate::time::JulianDate;
use affn::conic::{
    ClassifiedSemiMajorAxisParam, ConicOrientation, Elliptic, SemiMajorAxisParam,
    TypedSemiMajorAxisParam,
};
use affn::frames::EclipticMeanJ2000;
use crate::qtty::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Keplerian orbital elements in semi-major-axis + mean-anomaly-at-epoch form.
///
/// This is the canonical orbital element set for elliptic Solar System bodies
/// (`0 ≤ e < 1`). For parabolic or hyperbolic trajectories, use
/// [`ConicOrbit`](crate::astro::conic::ConicOrbit) instead.
///
/// The elliptic constraint is enforced at construction time via `try_new`.
/// The infallible `new` constructor skips validation and is intended for
/// compile-time body constants with known-correct values.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct KeplerianOrbit<U: LengthUnit = AstronomicalUnit> {
    shape: TypedSemiMajorAxisParam<U, Elliptic>,
    orientation: ConicOrientation<EclipticMeanJ2000>,
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
            self.shape.semi_major_axis(),
            self.shape.eccentricity(),
            self.orientation.inclination(),
            self.orientation.longitude_of_ascending_node(),
            self.orientation.argument_of_periapsis(),
            self.mean_anomaly_at_epoch,
            self.epoch,
        )
    }
}

impl<U: LengthUnit> KeplerianOrbit<U> {
    /// Creates a new set of Keplerian orbital elements **without validation**.
    ///
    /// Intended for compile-time body constants where the values are known to
    /// be a valid elliptic orbit. For user-supplied data use
    /// [`try_new`](Self::try_new) instead.
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
            shape: TypedSemiMajorAxisParam::new_unchecked(SemiMajorAxisParam::new_unchecked(
                semi_major_axis,
                eccentricity,
            )),
            orientation: ConicOrientation::new(
                inclination,
                longitude_of_ascending_node,
                argument_of_periapsis,
            ),
            mean_anomaly_at_epoch,
            epoch,
        }
    }

    /// Creates a new set of Keplerian orbital elements, returning an error
    /// if the orbit is not a valid ellipse.
    ///
    /// Prefer this over [`new()`](Self::new) for user-supplied or deserialized
    /// data. The orbit is validated exactly once at construction time.
    pub fn try_new(
        semi_major_axis: Quantity<U>,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Result<Self, crate::astro::conic::ConicError> {
        use crate::astro::conic::map_validation_error;

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

        if !mean_anomaly_at_epoch.value().is_finite() {
            return Err(ConicError::InvalidMeanAnomaly);
        }
        if !epoch.value().is_finite() {
            return Err(ConicError::InvalidEpoch);
        }
        Ok(Self {
            shape: typed,
            orientation,
            mean_anomaly_at_epoch,
            epoch,
        })
    }

    /// The typed elliptic shape (semi-major axis + eccentricity).
    #[inline]
    pub fn shape(&self) -> &TypedSemiMajorAxisParam<U, Elliptic> {
        &self.shape
    }

    /// The 3-D orientation in the ecliptic mean J2000 frame.
    #[inline]
    pub fn orientation(&self) -> &ConicOrientation<EclipticMeanJ2000> {
        &self.orientation
    }
}

// =============================================================================
// PreparedOrbit — validated + precomputed for hot propagation
// =============================================================================

/// Precomputed orientation trig values for efficient repeated rotation to ecliptic.
///
/// These values are computed once from the three orbital orientation angles and
/// reused on every propagation call, avoiding redundant `sin_cos()` evaluations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct OrientationTrig {
    sin_i: f64,
    cos_i: f64,
    sin_omega: f64,
    cos_omega: f64,
    sin_node: f64,
    cos_node: f64,
}

impl OrientationTrig {
    #[inline]
    pub(crate) fn sin_i(&self) -> f64 {
        self.sin_i
    }
    #[inline]
    pub(crate) fn cos_i(&self) -> f64 {
        self.cos_i
    }
    #[inline]
    pub(crate) fn sin_omega(&self) -> f64 {
        self.sin_omega
    }
    #[inline]
    pub(crate) fn cos_omega(&self) -> f64 {
        self.cos_omega
    }
    #[inline]
    pub(crate) fn sin_node(&self) -> f64 {
        self.sin_node
    }
    #[inline]
    pub(crate) fn cos_node(&self) -> f64 {
        self.cos_node
    }

    pub(crate) fn from_orientation(o: &ConicOrientation<EclipticMeanJ2000>) -> Self {
        let (sin_i, cos_i) = o.inclination().to::<Radian>().value().sin_cos();
        let (sin_omega, cos_omega) = o.argument_of_periapsis().to::<Radian>().value().sin_cos();
        let (sin_node, cos_node) = o
            .longitude_of_ascending_node()
            .to::<Radian>()
            .value()
            .sin_cos();
        Self {
            sin_i,
            cos_i,
            sin_omega,
            cos_omega,
            sin_node,
            cos_node,
        }
    }
}

/// A validated, precomputed elliptic orbit optimized for repeated propagation.
///
/// `PreparedOrbit` guarantees at construction time that the underlying
/// [`KeplerianOrbit`] is a valid ellipse. It caches:
///
/// - Mean motion (radians/day) — avoids recomputing `k / a^{1.5}` every call.
/// - Orientation trig terms — avoids 3× `sin_cos()` per call.
/// - Mean anomaly at epoch (radians) — avoids degree→radian conversion per call.
///
/// Use [`TryFrom<KeplerianOrbit>`] or [`PreparedOrbit::try_from_elements`] to
/// construct.
///
/// # Example
///
/// ```rust
/// use siderust::astro::orbit::{KeplerianOrbit, PreparedOrbit};
/// use siderust::time::JulianDate;
/// use siderust::qtty::*;
///
/// let orbit = KeplerianOrbit::new(
///     1.0 * AU, 0.0167,
///     Degrees::new(0.00005), Degrees::new(-11.26064),
///     Degrees::new(102.94719), Degrees::new(100.46435),
///     JulianDate::J2000,
/// );
/// let prepared = PreparedOrbit::try_from(orbit).unwrap();
/// let pos = prepared.position_at(JulianDate::new(2459200.5));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedOrbit {
    /// The validated source orbit (kept for introspection / conversion back).
    elements: KeplerianOrbit,
    /// Precomputed mean motion in radians per day.
    mean_motion_rad_per_day: f64,
    /// Precomputed mean anomaly at epoch in radians.
    m0_rad: f64,
    /// Precomputed orientation trig.
    trig: OrientationTrig,
}

impl PreparedOrbit {
    /// Construct from individual orbital elements, validating and precomputing.
    pub fn try_from_elements(
        semi_major_axis: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_periapsis: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Result<Self, ConicError> {
        let orbit = KeplerianOrbit::try_new(
            semi_major_axis,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
            mean_anomaly_at_epoch,
            epoch,
        )?;
        Ok(Self::from_validated(orbit))
    }

    /// Access the underlying validated [`KeplerianOrbit`].
    #[inline]
    pub fn elements(&self) -> &KeplerianOrbit {
        &self.elements
    }

    /// Precomputed mean motion in radians per day.
    #[inline]
    pub fn mean_motion_rad_per_day(&self) -> f64 {
        self.mean_motion_rad_per_day
    }

    /// Precomputed mean anomaly at epoch in radians.
    #[inline]
    pub(crate) fn m0_rad(&self) -> f64 {
        self.m0_rad
    }

    /// Precomputed orientation trig values.
    #[inline]
    pub(crate) fn orientation_trig(&self) -> &OrientationTrig {
        &self.trig
    }

    /// Internal: build from an already-validated orbit.
    pub(crate) fn from_validated(orbit: KeplerianOrbit) -> Self {
        let a = orbit.shape().semi_major_axis().value();
        // Kepler's 3rd law: T [Gaussian years] = a [AU]^{3/2}
        let t_gaussian_years = a * a.sqrt();
        let period_days = GaussianYears::new(t_gaussian_years).to::<Day>().value();
        let mean_motion_rad_per_day = std::f64::consts::TAU / period_days;
        let m0_rad = orbit.mean_anomaly_at_epoch.to::<Radian>().value();
        let trig = OrientationTrig::from_orientation(orbit.orientation());
        Self {
            elements: orbit,
            mean_motion_rad_per_day,
            m0_rad,
            trig,
        }
    }
}

impl TryFrom<KeplerianOrbit> for PreparedOrbit {
    type Error = ConicError;

    fn try_from(orbit: KeplerianOrbit) -> Result<Self, Self::Error> {
        // KeplerianOrbit::try_new already validates; orbits produced by `new`
        // (const, body constants) are trusted. Build directly.
        Ok(Self::from_validated(orbit))
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for PreparedOrbit {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let orbit = KeplerianOrbit::deserialize(deserializer)?;
        Self::try_from(orbit).map_err(serde::de::Error::custom)
    }
}

#[cfg(feature = "serde")]
impl Serialize for PreparedOrbit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.elements.serialize(serializer)
    }
}
