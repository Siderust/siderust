// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Keplerian Orbit Model
//!
//! ## Scientific scope
//!
//! A **Keplerian orbit** is the idealized two-body solution of Newton's law of
//! gravitation.  The shape and orientation of the orbit are fully described by
//! the **six classical orbital elements**:
//!
//! | Element | Symbol | Unit | Meaning |
//! |---------|--------|------|---------|
//! | Semi-major axis | a | AU | Size of the orbit (half the long axis) |
//! | Eccentricity | e | — | Shape: 0 = circle, 0 < e < 1 = ellipse |
//! | Inclination | i | ° | Tilt relative to the ecliptic plane |
//! | Longitude of ascending node | Ω | ° | Where the orbit crosses the ecliptic going north |
//! | Argument of periapsis | ω | ° | Angle from node to periapsis |
//! | Mean anomaly at epoch | M₀ | ° | Angular position at reference epoch |
//!
//! This module restricts itself to **elliptic** orbits (e < 1).  For parabolic
//! and hyperbolic trajectories use [`crate::astro::conic::ConicOrbit`].
//!
//! ## Technical scope
//!
//! This module is the astronomy-facing wrapper over the lower-level
//! [`keplerian`](https://docs.rs/keplerian) crate. Pure Keplerian math such as
//! anomaly solving and two-body propagation lives in `keplerian`; `siderust`
//! keeps [`KeplerianOrbit<U>`] here because it binds that math to astronomy
//! semantics: TT-scale [`JulianDate`] epochs from `tempoch`,
//! [`EclipticMeanJ2000`] orientation from `affn`, and the heliocentric
//! Gaussian/GM conventions used by Solar-System elements.
//!
//! - [`KeplerianOrbit<U>`] — storage type for the six elements, generic over
//!   length unit `U`.  Constructed via the infallible `new` (trusted constants)
//!   or the validated `try_new` (user/deserialized data).
//! - [`PreparedOrbit`] — precomputed version for repeated propagation: caches
//!   mean motion as a typed [`AngularRate<Radian, Day>`], orientation trig, and
//!   M₀ in radians.  Mean motion is exposed through the typed accessor
//!   [`PreparedOrbit::mean_motion`]; callers requiring the raw `f64` should
//!   call `.value()` on the return value.
//! - `OrientationTrig` — cached sin/cos of inclination, Ω, ω (crate-internal).
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann-Bell.
//! - Standish, E. M. (1992). "Keplerian Elements for Approximate Positions of
//!   the Major Planets". *JPL Solar System Dynamics*.
//!   <https://ssd.jpl.nasa.gov/planets/approx_pos.html>

use crate::astro::conic::{elliptic_geometry_from_sma, ConicError};
use crate::astro::units::heliocentric_period_days;
use crate::qtty::angular_rate::AngularRate;
use crate::qtty::*;
use crate::time::JulianDate;
use affn::conic::{ConicOrientation, Elliptic, SemiMajorAxisParam, TypedSemiMajorAxisParam};
use affn::frames::EclipticMeanJ2000;

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
        let geometry = elliptic_geometry_from_sma(
            semi_major_axis,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
        )?;

        if !mean_anomaly_at_epoch.is_finite() {
            return Err(ConicError::InvalidMeanAnomaly);
        }
        if !epoch.raw().value().is_finite() {
            return Err(ConicError::InvalidEpoch);
        }
        Ok(Self {
            shape: *geometry.shape(),
            orientation: *geometry.orientation(),
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
        let (sin_i, cos_i) = o.inclination().sin_cos();
        let (sin_omega, cos_omega) = o.argument_of_periapsis().sin_cos();
        let (sin_node, cos_node) = o.longitude_of_ascending_node().sin_cos();
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
/// - Mean motion as [`AngularRate<Radian, Day>`] — avoids recomputing
///   `k / a^{1.5}` every call and carries the physical unit.
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
///     siderust::J2000,
/// );
/// let prepared = PreparedOrbit::try_from(orbit).unwrap();
/// let pos = prepared.position_at(siderust::JulianDate::new(2459200.5));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedOrbit {
    /// The validated source orbit (kept for introspection / conversion back).
    elements: KeplerianOrbit,
    /// Precomputed mean motion (rad/day), typed.
    mean_motion: AngularRate<Radian, Day>,
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

    /// Precomputed mean motion (radians per day), typed as [`AngularRate<Radian, Day>`].
    ///
    /// Callers that need a raw `f64` for math kernels should call `.value()`.
    #[inline]
    pub fn mean_motion(&self) -> AngularRate<Radian, Day> {
        self.mean_motion
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
        // Kepler's 3rd law in the AU-day system (heliocentric).
        let period_days = heliocentric_period_days(a);
        let mean_motion = AngularRate::<Radian, Day>::new(std::f64::consts::TAU / period_days);
        let m0_rad = orbit.mean_anomaly_at_epoch.to::<Radian>().value();
        let trig = OrientationTrig::from_orientation(orbit.orientation());
        Self {
            elements: orbit,
            mean_motion,
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
