// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Reference Centers Module
//!
//! This module defines astronomical reference centers (origins) for coordinate systems.
//! A reference center specifies the origin point from which positions are measured.
//!
//! ## Architecture
//!
//! All center types implement the [`ReferenceCenter`] trait from `affn`, which provides
//! a common interface. The trait itself is re-exported from `affn` for convenience.
//!
//! ## Predefined Centers
//!
//! The following reference centers are provided:
//!
//! - [`Barycentric`]: Center of mass of the solar system.
//! - [`Heliocentric`]: Center of the Sun.
//! - [`Geocentric`]: Center of the Earth.
//! - [`Topocentric`]: Observer's location on the surface of the Earth (parameterized by [`Geodetic<ECEF>`]).
//! - [`Bodycentric`]: Generic center for any orbiting celestial body (parameterized by [`BodycentricParams`]).
//!
//! ## Compile-Time vs Runtime Safety
//!
//! | Center            | `Params`              | Safety level                        |
//! |-------------------|-----------------------|-------------------------------------|
//! | `Barycentric`     | `()`                  | **Compile-time** — zero cost        |
//! | `Heliocentric`    | `()`                  | **Compile-time** — zero cost        |
//! | `Geocentric`      | `()`                  | **Compile-time** — zero cost        |
//! | `Topocentric`     | [`Geodetic<ECEF>`]    | **Runtime-parameterized** (see below) |
//! | `Bodycentric`     | [`BodycentricParams`] | **Runtime-parameterized** (see below) |
//!
//! For centers with `Params = ()`, the type system guarantees correctness at
//! compile time with no runtime overhead.
//!
//! For **parameterized centers**, the *center type* is still enforced at compile
//! time (you cannot mix `Topocentric` and `Geocentric` positions), but operations
//! that require *matching parameters* — such as `Position - Position` and
//! `distance_to` — verify parameter equality at **runtime**:
//!
//! - **Panicking API** — the default operators (`Sub`, `distance_to`) `assert!`
//!   in all build profiles.
//! - **Checked API** — [`Position::checked_sub`](affn::Position::checked_sub) and
//!   [`Position::try_distance_to`](affn::Position::try_distance_to) return
//!   `Result<_, CenterParamsMismatchError>`.
//!
//! ### Validated Construction
//!
//! Observatory catalog constants (e.g., in [`crate::observatories`]) are of
//! type [`Geodetic<ECEF>`] and can be used directly as `Topocentric` parameters.
//! Use [`Geodetic::<ECEF>::new`] to construct validated geodetic coordinates
//! (longitude and latitude are normalised automatically).
//!
//! ## Extending
//!
//! To define a new reference center, use the derive macro:
//!
//! ```rust
//! use affn::prelude::*;
//!
//! #[derive(Debug, Copy, Clone, ReferenceCenter)]
//! struct Lunarcentric;
//! assert_eq!(Lunarcentric::center_name(), "Lunarcentric");
//! ```
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::centers::{ReferenceCenter, Geocentric};
//!
//! let name = Geocentric::center_name();
//! assert_eq!(name, "Geocentric");
//! ```

use crate::astro::orbit::Orbit;
use qtty::*;
use std::fmt::Debug;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// Re-export core traits from affn
pub use affn::centers::ReferenceCenter;
pub use affn::{AffineCenter, CenterParamsMismatchError, NoCenter};
// Import derives from prelude for use in this module
use affn::prelude::ReferenceCenter as DeriveReferenceCenter;

use super::frames::ECEF;

/// A geodetic position: an ellipsoidal coordinate with a geocentric origin.
///
/// The frame `F` determines which ellipsoid is used via
/// [`HasEllipsoid`](affn::ellipsoid::HasEllipsoid):
///
/// - `Geodetic<ECEF>` — WGS84 ellipsoid
/// - `Geodetic<ITRF>` — GRS80 ellipsoid
pub type Geodetic<F, U = qtty::Meter> = affn::ellipsoidal::Position<Geocentric, F, U>;

// Required for Transform specialization
#[derive(Debug, Copy, Clone, DeriveReferenceCenter)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Heliocentric;

#[derive(Debug, Copy, Clone, DeriveReferenceCenter)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Barycentric;


// =============================================================================
// Topocentric Center (parameterized)
// =============================================================================

/// Observer's location on the surface of the Earth.
///
/// Unlike other reference centers, `Topocentric` is *parameterized*: coordinates
/// with this center carry a [`Geodetic<ECEF>`] that specifies the observer's
/// geographic location.  This allows horizontal coordinates to know their
/// observation site without external context.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::centers::{Topocentric, Geodetic, ReferenceCenter};
/// use siderust::coordinates::frames::ECEF;
/// use qtty::*;
///
/// // Topocentric coordinates require a geodetic site
/// let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
///
/// // The site is stored as Topocentric::Params
/// assert_eq!(std::mem::size_of::<<Topocentric as ReferenceCenter>::Params>(),
///            std::mem::size_of::<Geodetic<ECEF>>());
/// ```
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Topocentric;

impl ReferenceCenter for Topocentric {
    type Params = Geodetic<ECEF>;
    fn center_name() -> &'static str {
        "Topocentric"
    }
}
impl affn::AffineCenter for Topocentric {}

impl Topocentric {
    /// Creates a Topocentric Horizontal position with observer site.
    ///
    /// Angles are canonicalized:
    /// - `alt` is folded to `[-90°, +90°]`
    /// - `az` is normalized to `[0°, 360°)`
    #[inline]
    pub fn horizontal<U: qtty::LengthUnit, T: Into<qtty::Quantity<U>>>(
        site: Geodetic<ECEF>,
        alt: qtty::Degrees,
        az: qtty::Degrees,
        distance: T,
    ) -> affn::spherical::Position<Topocentric, super::frames::Horizontal, U> {
        affn::spherical::Position::new_raw_with_params(
            site,
            alt.wrap_quarter_fold(),
            az.normalize(),
            distance.into(),
        )
    }
}

#[derive(Debug, Copy, Clone, DeriveReferenceCenter)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Geocentric;

// =============================================================================
// Bodycentric: Generic center for any orbiting celestial body
// =============================================================================

/// Specifies the reference center for an orbit (where the orbit is defined relative to).
///
/// When transforming to/from body-centric coordinates, the orbit must be converted
/// to match the coordinate system being transformed. This enum indicates which
/// standard center the orbit is relative to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum OrbitReferenceCenter {
    /// Orbit is defined relative to the solar system barycenter.
    Barycentric,
    /// Orbit is defined relative to the Sun (most common for planets, asteroids, comets).
    #[default]
    Heliocentric,
    /// Orbit is defined relative to Earth (for artificial satellites, the Moon).
    Geocentric,
}

/// Parameters for a body-centered coordinate system.
///
/// This struct specifies the orbital elements of a celestial body and which
/// reference center the orbit is defined relative to. This allows computing
/// the body's position at any Julian date using Keplerian propagation.
///
/// # Fields
///
/// - `orbit`: The Keplerian orbital elements of the body.
/// - `orbit_center`: Which center the orbit is defined relative to.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::centers::{BodycentricParams, OrbitReferenceCenter};
/// use siderust::astro::orbit::Orbit;
/// use siderust::time::JulianDate;
/// use qtty::*;
///
/// // Mars-like orbit (heliocentric)
/// let mars_orbit = Orbit::new(
///     1.524 * AU,           // semi-major axis
///     0.0934,               // eccentricity
///     Degrees::new(1.85),   // inclination
///     Degrees::new(49.56),  // longitude of ascending node
///     Degrees::new(286.5),  // argument of perihelion
///     Degrees::new(19.41),  // mean anomaly at epoch
///     JulianDate::J2000,    // epoch
/// );
///
/// let mars_params = BodycentricParams::heliocentric(mars_orbit);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BodycentricParams {
    /// The Keplerian orbital elements of the body.
    pub orbit: Orbit,
    /// Which standard center the orbit is defined relative to.
    pub orbit_center: OrbitReferenceCenter,
}

impl BodycentricParams {
    /// Creates parameters for a body-centered coordinate system.
    ///
    /// # Arguments
    ///
    /// - `orbit`: The Keplerian orbital elements of the body.
    /// - `orbit_center`: Which center the orbit is defined relative to.
    pub const fn new(orbit: Orbit, orbit_center: OrbitReferenceCenter) -> Self {
        Self {
            orbit,
            orbit_center,
        }
    }

    /// Creates parameters for a body orbiting the Sun (heliocentric orbit).
    ///
    /// This is the most common case for planets, asteroids, and comets.
    pub const fn heliocentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Heliocentric)
    }

    /// Creates parameters for a body orbiting Earth (geocentric orbit).
    ///
    /// Use this for artificial satellites, the Moon, etc.
    pub const fn geocentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Geocentric)
    }

    /// Creates parameters for a body orbiting the solar system barycenter.
    pub const fn barycentric(orbit: Orbit) -> Self {
        Self::new(orbit, OrbitReferenceCenter::Barycentric)
    }
}

impl Default for BodycentricParams {
    /// Returns default parameters with a circular 1 AU heliocentric orbit.
    ///
    /// Note: This default is primarily for internal use. In practice, you should
    /// always provide meaningful orbital elements for body-centric calculations.
    fn default() -> Self {
        use crate::time::JulianDate;
        use qtty::AstronomicalUnits;

        Self {
            orbit: Orbit::new(
                AstronomicalUnits::new(1.0),
                0.0,
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                Degrees::new(0.0),
                JulianDate::J2000,
            ),
            orbit_center: OrbitReferenceCenter::Heliocentric,
        }
    }
}

/// Generic center for any orbiting celestial body.
///
/// This allows defining coordinate systems centered on satellites, planets,
/// moons, comets, asteroids, or any other body with known orbital elements.
/// The body's position at any time is computed via Keplerian propagation.
///
/// # Type Aliases
///
/// For convenience, you can use `Bodycentric` for any case, but conceptually:
/// - Satellites orbiting Earth: `Bodycentric` with `OrbitReferenceCenter::Geocentric`
/// - Planets/asteroids/comets: `Bodycentric` with `OrbitReferenceCenter::Heliocentric`
/// - Moons of other planets: Requires hierarchical orbit handling (future work)
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::centers::{Bodycentric, BodycentricParams, ReferenceCenter};
/// use siderust::coordinates::cartesian::Position;
/// use siderust::coordinates::frames;
/// use siderust::astro::orbit::Orbit;
/// use siderust::time::JulianDate;
/// use qtty::*;
///
/// // Create orbital parameters for an Earth-orbiting satellite
/// let satellite_orbit = Orbit::new(
///     0.0000426 * AU,       // ~6378 km (low Earth orbit) in AU
///     0.001,                // nearly circular
///     Degrees::new(51.6),   // ISS-like inclination
///     Degrees::new(0.0),
///     Degrees::new(0.0),
///     Degrees::new(0.0),
///     JulianDate::J2000,
/// );
///
/// let sat_params = BodycentricParams::geocentric(satellite_orbit);
/// ```
#[derive(Debug, Copy, Clone, DeriveReferenceCenter)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[center(params = BodycentricParams)]
pub struct Bodycentric;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn center_names_are_correct() {
        assert_eq!(Barycentric::center_name(), "Barycentric");
        assert_eq!(Heliocentric::center_name(), "Heliocentric");
        assert_eq!(Topocentric::center_name(), "Topocentric");
        assert_eq!(Geocentric::center_name(), "Geocentric");
        assert_eq!(Bodycentric::center_name(), "Bodycentric");
        assert_eq!(<() as ReferenceCenter>::center_name(), "");
    }

    #[test]
    fn standard_centers_have_unit_params() {
        // Verify that standard centers use () as Params (zero-cost)
        let _: <Barycentric as ReferenceCenter>::Params = ();
        let _: <Heliocentric as ReferenceCenter>::Params = ();
        let _: <Geocentric as ReferenceCenter>::Params = ();
        let _: <() as ReferenceCenter>::Params = ();

        // Verify zero size
        assert_eq!(
            std::mem::size_of::<<Barycentric as ReferenceCenter>::Params>(),
            0
        );
        assert_eq!(
            std::mem::size_of::<<Heliocentric as ReferenceCenter>::Params>(),
            0
        );
        assert_eq!(
            std::mem::size_of::<<Geocentric as ReferenceCenter>::Params>(),
            0
        );
    }

    #[test]
    fn topocentric_has_geodetic_coord_params() {
        // Verify Topocentric uses Geodetic<ECEF> as Params
        let site = Geodetic::<ECEF>::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
        let _: <Topocentric as ReferenceCenter>::Params = site;

        // Verify non-zero size (stores actual data)
        assert!(std::mem::size_of::<<Topocentric as ReferenceCenter>::Params>() > 0);
    }

    #[test]
    fn geodetic_coord_default() {
        let coord = Geodetic::<ECEF>::default();
        assert_eq!(coord.lon, 0.0);
        assert_eq!(coord.lat, 0.0);
        assert_eq!(coord.height, 0.0);
    }

    #[test]
    fn geodetic_coord_equality() {
        let c1 = Geodetic::<ECEF>::new(10.0 * DEG, 20.0 * DEG, 100.0 * M);
        let c2 = Geodetic::<ECEF>::new(10.0 * DEG, 20.0 * DEG, 100.0 * M);
        let c3 = Geodetic::<ECEF>::new(10.0 * DEG, 20.0 * DEG, 200.0 * M);
        assert_eq!(c1, c2);
        assert_ne!(c1, c3);
    }

    #[test]
    fn bodycentric_has_params() {
        use crate::time::JulianDate;

        // Create a simple orbit
        let orbit = Orbit::new(
            1.524 * AU,
            0.0934,
            Degrees::new(1.85),
            Degrees::new(49.56),
            Degrees::new(286.5),
            Degrees::new(19.41),
            JulianDate::J2000,
        );

        let params = BodycentricParams::heliocentric(orbit);
        let _: <Bodycentric as ReferenceCenter>::Params = params;

        // Verify non-zero size (stores actual data)
        assert!(std::mem::size_of::<<Bodycentric as ReferenceCenter>::Params>() > 0);
    }

    #[test]
    fn bodycentric_params_constructors() {
        use crate::time::JulianDate;

        let orbit = Orbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let helio = BodycentricParams::heliocentric(orbit);
        assert_eq!(helio.orbit_center, OrbitReferenceCenter::Heliocentric);

        let geo = BodycentricParams::geocentric(orbit);
        assert_eq!(geo.orbit_center, OrbitReferenceCenter::Geocentric);

        let bary = BodycentricParams::barycentric(orbit);
        assert_eq!(bary.orbit_center, OrbitReferenceCenter::Barycentric);
    }

    #[test]
    fn bodycentric_params_default() {
        let params = BodycentricParams::default();
        assert_eq!(params.orbit_center, OrbitReferenceCenter::Heliocentric);
        assert_eq!(params.orbit.eccentricity, 0.0);
    }

    #[test]
    fn bodycentric_params_equality() {
        use crate::time::JulianDate;

        let orbit1 = Orbit::new(
            1.0 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );
        let orbit2 = Orbit::new(
            2.0 * AU, // Different semi-major axis
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let params1 = BodycentricParams::heliocentric(orbit1);
        let params2 = BodycentricParams::heliocentric(orbit1);
        let params3 = BodycentricParams::heliocentric(orbit2);

        assert_eq!(params1, params2);
        assert_ne!(params1, params3);
    }

    #[test]
    fn center_params_mismatch_error_reexported() {
        // Verify CenterParamsMismatchError is accessible through centers module
        let err = CenterParamsMismatchError { operation: "test" };
        let _: &dyn std::error::Error = &err;
    }
}
