//! Defines spherical coordinate types in the horizontal frame, using altitude (Alt), azimuth (Az), and distance (d).
//!
//! These coordinates are commonly used in observational astronomy to describe the position of celestial objects
//! relative to an observer's local horizon.
//!
//! # Coordinate Convention
//! The `Position<Center, Horizontal>` type uses:
//!
//! - **Altitude (Alt)** → `polar`: angle above or below the horizon, in degrees.
//! - **Azimuth (Az)**   → `azimuth`: angle measured clockwise from the north, in degrees.
//! - **Radial distance (d)** → distance from the reference center, typically in astronomical units (AstronomicalUnits).
//!
//! Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
//! Time is represented as Julian Date (JD).
//!
//! # Topocentric Horizontal Coordinates
//!
//! When using `Topocentric` as the center, the coordinate carries the observer's location
//! as [`ObserverSite`] inside the value. This enables transformations to use the site
//! information without external context.
//!
//! ```rust
//! use siderust::coordinates::centers::{Topocentric, ObserverSite};
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::frames::Horizontal;
//! use qtty::*;
//!
//! let greenwich = ObserverSite::new(0.0 * DEG, 51.4769 * DEG, 0.0 * M);
//! let horizontal = Position::<Topocentric, Horizontal, AstronomicalUnit>::with_site(
//!     greenwich, 45.0 * DEG, 180.0 * DEG, 1.0
//! );
//! assert_eq!(*horizontal.center_params(), greenwich);
//! ```
//!
//! # Example (direction)
//! ```rust
//! use siderust::coordinates::spherical::direction::Horizontal;
//! use qtty::*;
//!
//! let horizontal = Horizontal::new(30.0 * DEG, 90.0 * DEG);
//! println!("alt = {}, az = {}", horizontal.alt(), horizontal.az());
//! ```

use super::*;
use crate::coordinates::{centers::*, frames::*};
use qtty::*;

// =============================================================================
// Direction constructors (frame-only)
// =============================================================================

impl direction::Direction<Horizontal> {
    /// Constructs a new horizontal direction with normalized input angles.
    ///
    /// Altitude is normalized to [-90°, 90°], azimuth to [0°, 360°].
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    pub fn new_horizontal(alt: Degrees, az: Degrees) -> Self {
        Self::new(alt.wrap_quarter_fold(), az.normalize())
    }

    /// Returns the Altitude in degrees.
    pub fn alt(&self) -> Degrees {
        self.polar
    }

    /// Returns the Azimuth in degrees.
    pub fn az(&self) -> Degrees {
        self.azimuth
    }
}

// =============================================================================
// Position constructors for centers with Params = ()
// =============================================================================

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> Position<C, Horizontal, U> {
    /// Creates a new horizontal spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `distance`: LengthUnit to the object, typically in astronomical units (AstronomicalUnits).
    ///
    /// # Returns
    /// A new `Position` in the horizontal frame.
    pub const fn new_const(alt: Degrees, az: Degrees, distance: Quantity<U>) -> Self {
        Self::new_raw(alt, az, distance)
    }

    /// Constructs a new horizontal spherical coordinate with normalized input angular.
    ///
    /// Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `distance`: LengthUnit to the object, typically in astronomical units (AstronomicalUnits).
    ///
    /// # Returns
    /// A `Position` in the horizontal frame.
    pub fn new<A, T>(alt: A, az: A, distance: T) -> Self
    where
        A: Into<Degrees>,
        T: Into<Quantity<U>>,
    {
        Self::new_const(
            alt.into().wrap_quarter_fold(),
            az.into().normalize(),
            distance.into(),
        )
    }
}

// =============================================================================
// Topocentric-specific constructors (with ObserverSite)
// =============================================================================

impl<U: LengthUnit> Position<Topocentric, Horizontal, U> {
    /// Creates a new topocentric horizontal position with observer site.
    ///
    /// # Arguments
    /// - `site`: The observer's geographic location.
    /// - `alt`: Altitude, in degrees.
    /// - `az`: Azimuth (from north, clockwise), in degrees.
    /// - `distance`: Distance to the object.
    ///
    /// # Returns
    /// A `Position` in the horizontal frame with embedded observer site.
    pub const fn with_site_const(
        site: ObserverSite,
        alt: Degrees,
        az: Degrees,
        distance: Quantity<U>,
    ) -> Self {
        Self::new_raw_with_params(site, alt, az, distance)
    }

    /// Creates a new topocentric horizontal position with observer site and normalized angles.
    ///
    /// Altitude is normalized to [-90°, 90°], azimuth to [0°, 360°].
    ///
    /// # Arguments
    /// - `site`: The observer's geographic location.
    /// - `alt`: Altitude, in degrees.
    /// - `az`: Azimuth (from north, clockwise), in degrees.
    /// - `distance`: Distance to the object.
    pub fn with_site<A, T>(site: ObserverSite, alt: A, az: A, distance: T) -> Self
    where
        A: Into<Degrees>,
        T: Into<Quantity<U>>,
    {
        Self::with_site_const(
            site,
            alt.into().wrap_quarter_fold(),
            az.into().normalize(),
            distance.into(),
        )
    }
}

// =============================================================================
// Common accessors for all centers
// =============================================================================

impl<C: ReferenceCenter, U: Unit> SphericalCoord<C, Horizontal, U> {
    /// Returns the Altitude (α) in degrees.
    pub fn alt(&self) -> Degrees {
        self.polar
    }

    /// Returns the Azimuth (θ) in degrees.
    pub fn az(&self) -> Degrees {
        self.azimuth
    }
}
