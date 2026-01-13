//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! This module provides siderust-owned wrapper types around `affn` with
//! frame-specific inherent constructors following IAU conventions.
//!
//! # Usage
//!
//! No trait imports needed—just use the types directly:
//!
//! ```rust
//! use siderust::coordinates::spherical::direction;
//! use qtty::*;
//!
//! // Inherent constructor with astronomical argument order
//! let dir = direction::ICRS::new(120.0 * DEG, 45.0 * DEG);  // (ra, dec)
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! assert_eq!(dir.dec(), 45.0 * DEG);
//! ```
//!
//! Or use generic types directly:
//!
//! ```rust
//! use siderust::coordinates::spherical::Direction;
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use qtty::*;
//!
//! let dir = Direction::<EquatorialMeanJ2000>::new(120.0 * DEG, 45.0 * DEG);
//! ```

use crate::coordinates::{centers, frames};
use std::ops::Deref;

use qtty::{Degrees, LengthUnit, Quantity};

// =============================================================================
// Helper: Angle Canonicalization
// =============================================================================

/// Normalizes azimuth to [0°, 360°).
#[inline]
fn normalize_azimuth(az: Degrees) -> Degrees {
    az.normalize()
}

/// Clamps polar angle to [-90°, +90°] with proper wrap handling.
#[inline]
fn clamp_polar(polar: Degrees) -> Degrees {
    polar.wrap_quarter_fold()
}

// =============================================================================
// Direction Wrapper Type
// =============================================================================

/// A spherical direction (frame-only, no center or distance).
///
/// This is a siderust-owned wrapper around `affn::spherical::Direction<F>`,
/// providing frame-specific inherent constructors with astronomical conventions.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `EquatorialMeanJ2000`, `Ecliptic`, `Horizontal`)
///
/// # Angular Conventions
///
/// The interpretation of angles depends on the frame:
///
/// | Frame       | First Arg         | Second Arg      | Constructor       |
/// |-------------|-------------------|-----------------|-------------------|
/// | ICRS        | Right Ascension α | Declination δ   | `new(ra, dec)`    |
/// | Equatorial* | Right Ascension α | Declination δ   | `new(ra, dec)`    |
/// | Ecliptic    | Longitude λ       | Latitude β      | `new(lon, lat)`   |
/// | Horizontal  | Altitude Alt      | Azimuth Az      | `new(alt, az)`    |
/// | Geographic  | Longitude λ       | Latitude φ      | `new(lon, lat)`   |
///
/// *Equatorial* refers to `EquatorialMeanJ2000`, `EquatorialMeanOfDate`, or
/// `EquatorialTrueOfDate`, which share the same angular convention.
#[derive(Debug, Clone, Copy)]
#[repr(transparent)]
pub struct Direction<F: frames::ReferenceFrame> {
    inner: affn::spherical::Direction<F>,
}

impl<F: frames::ReferenceFrame> Direction<F> {
    /// Creates a direction from the underlying `affn` type.
    #[inline]
    pub const fn from_inner(inner: affn::spherical::Direction<F>) -> Self {
        Self { inner }
    }

    /// Returns the underlying `affn` direction.
    #[inline]
    pub const fn into_inner(self) -> affn::spherical::Direction<F> {
        self.inner
    }

    /// Returns a reference to the underlying `affn` direction.
    #[inline]
    pub const fn as_inner(&self) -> &affn::spherical::Direction<F> {
        &self.inner
    }

    /// Returns the polar angle (latitude, declination, or altitude) in degrees.
    #[inline]
    pub fn polar(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuthal angle (longitude, right ascension, or azimuth) in degrees.
    #[inline]
    pub fn azimuth(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Promotes this direction to a position with the given distance.
    #[inline]
    pub fn position<C, U>(self, distance: Quantity<U>) -> Position<C, F, U>
    where
        C: centers::ReferenceCenter<Params = ()>,
        U: LengthUnit,
    {
        Position::from_inner(self.inner.position::<C, U>(distance))
    }

    /// Promotes this direction to a position with the given distance and center parameters.
    #[inline]
    pub fn position_with_params<C, U>(
        self,
        center_params: C::Params,
        distance: Quantity<U>,
    ) -> Position<C, F, U>
    where
        C: centers::ReferenceCenter,
        U: LengthUnit,
    {
        Position::from_inner(
            self.inner
                .position_with_params::<C, U>(center_params, distance),
        )
    }

    /// Converts to Cartesian direction.
    #[inline]
    pub fn to_cartesian(&self) -> crate::coordinates::cartesian::Direction<F> {
        self.inner.to_cartesian()
    }

    /// Constructs from a Cartesian direction.
    #[inline]
    pub fn from_cartesian(cart: &crate::coordinates::cartesian::Direction<F>) -> Self {
        Self::from_inner(affn::spherical::Direction::from_cartesian(cart))
    }
}

impl<F: frames::ReferenceFrame> From<affn::spherical::Direction<F>> for Direction<F> {
    #[inline]
    fn from(inner: affn::spherical::Direction<F>) -> Self {
        Self { inner }
    }
}

impl<F: frames::ReferenceFrame> From<Direction<F>> for affn::spherical::Direction<F> {
    #[inline]
    fn from(wrapper: Direction<F>) -> Self {
        wrapper.inner
    }
}

impl<F: frames::ReferenceFrame> Deref for Direction<F> {
    type Target = affn::spherical::Direction<F>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<F: frames::ReferenceFrame> std::fmt::Display for Direction<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.inner)
    }
}

// =============================================================================
// Position Wrapper Type
// =============================================================================

/// A spherical position (center + frame + distance).
///
/// This is a siderust-owned wrapper around `affn::spherical::Position<C, F, U>`,
/// providing frame-specific inherent constructors with astronomical conventions.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`, `Geocentric`)
/// - `F`: The reference frame (e.g., `ICRS`, `EquatorialMeanJ2000`, `Ecliptic`)
/// - `U`: The length unit (e.g., `AstronomicalUnit`, `Kilometer`)
#[derive(Debug, Clone)]
#[repr(transparent)]
pub struct Position<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: LengthUnit> {
    inner: affn::spherical::Position<C, F, U>,
}

impl<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: LengthUnit> Position<C, F, U> {
    /// Creates a position from the underlying `affn` type.
    #[inline]
    pub const fn from_inner(inner: affn::spherical::Position<C, F, U>) -> Self {
        Self { inner }
    }

    /// Returns the underlying `affn` position.
    #[inline]
    pub fn into_inner(self) -> affn::spherical::Position<C, F, U> {
        self.inner
    }

    /// Returns a reference to the underlying `affn` position.
    #[inline]
    pub const fn as_inner(&self) -> &affn::spherical::Position<C, F, U> {
        &self.inner
    }

    /// Returns the polar angle (latitude, declination, or altitude) in degrees.
    #[inline]
    pub fn polar(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuthal angle (longitude, right ascension, or azimuth) in degrees.
    #[inline]
    pub fn azimuth(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the radial distance.
    #[inline]
    pub fn distance(&self) -> Quantity<U> {
        self.inner.distance
    }

    /// Returns a reference to the center parameters.
    #[inline]
    pub fn center_params(&self) -> &C::Params {
        self.inner.center_params()
    }

    /// Extracts the direction (discarding distance).
    #[inline]
    pub fn direction(&self) -> Direction<F> {
        Direction::from_inner(self.inner.direction())
    }

    /// Converts to Cartesian position.
    #[inline]
    pub fn to_cartesian(&self) -> crate::coordinates::cartesian::Position<C, F, U> {
        self.inner.to_cartesian()
    }

    /// Constructs from a Cartesian position.
    #[inline]
    pub fn from_cartesian(cart: &crate::coordinates::cartesian::Position<C, F, U>) -> Self {
        Self::from_inner(affn::spherical::Position::from_cartesian(cart))
    }

    /// Calculates the angular separation between this position and another.
    #[inline]
    pub fn angular_separation(&self, other: Self) -> Degrees {
        self.inner.angular_separation(other.inner)
    }

    /// Euclidean distance to another position in the same center & frame.
    #[inline]
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        self.inner.distance_to(&other.inner)
    }

    /// Creates a position with raw (unchecked) coordinates.
    ///
    /// This is a const-compatible constructor for compile-time catalog entries.
    /// **No normalization is applied** — the caller must ensure coordinates are valid.
    ///
    /// # Arguments
    /// - `polar`: The polar angle (declination/latitude) — should be in [-90°, +90°]
    /// - `azimuth`: The azimuthal angle (RA/longitude) — should be in [0°, 360°)
    /// - `distance`: The radial distance
    #[inline]
    pub const fn new_raw(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self
    where
        C: centers::ReferenceCenter<Params = ()>,
    {
        Self::from_inner(affn::spherical::Position::new_raw(polar, azimuth, distance))
    }
}

impl<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: LengthUnit>
    From<affn::spherical::Position<C, F, U>> for Position<C, F, U>
{
    #[inline]
    fn from(inner: affn::spherical::Position<C, F, U>) -> Self {
        Self { inner }
    }
}

impl<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: LengthUnit> From<Position<C, F, U>>
    for affn::spherical::Position<C, F, U>
{
    #[inline]
    fn from(wrapper: Position<C, F, U>) -> Self {
        wrapper.inner
    }
}

impl<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: LengthUnit> Deref
    for Position<C, F, U>
{
    type Target = affn::spherical::Position<C, F, U>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<C, F, U> std::fmt::Display for Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.inner)
    }
}

// =============================================================================
// Frame-Specific Direction Constructors: ICRS
// =============================================================================

impl Direction<frames::ICRS> {
    /// Creates an ICRS direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Direction Constructors: Equatorial
// =============================================================================

impl Direction<frames::EquatorialMeanJ2000> {
    /// Creates a mean-J2000 equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl Direction<frames::EquatorialMeanOfDate> {
    /// Creates a mean-of-date equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl Direction<frames::EquatorialTrueOfDate> {
    /// Creates a true-of-date equatorial direction from Right Ascension and Declination.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), eastward from vernal equinox, normalized to [0°, 360°)
    /// - `dec`: Declination (δ), from celestial equator, clamped to [-90°, +90°]
    #[inline]
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(dec),
            normalize_azimuth(ra),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Direction Constructors: Ecliptic
// =============================================================================

impl Direction<frames::Ecliptic> {
    /// Creates an Ecliptic direction from ecliptic longitude and latitude.
    ///
    /// # Arguments
    /// - `lon`: Ecliptic longitude (λ), normalized to [0°, 360°)
    /// - `lat`: Ecliptic latitude (β), clamped to [-90°, +90°]
    #[inline]
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(lat),
            normalize_azimuth(lon),
        ))
    }

    /// Returns the ecliptic longitude (λ) in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the ecliptic latitude (β) in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Direction Constructors: Horizontal
// =============================================================================

impl Direction<frames::Horizontal> {
    /// Creates a Horizontal direction from altitude and azimuth (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `alt`: Altitude (elevation) above the horizon, clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    #[inline]
    pub fn new(alt: Degrees, az: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(alt),
            normalize_azimuth(az),
        ))
    }

    /// Returns the altitude (elevation) in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth (from North through East) in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

// =============================================================================
// Frame-Specific Direction Constructors: Geographic (ECEF)
// =============================================================================

impl Direction<frames::ECEF> {
    /// Creates a Geographic direction from longitude and latitude.
    ///
    /// # Arguments
    /// - `lon`: Geodetic longitude, positive eastward, normalized to [0°, 360°)
    /// - `lat`: Geodetic latitude, positive northward, clamped to [-90°, +90°]
    #[inline]
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::from_inner(affn::spherical::Direction::new(
            clamp_polar(lat),
            normalize_azimuth(lon),
        ))
    }

    /// Returns the geodetic longitude in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the geodetic latitude in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Position Constructors: ICRS (for centers with Params = ())
// =============================================================================

impl<C, U> Position<C, frames::ICRS, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates an ICRS position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Position Constructors: Equatorial (for centers with Params = ())
// =============================================================================

impl<C, U> Position<C, frames::EquatorialMeanJ2000, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a mean-J2000 equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl<C, U> Position<C, frames::EquatorialMeanOfDate, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a mean-of-date equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

impl<C, U> Position<C, frames::EquatorialTrueOfDate, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a true-of-date equatorial position from Right Ascension, Declination, and distance.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), normalized to [0°, 360°)
    /// - `dec`: Declination (δ), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(ra: Degrees, dec: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(dec),
            normalize_azimuth(ra),
            distance.into(),
        ))
    }

    /// Returns the Right Ascension (α) in degrees.
    #[inline]
    pub fn ra(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the Declination (δ) in degrees.
    #[inline]
    pub fn dec(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Position Constructors: Ecliptic (for centers with Params = ())
// =============================================================================

impl<C, U> Position<C, frames::Ecliptic, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates an Ecliptic position from longitude, latitude, and distance.
    ///
    /// # Arguments
    /// - `lon`: Ecliptic longitude (λ), normalized to [0°, 360°)
    /// - `lat`: Ecliptic latitude (β), clamped to [-90°, +90°]
    /// - `distance`: Radial distance from the center
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(lon: Degrees, lat: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(lat),
            normalize_azimuth(lon),
            distance.into(),
        ))
    }

    /// Returns the ecliptic longitude (λ) in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the ecliptic latitude (β) in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Frame-Specific Position Constructors: Horizontal (for centers with Params = ())
// =============================================================================

impl<C, U> Position<C, frames::Horizontal, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a Horizontal position from altitude, azimuth, and distance (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `alt`: Altitude (elevation), clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    /// - `distance`: Radial distance from the observer
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(alt: Degrees, az: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(alt),
            normalize_azimuth(az),
            distance.into(),
        ))
    }

    /// Returns the altitude in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

// =============================================================================
// Frame-Specific Position Accessors: Horizontal with Topocentric center
// =============================================================================

impl<U: LengthUnit> Position<centers::Topocentric, frames::Horizontal, U> {
    /// Creates a Topocentric Horizontal position with observer site parameters (IAU Alt-Az convention).
    ///
    /// # Arguments
    /// - `site`: The observer site parameters
    /// - `alt`: Altitude (elevation), clamped to [-90°, +90°]
    /// - `az`: Azimuth, measured from North through East, normalized to [0°, 360°)
    /// - `distance`: Radial distance from the observer
    #[inline]
    pub fn new_with_site<T: Into<Quantity<U>>>(
        site: centers::ObserverSite,
        alt: Degrees,
        az: Degrees,
        distance: T,
    ) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw_with_params(
            site,
            clamp_polar(alt),
            normalize_azimuth(az),
            distance.into(),
        ))
    }

    /// Returns the altitude in degrees.
    #[inline]
    pub fn alt(&self) -> Degrees {
        self.inner.polar
    }

    /// Returns the azimuth in degrees.
    #[inline]
    pub fn az(&self) -> Degrees {
        self.inner.azimuth
    }
}

// =============================================================================
// Frame-Specific Position Constructors: Geographic/ECEF (for centers with Params = ())
// =============================================================================

impl<C, U> Position<C, frames::ECEF, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    U: LengthUnit,
{
    /// Creates a Geographic position from longitude, latitude, and altitude.
    ///
    /// # Arguments
    /// - `lon`: Geodetic longitude, positive eastward, normalized to [0°, 360°)
    /// - `lat`: Geodetic latitude, positive northward, clamped to [-90°, +90°]
    /// - `distance`: Distance from the center (e.g., Earth radius + altitude)
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(lon: Degrees, lat: Degrees, distance: T) -> Self {
        Self::from_inner(affn::spherical::Position::new_raw(
            clamp_polar(lat),
            normalize_azimuth(lon),
            distance.into(),
        ))
    }

    /// Returns the geodetic longitude in degrees.
    #[inline]
    pub fn lon(&self) -> Degrees {
        self.inner.azimuth
    }

    /// Returns the geodetic latitude in degrees.
    #[inline]
    pub fn lat(&self) -> Degrees {
        self.inner.polar
    }
}

// =============================================================================
// Direction type aliases (frame-only, no center)
// =============================================================================

pub mod direction {
    //! Frame-specific direction type aliases.
    //!
    //! These provide convenient shorthand for common direction types.

    use super::frames;
    pub use super::Direction;

    /// **Ecliptic** direction (longitude *λ*, latitude *β*).
    pub type Ecliptic = Direction<frames::Ecliptic>;
    /// **Equatorial mean J2000** direction (right‑ascension *α*, declination *δ*).
    pub type EquatorialMeanJ2000 = Direction<frames::EquatorialMeanJ2000>;
    /// **Equatorial mean of date** direction (right‑ascension *α*, declination *δ*).
    pub type EquatorialMeanOfDate = Direction<frames::EquatorialMeanOfDate>;
    /// **Equatorial true of date** direction (right‑ascension *α*, declination *δ*).
    pub type EquatorialTrueOfDate = Direction<frames::EquatorialTrueOfDate>;
    /// **Horizontal** direction (altitude *Alt*, azimuth *Az*).
    pub type Horizontal = Direction<frames::Horizontal>;
    /// **ICRS** direction.
    pub type ICRS = Direction<frames::ICRS>;
    /// **Geographic** (ECEF) direction: longitude, latitude.
    pub type Geographic = Direction<frames::ECEF>;
}

// =============================================================================
// Position type aliases (center + frame + distance)
// =============================================================================

pub mod position {
    //! Frame and center-specific position type aliases.
    //!
    //! These provide convenient shorthand for common position types.

    pub use super::Position;
    use super::{centers, frames};
    use qtty::Kilometer;

    /// **Heliocentric Ecliptic** coordinates *(λ, β, R)*.
    ///
    /// * `λ` – ecliptic longitude, degrees in `[0, 360)`
    /// * `β` – ecliptic latitude,  degrees in `[-90, 90]`
    /// * `R` – heliocentric distance in unit `U` (e.g. `AstronomicalUnit`)
    pub type Ecliptic<U> = Position<centers::Heliocentric, frames::Ecliptic, U>;

    /// **Geocentric Equatorial mean J2000** coordinates *(α, δ, d)*.
    ///
    /// * `α` – right‑ascension, degrees in `[0, 360)`
    /// * `δ` – declination, degrees in `[-90, 90]`
    /// * `d` – geocentric distance in unit `U` (e.g. `Kilometer`)
    pub type EquatorialMeanJ2000<U> =
        Position<centers::Geocentric, frames::EquatorialMeanJ2000, U>;

    /// **Geocentric Equatorial mean of date** coordinates *(α, δ, d)*.
    pub type EquatorialMeanOfDate<U> =
        Position<centers::Geocentric, frames::EquatorialMeanOfDate, U>;

    /// **Geocentric Equatorial true of date** coordinates *(α, δ, d)*.
    pub type EquatorialTrueOfDate<U> =
        Position<centers::Geocentric, frames::EquatorialTrueOfDate, U>;

    /// **Topocentric Horizontal** coordinates *(Alt, Az, d)*.
    ///
    /// * `Alt` – altitude above the horizon, degrees in `[-90, 90]`
    /// * `Az`  – azimuth from the north, degrees in `[0, 360)`
    /// * `d`   – straight‑line distance from the observer in unit `U`
    pub type Horizontal<U> = Position<centers::Topocentric, frames::Horizontal, U>;

    /// **Barycentric ICRS** coordinates.
    pub type ICRS<U> = Position<centers::Barycentric, frames::ICRS, U>;
    /// **Heliocentric ICRS** coordinates.
    pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
    /// **Geocentric ICRS** coordinates.
    pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;

    /// **Geographic (ECEF)** position: longitude, latitude, altitude (km).
    pub type Geographic = Position<centers::Geocentric, frames::ECEF, Kilometer>;
}
