//! Core Position type and generic implementations.

use super::direction_core::Direction;
use crate::coordinates::{centers, frames};
use qtty::{Degrees, LengthUnit, Quantity};
use std::ops::Deref;

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
    pub(super) inner: affn::spherical::Position<C, F, U>,
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
