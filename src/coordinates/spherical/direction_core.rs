//! Core Direction type and generic implementations.

use super::position_core::Position;
use crate::coordinates::{centers, frames};
use qtty::{Degrees, LengthUnit, Quantity};
use std::ops::Deref;

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
    pub(super) inner: affn::spherical::Direction<F>,
}

impl<F: frames::ReferenceFrame> Direction<F> {
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
        self.inner.position::<C, U>(distance).into()
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
        self.inner
            .position_with_params::<C, U>(center_params, distance)
            .into()
    }

    /// Converts to Cartesian direction.
    #[inline]
    pub fn to_cartesian(&self) -> crate::coordinates::cartesian::Direction<F> {
        self.inner.to_cartesian()
    }

    /// Constructs from a Cartesian direction.
    #[inline]
    pub fn from_cartesian(cart: &crate::coordinates::cartesian::Direction<F>) -> Self {
        affn::spherical::Direction::from_cartesian(cart).into()
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
