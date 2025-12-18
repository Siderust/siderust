//! # Cartesian Direction (Free Unit Vectors)
//!
//! This module defines direction types - unit vectors that represent orientations
//! in space without any spatial origin (center).
//!
//! ## Mathematical Foundations
//!
//! Directions are **free vectors**: they live in the vector space, not in affine
//! space. Unlike positions, directions are translation-invariant:
//!
//! - Rotating a direction is valid (frame transformation)
//! - "Translating" a direction to a different origin is **meaningless**
//!
//! Therefore, `Direction<F>` has no center parameter, only a frame `F`.
//!
//! ## Usage
//!
//! To get an observer-dependent direction (line of sight to a target), use
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) to compute
//! it from two positions.
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::cartesian::Direction;
//! use siderust::coordinates::frames::Ecliptic;
//!
//! // Create a normalized direction in the ecliptic frame
//! let dir = Direction::<Ecliptic>::normalize(1.0, 2.0, 2.0);
//! ```

use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::frames;
use qtty::{LengthUnit, Quantity};

use nalgebra::Vector3;
use std::marker::PhantomData;

/// A unit vector representing a direction in 3D space.
///
/// Directions are frame-dependent but center-independent (free vectors).
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
#[derive(Debug, Clone, Copy)]
pub struct Direction<F: frames::ReferenceFrame> {
    xyz: Vector3<Quantity<DirectionUnit>>,
    _frame: PhantomData<F>,
}

impl<F: frames::ReferenceFrame> Direction<F> {
    /// Creates a direction from components (must be pre-normalized or will be normalized).
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        let norm = Vector3::<f64>::new(x, y, z).normalize();
        Self {
            xyz: Vector3::new(
                Quantity::<DirectionUnit>::new(norm.x),
                Quantity::<DirectionUnit>::new(norm.y),
                Quantity::<DirectionUnit>::new(norm.z),
            ),
            _frame: PhantomData,
        }
    }

    /// Creates a direction from components, normalizing to unit length.
    pub fn normalize(x: f64, y: f64, z: f64) -> Self {
        Self::new(x, y, z)
    }

    /// Creates a direction from a nalgebra Vector3 of Quantity<DirectionUnit>.
    pub fn from_vec3(vec3: Vector3<Quantity<DirectionUnit>>) -> Self {
        Self {
            xyz: vec3,
            _frame: PhantomData,
        }
    }

    /// Returns the underlying nalgebra Vector3.
    pub const fn as_vec3(&self) -> Vector3<Quantity<DirectionUnit>> {
        self.xyz
    }

    /// Gets the x-component.
    pub fn x(&self) -> Quantity<DirectionUnit> {
        self.xyz[0]
    }

    /// Gets the y-component.
    pub fn y(&self) -> Quantity<DirectionUnit> {
        self.xyz[1]
    }

    /// Gets the z-component.
    pub fn z(&self) -> Quantity<DirectionUnit> {
        self.xyz[2]
    }

    /// Returns a Position vector in the same direction, scaled by the given magnitude.
    ///
    /// # Type Parameters
    /// - `C`: The center for the resulting position
    /// - `U`: The length unit for the magnitude
    pub fn position<C, U>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U>
    where
        C: crate::coordinates::centers::ReferenceCenter<Params = ()>,
        U: LengthUnit,
    {
        super::Position::new(
            magnitude * self.x().value(),
            magnitude * self.y().value(),
            magnitude * self.z().value(),
        )
    }

    /// Returns a Position vector in the same direction, scaled by the given magnitude,
    /// with explicit center parameters.
    pub fn position_with_params<C, U>(
        &self,
        center_params: C::Params,
        magnitude: Quantity<U>,
    ) -> super::Position<C, F, U>
    where
        C: crate::coordinates::centers::ReferenceCenter,
        U: LengthUnit,
    {
        super::Position::new_with_params(
            center_params,
            magnitude * self.x().value(),
            magnitude * self.y().value(),
            magnitude * self.z().value(),
        )
    }

    /// Converts this cartesian direction to a spherical direction.
    pub fn to_spherical(&self) -> crate::coordinates::spherical::Direction<F> {
        use qtty::Degrees;
        
        let x = self.x().value();
        let y = self.y().value();
        let z = self.z().value();
        let r = (x * x + y * y + z * z).sqrt();
        
        if r == 0.0 {
            return crate::coordinates::spherical::Direction::<F>::new_raw(
                Degrees::new(0.0),
                Degrees::new(0.0),
            );
        }
        
        let polar = Degrees::new((z / r).asin().to_degrees());
        let azimuth = Degrees::new(y.atan2(x).to_degrees());
        
        crate::coordinates::spherical::Direction::<F>::new_raw(polar, azimuth)
    }

    /// Returns a formatted string representation of this direction vector.
    pub fn display(&self) -> String {
        format!(
            "Frame: {}, X: {:.6}, Y: {:.6}, Z: {:.6}",
            F::frame_name(),
            self.x().value(),
            self.y().value(),
            self.z().value()
        )
    }
}

impl<F: frames::ReferenceFrame> std::fmt::Display for Direction<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Frame: {}, X: {:.6}, Y: {:.6}, Z: {:.6}",
            F::frame_name(),
            self.x().value(),
            self.y().value(),
            self.z().value()
        )
    }
}

/// Type aliases for common direction systems.
pub type Ecliptic = Direction<frames::Ecliptic>;
pub type Equatorial = Direction<frames::Equatorial>;
pub type Horizontal = Direction<frames::Horizontal>;
pub type Geographic = Direction<frames::ECEF>;
pub type ICRS = Direction<frames::ICRS>;
