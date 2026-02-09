// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Coordinate Extension Traits
//!
//! This module provides extension traits that add ergonomic transformation
//! methods to `affn` coordinate types. These traits enable method-chaining
//! style transformations with compile-time type safety.
//!
//! ## Design
//!
//! The extension traits wrap the provider-based transformation system with
//! a clean API:
//!
//! - `to_frame::<F2>(jd, ctx)` - Rotate to a new reference frame.
//! - `to_center::<C2>(jd, ctx)` - Translate to a new reference center.
//! - `to::<C2, F2>(jd, ctx)` - Combined center and frame transformation.
//!
//! ## Transformation Order
//!
//! For combined transformations (`to`), the order is:
//! 1. **Center first** (in source frame): Translate the position.
//! 2. **Then frame**: Rotate to the target frame.
//!
//! This order is chosen because:
//! - Center shifts depend on body positions which are frame-dependent.
//! - Shifting in the source frame before rotating is more intuitive.
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::transform::ext::PositionAstroExt;
//! use siderust::coordinates::transform::context::AstroContext;
//! use siderust::coordinates::cartesian::Position;
//! use siderust::coordinates::centers::{Barycentric, Geocentric};
//! use siderust::coordinates::frames::{Ecliptic, ICRS};
//! use siderust::time::JulianDate;
//! use qtty::AstronomicalUnit;
//!
//! let pos = Position::<Barycentric, Ecliptic, AstronomicalUnit>::new(1.0, 0.5, 0.2);
//! let ctx = AstroContext::new();
//! let jd = JulianDate::J2000;
//!
//! // Transform to Geocentric ICRS
//! let geo_icrs: Position<Geocentric, ICRS, AstronomicalUnit> = pos.to(&jd, &ctx);
//! ```

use crate::coordinates::cartesian::{Direction, Position, Vector};
use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::{CenterShiftProvider, FrameRotationProvider};
use crate::time::JulianDate;
use affn::Rotation3;
use qtty::{LengthUnit, Quantity, Unit};

// =============================================================================
// DirectionAstroExt - Extension trait for Direction<F>
// =============================================================================

/// Extension trait for `Direction<F>` providing frame transformations.
///
/// Directions are unit vectors (translation-invariant), so only frame
/// rotations apply. Center transformations are not meaningful for directions.
pub trait DirectionAstroExt<F: ReferenceFrame> {
    /// Rotates this direction to a new reference frame.
    ///
    /// # Type Parameters
    ///
    /// - `F2`: The target reference frame.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date for time-dependent rotations.
    /// - `ctx`: The astronomical context with model configuration.
    ///
    /// # Returns
    ///
    /// A new `Direction<F2>` representing the same physical direction
    /// expressed in frame `F2`.
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate, ctx: &AstroContext) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>;
}

impl<F: ReferenceFrame> DirectionAstroExt<F> for Direction<F> {
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate, ctx: &AstroContext) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>,
    {
        let rot: Rotation3 = <() as FrameRotationProvider<F, F2>>::rotation(*jd, ctx);
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        // The result is still normalized (rotations preserve length)
        Direction::new_unchecked(x, y, z)
    }
}

// =============================================================================
// VectorAstroExt - Extension trait for Vector<F, U>
// =============================================================================

/// Extension trait for `Vector<F, U>` providing frame transformations.
///
/// Vectors (displacements, velocities) are free vectors that are
/// translation-invariant. Only frame rotations apply.
pub trait VectorAstroExt<F: ReferenceFrame, U: Unit> {
    /// Rotates this vector to a new reference frame.
    ///
    /// # Type Parameters
    ///
    /// - `F2`: The target reference frame.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date for time-dependent rotations.
    /// - `ctx`: The astronomical context with model configuration.
    ///
    /// # Returns
    ///
    /// A new `Vector<F2, U>` representing the same physical vector
    /// expressed in frame `F2`.
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate, ctx: &AstroContext) -> Vector<F2, U>
    where
        (): FrameRotationProvider<F, F2>;
}

impl<F: ReferenceFrame, U: Unit> VectorAstroExt<F, U> for Vector<F, U> {
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate, ctx: &AstroContext) -> Vector<F2, U>
    where
        (): FrameRotationProvider<F, F2>,
    {
        let rot: Rotation3 = <() as FrameRotationProvider<F, F2>>::rotation(*jd, ctx);
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Vector::new(x, y, z)
    }
}

// =============================================================================
// PositionAstroExt - Extension trait for Position<C, F, U>
// =============================================================================

/// Extension trait for `Position<C, F, U>` providing coordinate transformations.
///
/// Positions are affine points that can undergo both frame rotations and
/// center translations.
pub trait PositionAstroExt<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Rotates this position to a new reference frame.
    ///
    /// The center remains unchanged; only the frame orientation changes.
    ///
    /// # Type Parameters
    ///
    /// - `F2`: The target reference frame.
    fn to_frame<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>;

    /// Translates this position to a new reference center.
    ///
    /// The frame remains unchanged; only the origin changes.
    ///
    /// # Type Parameters
    ///
    /// - `C2`: The target reference center.
    fn to_center<C2: ReferenceCenter<Params = ()>>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F, U>
    where
        (): CenterShiftProvider<C, C2, F>;

    /// Transforms this position to a new center and frame.
    ///
    /// # Transformation Order
    ///
    /// 1. Center shift (in source frame F): `C → C2`
    /// 2. Frame rotation: `F → F2`
    ///
    /// This order is documented and consistent. The inverse operation would
    /// apply the inverse rotation first, then the inverse translation.
    ///
    /// # Type Parameters
    ///
    /// - `C2`: The target reference center.
    /// - `F2`: The target reference frame.
    fn to<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>;
}

impl<C, F, U> PositionAstroExt<C, F, U> for Position<C, F, U>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    U: LengthUnit,
{
    fn to_frame<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>,
    {
        let rot: Rotation3 = <() as FrameRotationProvider<F, F2>>::rotation(*jd, ctx);
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Position::new(x, y, z)
    }

    fn to_center<C2: ReferenceCenter<Params = ()>>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F, U>
    where
        (): CenterShiftProvider<C, C2, F>,
    {
        let shift = <() as CenterShiftProvider<C, C2, F>>::shift(*jd, ctx);

        // The shift is in AU; we need to convert if U is different.
        // For now, assume AU and let the type system handle it.
        // TODO: Add unit conversion if U != AstronomicalUnit
        let shift_x = Quantity::<U>::new(shift[0]);
        let shift_y = Quantity::<U>::new(shift[1]);
        let shift_z = Quantity::<U>::new(shift[2]);

        Position::new(self.x() + shift_x, self.y() + shift_y, self.z() + shift_z)
    }

    fn to<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>,
    {
        // Order: center first (in source frame), then rotate
        self.to_center::<C2>(jd, ctx).to_frame::<F2>(jd, ctx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::{Barycentric, Geocentric};
    use crate::coordinates::frames::{Ecliptic, ICRS};
    use qtty::AstronomicalUnit;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_direction_frame_transform() {
        let dir = Direction::<ICRS>::new(1.0, 0.0, 0.0);
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        // ICRS to Ecliptic includes a small frame-bias; don't assume exact axis invariance.
        let dir_ecl: Direction<Ecliptic> = dir.to_frame(&jd, &ctx);

        // Must be finite and length-preserving.
        assert!(dir_ecl.x().is_finite() && dir_ecl.y().is_finite() && dir_ecl.z().is_finite());
        let n0 = (dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z()).sqrt();
        let n1 =
            (dir_ecl.x() * dir_ecl.x() + dir_ecl.y() * dir_ecl.y() + dir_ecl.z() * dir_ecl.z())
                .sqrt();
        assert!((n0 - n1).abs() < 1e-12);
    }

    #[test]
    fn test_direction_frame_roundtrip() {
        let dir = Direction::<ICRS>::new(1.0, 2.0, 3.0);
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let dir_ecl: Direction<Ecliptic> = dir.to_frame(&jd, &ctx);
        let dir_back: Direction<ICRS> = dir_ecl.to_frame(&jd, &ctx);

        assert!((dir_back.x() - dir.x()).abs() < EPSILON);
        assert!((dir_back.y() - dir.y()).abs() < EPSILON);
        assert!((dir_back.z() - dir.z()).abs() < EPSILON);
    }

    #[test]
    fn test_position_frame_transform() {
        let pos = Position::<Barycentric, ICRS, AstronomicalUnit>::new(1.0, 0.0, 0.0);
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let pos_ecl: Position<Barycentric, Ecliptic, AstronomicalUnit> = pos.to_frame(&jd, &ctx);

        assert!(pos_ecl.x().is_finite() && pos_ecl.y().is_finite() && pos_ecl.z().is_finite());

        // Length must be preserved under pure rotation.
        let n0 = (pos.x() * pos.x() + pos.y() * pos.y() + pos.z() * pos.z()).sqrt();
        let n1 =
            (pos_ecl.x() * pos_ecl.x() + pos_ecl.y() * pos_ecl.y() + pos_ecl.z() * pos_ecl.z())
                .sqrt();
        assert!((n0 - n1).abs() < 1e-12);
    }

    #[test]
    fn test_position_center_transform() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        // A point at the Geocentric origin should map to Earth's position in Barycentric
        let geo_origin = Position::<Geocentric, Ecliptic, AstronomicalUnit>::new(0.0, 0.0, 0.0);
        let bary: Position<Barycentric, Ecliptic, AstronomicalUnit> =
            geo_origin.to_center(&jd, &ctx);

        // Should be non-zero (Earth is ~1 AU from barycenter)
        let dist =
            (bary.x().value().powi(2) + bary.y().value().powi(2) + bary.z().value().powi(2)).sqrt();
        assert!(
            dist > 0.9 && dist < 1.1,
            "Earth should be ~1 AU from barycenter, got {}",
            dist
        );
    }

    #[test]
    fn test_position_combined_transform() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let pos = Position::<Barycentric, Ecliptic, AstronomicalUnit>::new(1.0, 0.5, 0.2);

        // Combined transform: Barycentric Ecliptic -> Geocentric ICRS
        let result: Position<Geocentric, ICRS, AstronomicalUnit> = pos.to(&jd, &ctx);

        // Verify it's not the same as the input (transformation happened)
        assert!(
            (result.x().value() - pos.x().value()).abs() > EPSILON
                || (result.y().value() - pos.y().value()).abs() > EPSILON
                || (result.z().value() - pos.z().value()).abs() > EPSILON
        );
    }

    #[test]
    fn test_position_identity_transforms() {
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let pos = Position::<Barycentric, ICRS, AstronomicalUnit>::new(1.5, 2.5, 3.5);

        // Identity frame transform
        let same_frame: Position<Barycentric, ICRS, AstronomicalUnit> = pos.to_frame(&jd, &ctx);
        assert!((same_frame.x().value() - pos.x().value()).abs() < EPSILON);
        assert!((same_frame.y().value() - pos.y().value()).abs() < EPSILON);
        assert!((same_frame.z().value() - pos.z().value()).abs() < EPSILON);

        // Identity center transform
        let same_center: Position<Barycentric, ICRS, AstronomicalUnit> = pos.to_center(&jd, &ctx);
        assert!((same_center.x().value() - pos.x().value()).abs() < EPSILON);
        assert!((same_center.y().value() - pos.y().value()).abs() < EPSILON);
        assert!((same_center.z().value() - pos.z().value()).abs() < EPSILON);
    }
}
