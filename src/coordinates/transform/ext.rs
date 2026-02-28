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
//! The default API uses IAU models with no context argument required:
//!
//! - `to_frame::<F2>(jd_tt)` — Rotate to a new reference frame.
//! - `to::<C2, F2>(jd_tt)` — Combined center and frame transformation.
//!
//! Center-only transforms are exposed on the [`TransformCenter`] trait:
//!
//! - `pos.to_center(params, jd)` — Shift to a new reference center (all variants).
//! - `pos.to_center_with(params, jd, &ctx)` — Same with a custom context.
//!
//! For frame-only expert overrides, a `_with` suffix variant accepts an [`AstroContext`]:
//!
//! - `to_frame_with::<F2>(jd_tt, &ctx)` — Frame rotation with custom context.
//! - `to_with::<C2, F2>(jd_tt, &ctx)` — Combined transform with custom context.
//!
//! Alternatively, wrap a coordinate with a custom context using
//! [`WithEngine`] and use the same method names:
//!
//! ```rust,ignore
//! coord.using(&engine).to_frame::<F2>(&jd);
//! ```
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
//! use siderust::coordinates::cartesian::Position;
//! use siderust::coordinates::centers::{Barycentric, Geocentric};
//! use siderust::coordinates::frames::{EclipticMeanJ2000, ICRS};
//! use siderust::time::JulianDate;
//! use qtty::AstronomicalUnit;
//!
//! let pos = Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(1.0, 0.5, 0.2);
//! let jd = JulianDate::J2000;
//!
//! // Transform to Geocentric ICRS — no context needed
//! let geo_icrs: Position<Geocentric, ICRS, AstronomicalUnit> = pos.to(&jd);
//! ```

use crate::coordinates::cartesian::{Direction, Position, Vector};
use crate::coordinates::centers::{Geodetic, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, ECEF};
use crate::coordinates::spherical;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::{CenterShiftProvider, FrameRotationProvider};
use crate::time::JulianDate;
use affn::Rotation3;
use qtty::{LengthUnit, Unit};

// =============================================================================
// DirectionAstroExt - Extension trait for Direction<F>
// =============================================================================

/// Extension trait for `Direction<F>` providing frame transformations.
///
/// Directions are unit vectors (translation-invariant), so only frame
/// rotations apply. Center transformations are not meaningful for directions.
pub trait DirectionAstroExt<F: ReferenceFrame> {
    /// Rotates this direction to a new reference frame using IAU defaults.
    ///
    /// # Type Parameters
    ///
    /// - `F2`: The target reference frame.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date (TT) for time-dependent rotations.
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>;

    /// Rotates this direction to a new reference frame with a custom context.
    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>;

    /// Converts this direction to ecliptic-of-date coordinates (convenience).
    ///
    /// Available for ICRS and GCRS directions via the provider system.
    fn to_ecliptic_of_date(
        &self,
        jd_tt: &JulianDate,
    ) -> Direction<crate::coordinates::frames::EclipticTrueOfDate>
    where
        Self: crate::coordinates::transform::ecliptic_of_date::ToEclipticTrueOfDate,
    {
        crate::coordinates::transform::ecliptic_of_date::ToEclipticTrueOfDate::to_ecliptic_of_date(
            self, jd_tt,
        )
    }

    /// Converts this direction to horizontal coordinates using TT only.
    ///
    /// UT1 is inferred from TT via the built-in ΔT model.
    fn to_horizontal(
        &self,
        jd_tt: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<crate::coordinates::frames::Horizontal>
    where
        Self: crate::coordinates::transform::horizontal::ToHorizontal,
    {
        let jd_ut1 = crate::astro::earth_rotation::jd_ut1_from_tt(*jd_tt);
        crate::coordinates::transform::horizontal::ToHorizontal::to_horizontal(
            self, &jd_ut1, jd_tt, site,
        )
    }

    /// Converts this direction to horizontal coordinates with explicit UT1+TT.
    ///
    /// Use this when you have a precise UT1 value (e.g. from IERS EOP).
    fn to_horizontal_precise(
        &self,
        jd_tt: &JulianDate,
        jd_ut1: &JulianDate,
        site: &Geodetic<ECEF>,
    ) -> Direction<crate::coordinates::frames::Horizontal>
    where
        Self: crate::coordinates::transform::horizontal::ToHorizontal,
    {
        crate::coordinates::transform::horizontal::ToHorizontal::to_horizontal(
            self, jd_ut1, jd_tt, site,
        )
    }
}

impl<F: ReferenceFrame> DirectionAstroExt<F> for Direction<F> {
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.to_frame_with(jd, &AstroContext::default())
    }

    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Direction<F2>
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
// SphericalDirectionAstroExt - Extension trait for spherical::Direction<F>
// =============================================================================

/// Extension trait for `spherical::Direction<F>` providing time-dependent
/// frame transformations via the provider system.
///
/// This is the spherical counterpart of [`DirectionAstroExt`]. Internally,
/// it converts to a cartesian [`Direction`], applies the rotation, and converts
/// back.
pub trait SphericalDirectionAstroExt<F: ReferenceFrame> {
    /// Rotates this spherical direction to a new reference frame (IAU defaults).
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> spherical::Direction<F2>
    where
        (): FrameRotationProvider<F, F2>;

    /// Rotates this spherical direction to a new reference frame with custom context.
    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> spherical::Direction<F2>
    where
        (): FrameRotationProvider<F, F2>;
}

impl<F: ReferenceFrame> SphericalDirectionAstroExt<F> for spherical::Direction<F> {
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> spherical::Direction<F2>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.to_frame_with(jd, &AstroContext::default())
    }

    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> spherical::Direction<F2>
    where
        (): FrameRotationProvider<F, F2>,
    {
        let cart: Direction<F> = self.to_cartesian();
        let cart_f2: Direction<F2> = cart.to_frame_with(jd, ctx);
        spherical::Direction::from_cartesian(&cart_f2)
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
    /// Rotates this vector to a new reference frame (IAU defaults).
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Vector<F2, U>
    where
        (): FrameRotationProvider<F, F2>;

    /// Rotates this vector to a new reference frame with custom context.
    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Vector<F2, U>
    where
        (): FrameRotationProvider<F, F2>;
}

impl<F: ReferenceFrame, U: Unit> VectorAstroExt<F, U> for Vector<F, U> {
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Vector<F2, U>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.to_frame_with(jd, &AstroContext::default())
    }

    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Vector<F2, U>
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

/// Extension trait for `Position<C, F, U>` providing frame transformations
/// and combined center+frame transformations.
///
/// Positions are affine points that can undergo both frame rotations and
/// center translations.
///
/// **Center-only transforms** are provided by [`TransformCenter`] (use
/// `pos.to_center(params, jd)`).
///
/// Default methods use IAU models with no context argument.
/// `_with` variants accept an [`AstroContext`] for expert overrides.
pub trait PositionAstroExt<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Rotates this position to a new reference frame (IAU defaults).
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>;

    /// Rotates this position to a new reference frame with custom context.
    fn to_frame_with<F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>;

    /// Transforms this position to a new center and frame (IAU defaults).
    ///
    /// # Transformation Order
    ///
    /// 1. Center shift (in source frame F): `C → C2`
    /// 2. Frame rotation: `F → F2`
    fn to<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>;

    /// Transforms this position to a new center and frame with custom context.
    fn to_with<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
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
    fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.to_frame_with(jd, &AstroContext::default())
    }

    fn to_frame_with<F2: ReferenceFrame>(
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

    fn to<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>,
    {
        self.to_with(jd, &AstroContext::default())
    }

    fn to_with<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
        ctx: &AstroContext,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>,
    {
        // Order: center first (in source frame), then rotate
        <Self as TransformCenter<C2, F, U>>::to_center_with(self, (), *jd, ctx)
            .to_frame_with::<F2>(jd, ctx)
    }
}

// =============================================================================
// WithEngine - Builder for custom context
// =============================================================================

/// A wrapper that pairs a coordinate reference with a custom [`AstroContext`],
/// enabling `.using(&engine).to_frame::<F2>(&jd)` style calls.
///
/// # Example
///
/// ```rust,ignore
/// let engine = AstroContext::new();
/// let result = direction.using(&engine).to_frame::<EclipticMeanJ2000>(&jd);
/// ```
pub struct WithEngine<'a, T> {
    inner: &'a T,
    ctx: &'a AstroContext,
}

/// Helper trait to create [`WithEngine`] wrappers.
pub trait UsingEngine: Sized {
    /// Wrap this coordinate with a custom [`AstroContext`] for the next
    /// transformation call.
    fn using<'a>(&'a self, engine: &'a AstroContext) -> WithEngine<'a, Self> {
        WithEngine {
            inner: self,
            ctx: engine,
        }
    }
}

// Blanket impl: every type gets `.using()`
impl<T> UsingEngine for T {}

// --- WithEngine impls for Direction<F> ---

impl<'a, F: ReferenceFrame> WithEngine<'a, Direction<F>> {
    /// Rotates this direction to a new reference frame using the wrapped context.
    pub fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Direction<F2>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.inner.to_frame_with(jd, self.ctx)
    }
}

// --- WithEngine impls for Position<C, F, U> ---

impl<'a, C, F, U> WithEngine<'a, Position<C, F, U>>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    U: LengthUnit,
{
    /// Rotates this position to a new reference frame using the wrapped context.
    pub fn to_frame<F2: ReferenceFrame>(&self, jd: &JulianDate) -> Position<C, F2, U>
    where
        (): FrameRotationProvider<F, F2>,
    {
        self.inner.to_frame_with(jd, self.ctx)
    }

    /// Translates this position to a new reference center using the wrapped context.
    pub fn to_center<C2: ReferenceCenter<Params = ()>>(&self, jd: &JulianDate) -> Position<C2, F, U>
    where
        (): CenterShiftProvider<C, C2, F>,
    {
        <Position<C, F, U> as TransformCenter<C2, F, U>>::to_center_with(
            self.inner,
            (),
            *jd,
            self.ctx,
        )
    }

    /// Combined center + frame transform using the wrapped context.
    pub fn to<C2: ReferenceCenter<Params = ()>, F2: ReferenceFrame>(
        &self,
        jd: &JulianDate,
    ) -> Position<C2, F2, U>
    where
        (): CenterShiftProvider<C, C2, F>,
        (): FrameRotationProvider<F, F2>,
    {
        self.inner.to_with(jd, self.ctx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::{Barycentric, Geocentric};
    use crate::coordinates::frames::{EclipticMeanJ2000, ICRS};
    use qtty::AstronomicalUnit;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_direction_frame_transform() {
        let dir = Direction::<ICRS>::new(1.0, 0.0, 0.0);
        let jd = JulianDate::J2000;

        // ICRS to EclipticMeanJ2000 includes a small frame-bias; don't assume exact axis invariance.
        let dir_ecl: Direction<EclipticMeanJ2000> = dir.to_frame(&jd);

        // Must be finite and length-preserving.
        assert!(dir_ecl.x().is_finite() && dir_ecl.y().is_finite() && dir_ecl.z().is_finite());
        let n0 = (dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z()).sqrt();
        let n1 =
            (dir_ecl.x() * dir_ecl.x() + dir_ecl.y() * dir_ecl.y() + dir_ecl.z() * dir_ecl.z())
                .sqrt();
        assert!((n0 - n1).abs() < 1e-12);
    }

    #[test]
    fn test_direction_frame_transform_with_ctx() {
        let dir = Direction::<ICRS>::new(1.0, 0.0, 0.0);
        let ctx = AstroContext::default();
        let jd = JulianDate::J2000;

        let dir_ecl: Direction<EclipticMeanJ2000> = dir.to_frame_with(&jd, &ctx);
        let dir_ecl_default: Direction<EclipticMeanJ2000> = dir.to_frame(&jd);

        assert!((dir_ecl.x() - dir_ecl_default.x()).abs() < 1e-15);
        assert!((dir_ecl.y() - dir_ecl_default.y()).abs() < 1e-15);
        assert!((dir_ecl.z() - dir_ecl_default.z()).abs() < 1e-15);
    }

    #[test]
    fn test_direction_frame_roundtrip() {
        let dir = Direction::<ICRS>::new(1.0, 2.0, 3.0);
        let jd = JulianDate::J2000;

        let dir_ecl: Direction<EclipticMeanJ2000> = dir.to_frame(&jd);
        let dir_back: Direction<ICRS> = dir_ecl.to_frame(&jd);

        assert!((dir_back.x() - dir.x()).abs() < EPSILON);
        assert!((dir_back.y() - dir.y()).abs() < EPSILON);
        assert!((dir_back.z() - dir.z()).abs() < EPSILON);
    }

    #[test]
    fn test_position_frame_transform() {
        let pos = Position::<Barycentric, ICRS, AstronomicalUnit>::new(1.0, 0.0, 0.0);
        let jd = JulianDate::J2000;

        let pos_ecl: Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> = pos.to_frame(&jd);

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
        let jd = JulianDate::J2000;

        // A point at the Geocentric origin should map to Earth's position in Barycentric
        let geo_origin =
            Position::<Geocentric, EclipticMeanJ2000, AstronomicalUnit>::new(0.0, 0.0, 0.0);
        let bary: Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> =
            geo_origin.to_center(jd);

        // Should be non-zero (Earth is ~1 AU from barycenter)
        let dist = bary.distance();
        assert!(
            dist > 0.9 && dist < 1.1,
            "Earth should be ~1 AU from barycenter, got {}",
            dist
        );
    }

    #[test]
    fn test_position_combined_transform() {
        let jd = JulianDate::J2000;

        let pos = Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(1.0, 0.5, 0.2);

        // Combined transform: Barycentric EclipticMeanJ2000 -> Geocentric ICRS
        let result: Position<Geocentric, ICRS, AstronomicalUnit> = pos.to(&jd);

        // Verify it's not the same as the input (transformation happened)
        assert!(
            (result.x() - pos.x()).abs() > EPSILON
                || (result.y() - pos.y()).abs() > EPSILON
                || (result.z() - pos.z()).abs() > EPSILON
        );
    }

    #[test]
    fn test_position_identity_transforms() {
        let jd = JulianDate::J2000;

        let pos = Position::<Barycentric, ICRS, AstronomicalUnit>::new(1.5, 2.5, 3.5);

        // Identity frame transform
        let same_frame: Position<Barycentric, ICRS, AstronomicalUnit> = pos.to_frame(&jd);
        assert!((same_frame.x() - pos.x()).abs() < EPSILON);
        assert!((same_frame.y() - pos.y()).abs() < EPSILON);
        assert!((same_frame.z() - pos.z()).abs() < EPSILON);

        // Identity center transform (via ShiftCenter)
        let same_center: Position<Barycentric, ICRS, AstronomicalUnit> = pos.to_center(jd);
        assert!((same_center.x() - pos.x()).abs() < EPSILON);
        assert!((same_center.y() - pos.y()).abs() < EPSILON);
        assert!((same_center.z() - pos.z()).abs() < EPSILON);
    }

    #[test]
    fn test_using_engine() {
        let dir = Direction::<ICRS>::new(1.0, 0.0, 0.0);
        let engine = AstroContext::default();
        let jd = JulianDate::J2000;

        let dir_ecl: Direction<EclipticMeanJ2000> = dir.using(&engine).to_frame(&jd);
        let dir_ecl_direct: Direction<EclipticMeanJ2000> = dir.to_frame(&jd);

        assert!((dir_ecl.x() - dir_ecl_direct.x()).abs() < 1e-15);
        assert!((dir_ecl.y() - dir_ecl_direct.y()).abs() < 1e-15);
        assert!((dir_ecl.z() - dir_ecl_direct.z()).abs() < 1e-15);
    }
}
