//! # Cartesian Coordinates Module
//!
//! This module provides strongly-typed Cartesian coordinate types for astronomical
//! calculations, with mathematical correctness enforced through the type system.
//!
//! ## Mathematical Model
//!
//! The module implements a rigorous distinction between three fundamentally
//! different mathematical objects:
//!
//! ### [`Direction<F>`] — Unit Vector (Orientation)
//!
//! A dimensionless unit vector representing pure orientation in space.
//!
//! - **Frame-dependent**: Orientation is relative to frame `F`
//! - **Center-independent**: Directions are translation-invariant
//! - **Dimensionless**: Magnitude is always 1
//! - **Valid operations**: Rotation (frame transform), dot/cross products
//!
//! ### [`Displacement<F, U>`] — Free Displacement Vector
//!
//! A free vector representing directed magnitude (displacement) in space.
//!
//! - **Frame-dependent**: Direction is relative to frame `F`
//! - **Center-independent**: Displacements are translation-invariant
//! - **Dimensioned**: Has length unit `U`
//! - **Valid operations**: Addition, subtraction, scaling, normalization
//!
//! ### [`Position<C, F, U>`] — Affine Point
//!
//! A point in affine space, representing a location relative to an origin.
//!
//! - **Center-dependent**: Position is measured from center `C`
//! - **Frame-dependent**: Coordinates are relative to frame `F`
//! - **Dimensioned**: Has length unit `U`
//! - **Valid operations**: Subtraction (yields Displacement), translation by Displacement
//!
//! ## Algebraic Rules
//!
//! The type system enforces these mathematical constraints:
//!
//! | Operation | Result | Meaning |
//! |-----------|--------|---------|
//! | `Position - Position` | `Displacement` | Displacement between points |
//! | `Position + Displacement` | `Position` | Translate point |
//! | `Displacement + Displacement` | `Displacement` | Add displacements |
//! | `Direction * Length` | `Displacement` | Scale direction |
//! | `normalize(Displacement)` | `Direction` | Extract orientation |
//!
//! ## Forbidden Operations (compile errors)
//!
//! - `Position + Position` — Adding points has no meaning
//! - `Direction + anything` — Unit vectors aren't additive
//! - Center transform on `Direction` — Directions have no center
//!
//! ## Line of Sight
//!
//! To compute the direction from an observer to a target, use [`line_of_sight`]:
//!
//! ```rust
//! use siderust::coordinates::algebra::cartesian::{line_of_sight, Position};
//! use siderust::coordinates::algebra::centers::Heliocentric;
//! use siderust::coordinates::algebra::frames::Ecliptic;
//! use qtty::*;
//!
//! let observer = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(0.0, 0.0, 0.0);
//! let target = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(1.0, 1.0, 0.0);
//!
//! let direction = line_of_sight(&observer, &target);
//! ```
//!
//! ## Architecture
//!
//! All types share a common internal storage [`XYZ<T>`](xyz::XYZ) that implements
//! component-wise operations once. The semantic types are thin wrappers with
//! `PhantomData` markers for type safety. This provides:
//!
//! - **Zero-cost abstractions**: `#[repr(transparent)]` where applicable
//! - **No code duplication**: Math is centralized in `XYZ<T>`
//! - **Type safety**: Invalid operations are compile errors

// Internal shared storage
mod xyz;

// Semantic types
mod position;      // Position<C, F, U> - Affine point
mod displacement;  // Displacement<F, U> - Free displacement vector
mod direction;     // Direction<F> - Unit vector
mod velocity;      // Velocity<F, U> - Rate of change

// =============================================================================
// Public Re-exports
// =============================================================================

pub use position::Position;
pub use displacement::Displacement;
pub use direction::Direction;
pub use velocity::Velocity;

/// Legacy type alias for backward compatibility.
///
/// In earlier versions, `Vector<C, F, U>` was the main Cartesian type.
/// New code should use:
/// - [`Position<C, F, U>`] for affine points (locations with a center)
/// - [`Displacement<F, U>`] for free vectors (center-independent)
#[deprecated(
    since = "0.4.0",
    note = "Use Position<C, F, U> for positions or Displacement<F, U> for displacements"
)]
pub type Vector<C, F, U> = Position<C, F, U>;

// =============================================================================
// Line of Sight Functions
// =============================================================================

use crate::coordinates::algebra::centers::ReferenceCenter;
use crate::coordinates::algebra::frames::ReferenceFrame;
use qtty::{LengthUnit, Quantity};

/// Computes the line-of-sight direction from an observer to a target.
///
/// This is the mathematically correct way to obtain an observer-dependent
/// direction. The result is a unit direction (free vector) pointing from
/// the observer position toward the target position.
///
/// # Mathematical Definition
///
/// ```text
/// displacement = target - observer
/// direction = normalize(displacement)
/// ```
///
/// # Requirements
///
/// Both positions must be in the **same center and frame**. If they are not,
/// convert them first using center and/or frame transformations.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::algebra::cartesian::{line_of_sight, Position, Direction};
/// use siderust::coordinates::algebra::centers::Heliocentric;
/// use siderust::coordinates::algebra::frames::Ecliptic;
/// use qtty::*;
///
/// let observer = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(0.0, 0.0, 0.0);
/// let target = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(1.0, 1.0, 0.0);
///
/// let los: Direction<Ecliptic> = line_of_sight(&observer, &target);
/// ```
///
/// # Panics
///
/// Panics if the observer and target positions are identical (zero displacement).
#[inline]
pub fn line_of_sight<C, F, U>(observer: &Position<C, F, U>, target: &Position<C, F, U>) -> Direction<F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    // Position - Position -> Displacement
    let displacement: Displacement<F, U> = target.sub(observer);
    
    // normalize(Displacement) -> Direction
    displacement
        .normalize()
        .expect("line_of_sight requires distinct observer and target positions")
}

/// Computes the line-of-sight direction and distance from an observer to a target.
///
/// Returns both the unit direction and the distance (magnitude of the displacement).
/// This is useful when you need both the pointing direction and the range.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::algebra::cartesian::{line_of_sight_with_distance, Position};
/// use siderust::coordinates::algebra::centers::Heliocentric;
/// use siderust::coordinates::algebra::frames::Ecliptic;
/// use qtty::*;
///
/// let observer = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(0.0, 0.0, 0.0);
/// let target = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(3.0, 4.0, 0.0);
///
/// let (dir, dist) = line_of_sight_with_distance(&observer, &target);
/// assert!((dist.value() - 5.0).abs() < 1e-10);
/// assert!((dir.x() - 0.6).abs() < 1e-10);
/// assert!((dir.y() - 0.8).abs() < 1e-10);
/// ```
///
/// # Panics
///
/// Panics if the observer and target positions are identical.
#[inline]
pub fn line_of_sight_with_distance<C, F, U>(
    observer: &Position<C, F, U>,
    target: &Position<C, F, U>,
) -> (Direction<F>, Quantity<U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    let displacement: Displacement<F, U> = target.sub(observer);
    let distance = displacement.magnitude();
    let direction = displacement
        .normalize()
        .expect("line_of_sight requires distinct observer and target positions");
    
    (direction, distance)
}

/// Attempts to compute line-of-sight, returning `None` if positions are identical.
///
/// This is the non-panicking version of [`line_of_sight`].
#[inline]
pub fn try_line_of_sight<C, F, U>(
    observer: &Position<C, F, U>,
    target: &Position<C, F, U>,
) -> Option<Direction<F>>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    let displacement: Displacement<F, U> = target.sub(observer);
    displacement.normalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::algebra::centers::Heliocentric;
    use crate::coordinates::algebra::frames::Ecliptic;
    use qtty::AstronomicalUnit;

    type PosAu = Position<Heliocentric, Ecliptic, AstronomicalUnit>;

    #[test]
    fn test_line_of_sight_basic() {
        let observer = PosAu::new(0.0, 0.0, 0.0);
        let target = PosAu::new(3.0, 4.0, 0.0);

        let los = line_of_sight(&observer, &target);
        assert!((los.x() - 0.6).abs() < 1e-12);
        assert!((los.y() - 0.8).abs() < 1e-12);
        assert!(los.z().abs() < 1e-12);
    }

    #[test]
    fn test_line_of_sight_with_distance() {
        let observer = PosAu::new(1.0, 1.0, 1.0);
        let target = PosAu::new(4.0, 5.0, 1.0);

        let (dir, dist) = line_of_sight_with_distance(&observer, &target);
        assert!((dist.value() - 5.0).abs() < 1e-12);
        assert!((dir.x() - 0.6).abs() < 1e-12);
        assert!((dir.y() - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_try_line_of_sight_same_position() {
        let pos = PosAu::new(1.0, 2.0, 3.0);
        assert!(try_line_of_sight(&pos, &pos).is_none());
    }

    #[test]
    fn test_position_displacement_roundtrip() {
        let a = PosAu::new(1.0, 2.0, 3.0);
        let b = PosAu::new(5.0, 7.0, 11.0);

        // a + (b - a) == b
        let displacement = b - a;
        let result = a + displacement;

        assert!((result.x().value() - b.x().value()).abs() < 1e-12);
        assert!((result.y().value() - b.y().value()).abs() < 1e-12);
        assert!((result.z().value() - b.z().value()).abs() < 1e-12);
    }

    #[test]
    fn test_direction_scale_to_displacement() {
        use qtty::AstronomicalUnit;
        
        let dir = Direction::<Ecliptic>::new(1.0, 0.0, 0.0);
        let disp: Displacement<Ecliptic, AstronomicalUnit> = dir * Quantity::new(3.0);

        assert!((disp.x().value() - 3.0).abs() < 1e-12);
        assert!(disp.y().value().abs() < 1e-12);
    }

    #[test]
    fn test_displacement_normalize_to_direction() {
        use qtty::AstronomicalUnit;

        let disp = Displacement::<Ecliptic, AstronomicalUnit>::new(3.0, 4.0, 0.0);
        let dir = disp.normalize().expect("non-zero displacement");

        // Check normalization: 3/5 = 0.6, 4/5 = 0.8
        assert!((dir.x() - 0.6).abs() < 1e-12);
        assert!((dir.y() - 0.8).abs() < 1e-12);
        assert!(dir.z().abs() < 1e-12);
    }
}
