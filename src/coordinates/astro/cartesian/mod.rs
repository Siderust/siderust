//! Cartesian coordinate type aliases for common astronomical systems.
//!
//! This module provides convenient type aliases that combine algebraic cartesian
//! types with standard astronomical reference frames and centers.
//!
//! ## Semantic Types
//!
//! - [`Position`]: Affine points with center, frame, and unit
//! - [`Displacement`]: Displacement vectors with frame and unit (center-independent)
//! - [`Direction`]: Unit vectors with frame only (dimensionless)
//! - [`Velocity`]: Rate-of-change vectors with frame and unit
//! - [`Vector`]: Generic free vector (base type for Displacement and Velocity)

use crate::coordinates::algebra::{centers, frames};

/// Re-export the algebraic types
pub use crate::coordinates::algebra::cartesian::{Direction, Displacement, Position, Vector, Velocity};

// =============================================================================
// Direction type aliases (frame-only, no center, dimensionless)
// =============================================================================

pub mod direction {
    pub use super::Direction;
    use super::frames;

    /// **Ecliptic** cartesian direction (unit vector).
    pub type Ecliptic = Direction<frames::Ecliptic>;
    /// **Equatorial** cartesian direction (unit vector).
    pub type Equatorial = Direction<frames::Equatorial>;
    /// **Horizontal** cartesian direction (unit vector).
    pub type Horizontal = Direction<frames::Horizontal>;
    /// **Geographic (ECEF)** cartesian direction (unit vector).
    pub type Geographic = Direction<frames::ECEF>;
    /// **ICRS** cartesian direction (unit vector).
    pub type ICRS = Direction<frames::ICRS>;
}

// =============================================================================
// Displacement type aliases (frame + unit, no center)
// =============================================================================

pub mod displacement {
    pub use super::Displacement;
    use super::frames;

    /// **Ecliptic** displacement vector.
    pub type Ecliptic<U> = Displacement<frames::Ecliptic, U>;
    /// **Equatorial** displacement vector.
    pub type Equatorial<U> = Displacement<frames::Equatorial, U>;
    /// **Horizontal** displacement vector.
    pub type Horizontal<U> = Displacement<frames::Horizontal, U>;
    /// **ICRS** displacement vector.
    pub type ICRS<U> = Displacement<frames::ICRS, U>;
}

// =============================================================================
// Position type aliases (center + frame + unit)
// =============================================================================

pub mod position {
    pub use super::Position;
    use super::{centers, frames};

    /// **Heliocentric Ecliptic** cartesian position.
    pub type Ecliptic<U, C = centers::Heliocentric> = Position<C, frames::Ecliptic, U>;
    /// **Geocentric Equatorial** cartesian position.
    pub type Equatorial<U, C = centers::Geocentric> = Position<C, frames::Equatorial, U>;
    /// **Topocentric Horizontal** cartesian position.
    pub type Horizontal<U, C = centers::Topocentric> = Position<C, frames::Horizontal, U>;
    /// **Geocentric Geographic (ECEF)** cartesian position.
    pub type Geographic<U, C = centers::Geocentric> = Position<C, frames::ECEF, U>;
    /// **Barycentric ICRS** cartesian position.
    pub type ICRS<U, C = centers::Barycentric> = Position<C, frames::ICRS, U>;
    /// **Heliocentric ICRS** cartesian position.
    pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
    /// **Geocentric ICRS** cartesian position.
    pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;
    /// **Topocentric ICRS** cartesian position.
    pub type TCRS<U> = Position<centers::Topocentric, frames::ICRS, U>;
}

// =============================================================================
// Velocity type aliases (frame + unit, no center)
// =============================================================================

pub mod velocity {
    pub use super::Velocity;
    use super::frames;

    /// **Ecliptic** cartesian velocity vector.
    pub type Ecliptic<U> = Velocity<frames::Ecliptic, U>;
    /// **Equatorial** cartesian velocity vector.
    pub type Equatorial<U> = Velocity<frames::Equatorial, U>;
    /// **Horizontal** cartesian velocity vector.
    pub type Horizontal<U> = Velocity<frames::Horizontal, U>;
    /// **ICRS** cartesian velocity vector.
    pub type ICRS<U> = Velocity<frames::ICRS, U>;
}
