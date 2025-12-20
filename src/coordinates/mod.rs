//! # Coordinates Module
//!
//! This module defines strongly typed spherical and cartesian coordinate systems used in astronomy.
//! The coordinate systems are implemented using Rust's type system with phantom types
//! to enforce compile-time safety. These phantom types represent the **frame** and **center** of
//! the coordinate system, ensuring that operations between incompatible coordinate systems are
//! disallowed unless explicitly converted. Moreover, thanks to the Unit module we can distinguish
//! the different vector types such as Directions (unitless), Position (Lenght Units) and Velocity
//! (Velocity Units), that enforce the compiler to validate any transformation of coordinates.

//! ## Key Concepts
//! - **Position, Direction and Velocity Types**: Both spherical and cartesian coordinates are parameterized
//!   by a reference center (e.g., `Heliocentric`, `Geocentric`), a reference frame (e.g., `Ecliptic`, `Equatorial`, `ICRS`),
//!   and a measure unit (`Unitless`, `LengthUnit`, `VelocityUnit`). This ensures that only compatible coordinates can be used together.
//! - **Phantom Types**: The `Center`, `Frame` and `Unit`types are zero-cost markers that encode coordinate semantics at compile time.
//! - **Type Safety**: Operations between coordinates are only allowed when their type parameters match, preventing accidental mixing of frames, centers or magnitude.
//! - **Conversions**: Seamless conversion between spherical and cartesian forms, and between different frames and centers, is provided via `From`/`Into` and the `Transform` trait.
//!
//! ## Architectural Separation
//!
//! The coordinate system is now organized into two main modules:
//!
//! ### `algebra` - Pure Mathematical Structures
//!
//! Contains the abstract algebraic coordinate types independent of physical context:
//! - **Reference frames and centers**: Trait definitions for orientation and origin
//! - **Cartesian types**: Vector, Direction, Position, Velocity
//! - **Spherical types**: Direction, Position (base implementations)
//!
//! ### `astro` - Physical/Astronomical Implementations
//!
//! Contains domain-specific coordinate systems with astronomical conventions:
//! - Frame-specific convenience constructors (e.g., `new_ecliptic(lon, lat)`)
//! - Astronomical naming conventions (RA/Dec, lon/lat, alt/az)
//! - Physical coordinate system extensions
//!
//! ### Legacy Compatibility
//!
//! For backward compatibility, the original module structure is maintained with
//! re-exports from the new `algebra` and `astro` modules:
//! - `coordinates::frames` → `algebra::frames`
//! - `coordinates::centers` → `algebra::centers`
//! - `coordinates::cartesian` → `algebra::cartesian`
//! - `coordinates::spherical` → `astro::spherical`
//!
//! ## Coordinate Transform Architecture
//!
//! The coordinate system maintains a clean separation of concerns:
//!
//! - **Center transforms** (translations): Apply only to positions. Moving from geocentric to
//!   heliocentric is a pure vector subtraction. No observation effects.
//!
//! - **Frame transforms** (rotations): Apply to positions, directions, and velocities.
//!   Changing from ecliptic to equatorial is a pure rotation matrix.
//!
//! - **Observation transforms** (in [`observation`] module): Observer-dependent effects like
//!   aberration. These require explicit `ObserverState` and produce directions with explicit
//!   observational state (`Astrometric` or `Apparent`).
//!
//! ## Supported Reference Frames and Centers
//! - **Frames**: `Equatorial`, `Ecliptic`, `Horizontal`, `ICRS`, `ECEF`
//! - **Centers**: `Heliocentric`, `Geocentric`, `Barycentric`, `Topocentric`, `Bodycentric`
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical;
//! use siderust::coordinates::cartesian;
//! use siderust::coordinates::algebra::frames::Ecliptic;
//! use qtty::*;
//!
//! // Create an ecliptic spherical direction (frame-only, no center)
//! let spherical = spherical::Direction::<Ecliptic>::new(
//!     45.0 * DEG, 7.0 * DEG
//! );
//!
//! // Convert to cartesian coordinates
//! let cartesian: cartesian::Direction<Ecliptic> = spherical.to_cartesian();
//!
//! // Convert back to spherical coordinates
//! let spherical_converted: spherical::Direction<Ecliptic> = cartesian.to_spherical();
//!
//! println!("Spherical -> Cartesian -> Spherical: {:?}", spherical_converted);
//! ```
//!
//! ## Submodules
//! - **algebra**: Pure mathematical coordinate structures (frames, centers, vectors)
//! - **astro**: Physical/astronomical coordinate implementations
//! - **transform**: Generic transformations between coordinate systems and frames
//! - **observation**: Observational state types (`Astrometric`, `Apparent`) and aberration
//!
//! ## Legacy Exports (for backward compatibility)
//! - **cartesian**: Re-export of `algebra::cartesian`
//! - **spherical**: Re-export of `astro::spherical`
//! - **frames**: Re-export of `algebra::frames`
//! - **centers**: Re-export of `algebra::centers`

pub mod algebra;
pub mod astro;
pub mod observation;
pub mod transform;

// Legacy re-exports for backward compatibility
pub mod cartesian {
    // Re-export algebraic types
    pub use crate::coordinates::algebra::cartesian::*;
    
    // Re-export astronomical type aliases from astro module
    pub use crate::coordinates::astro::cartesian::{direction, position, velocity};
}

pub mod spherical {
    pub use crate::coordinates::astro::spherical::*;
}

pub mod frames {
    pub use crate::coordinates::algebra::frames::*;
}

pub mod centers {
    pub use crate::coordinates::algebra::centers::*;
}
