// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Ephemeris Provider Abstraction
//!
//! This module defines the [`Ephemeris`] trait, which abstracts over different
//! ephemeris backends (VSOP87/ELP2000, DE440/DE441, etc.) for computing fundamental
//! solar-system body positions and velocities.
//!
//! ## Design
//!
//! - **Compile-time selection**: Backends are zero-sized marker types; all
//!   dispatch is monomorphized away.
//! - **Feature-gated**: VSOP87 is always available. DE440 and DE441 require
//!   the `de440`/`de441` Cargo features.
//! - **Body-centric user API preserved**: `Earth::vsop87e(jd)` etc. remain
//!   unchanged. The `Ephemeris` trait is used internally by the coordinate
//!   transform pipeline.
//!
//! ## Available Backends
//!
//! | Backend              | Feature   | Source                  |
//! |----------------------|-----------|-------------------------|
//! | [`Vsop87Ephemeris`]  | (always)  | VSOP87 + ELP2000-82B    |
//! | `De440Ephemeris`     | `de440`   | JPL DE440               |
//! | `De441Ephemeris`     | `de441`   | JPL DE441 (part-2 BSP)  |
//!
//! ## Usage
//!
//! The [`DefaultEphemeris`](crate::coordinates::transform::context::DefaultEphemeris)
//! type alias in [`AstroContext`](crate::coordinates::transform::context::AstroContext)
//! automatically selects the appropriate backend based on enabled features.
//!
//! Users who need explicit control can parameterize `AstroContext` directly:
//! ```rust,ignore
//! use siderust::calculus::ephemeris::Vsop87Ephemeris;
//! use siderust::coordinates::transform::context::AstroContext;
//!
//! let ctx: AstroContext<Vsop87Ephemeris> = AstroContext::with_types();
//! ```

mod vsop87_backend;

#[cfg(feature = "de440")]
mod de440_backend;
#[cfg(feature = "de441")]
mod de441_backend;

pub use vsop87_backend::Vsop87Ephemeris;

#[cfg(feature = "de440")]
pub use de440_backend::De440Ephemeris;
#[cfg(feature = "de441")]
pub use de441_backend::De441Ephemeris;

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::Ecliptic,
};
use crate::targets::Target;
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, Day, Kilometer};

/// Velocity measured in AU/day.
pub type AuPerDay = qtty::Per<AstronomicalUnit, Day>;

/// Trait abstracting over ephemeris backends for fundamental body states.
///
/// All methods are associated functions (no `&self` receiver) because backends
/// are zero-sized marker types — the coefficient data lives in static arrays
/// generated at compile time.
///
/// The five methods cover every ephemeris call site in the coordinate transform
/// pipeline, aberration, and observer-state modules.
pub trait Ephemeris {
    /// Sun's position in barycentric ecliptic coordinates (AU).
    ///
    /// Used by: Heliocentric ↔ Barycentric center shifts.
    fn sun_barycentric(jd: JulianDate)
        -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>;

    /// Earth's position in barycentric ecliptic coordinates (AU).
    ///
    /// Used by: Geocentric ↔ Barycentric center shifts.
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>;

    /// Earth's position in heliocentric ecliptic coordinates (AU).
    ///
    /// Used by: Geocentric ↔ Heliocentric center shifts.
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>>;

    /// Earth's velocity in barycentric ecliptic coordinates (AU/day).
    ///
    /// Used by: annual aberration, observer state velocity.
    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<Ecliptic, AuPerDay>;

    /// Moon's position in geocentric ecliptic coordinates (km).
    ///
    /// Used by: Moon topocentric/horizontal pipeline.
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, Ecliptic, Kilometer>;
}
