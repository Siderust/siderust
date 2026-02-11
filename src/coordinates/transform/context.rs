// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical Context
//!
//! This module provides the [`AstroContext`] type, which holds configuration
//! for astronomical transformations, including ephemeris selection, Earth
//! orientation parameters, and model choices.
//!
//! ## Design Philosophy
//!
//! The context is designed to be:
//! - **Generic**: Parameterized over ephemeris and model types for flexibility.
//! - **Zero-cost when using defaults**: With `Vsop87Ephemeris` (a ZST), there's no overhead.
//! - **Pluggable**: Custom ephemeris backends (e.g., JPL DE) can be used by changing the type.
//! - **Passed by reference**: Context is borrowed, not consumed, in transforms.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::transform::context::AstroContext;
//!
//! // Create a default context (uses VSOP87)
//! let ctx = AstroContext::new();
//!
//! // Use with transforms:
//! // position.to_frame::<Ecliptic>(&jd, &ctx);
//! ```
//!
//! ## Custom Ephemeris
//!
//! To use a custom ephemeris backend:
//!
//! ```rust,ignore
//! use siderust::coordinates::transform::context::AstroContext;
//! use my_crate::JplDeEphemeris;
//!
//! let jpl = JplDeEphemeris::load("de440.bsp")?;
//! let ctx = AstroContext::with_ephemeris(jpl);
//! ```

use std::marker::PhantomData;

use super::ephemeris::Vsop87Ephemeris;

/// Default ephemeris: uses built-in VSOP87 planetary theory.
///
/// This is a type alias for [`Vsop87Ephemeris`], which is a zero-sized type (ZST).
/// Using the default ephemeris has no runtime overhead.
pub type DefaultEphemeris = Vsop87Ephemeris;

/// Default Earth orientation model marker.
#[derive(Debug, Clone, Copy, Default)]
pub struct DefaultEop;

/// Default nutation/precession model marker.
#[derive(Debug, Clone, Copy, Default)]
pub struct DefaultNutationModel;

/// Astronomical context for coordinate transformations.
///
/// This structure holds configuration for time-dependent transformations,
/// including:
/// - Ephemeris source for planetary positions
/// - Earth Orientation Parameters (EOP) for polar motion and UT1-UTC
/// - Nutation/precession model selection
///
/// # Type Parameters
///
/// - `Eph`: Ephemeris provider type (default: [`DefaultEphemeris`] = [`Vsop87Ephemeris`]).
/// - `Eop`: Earth Orientation Parameters type (default: [`DefaultEop`]).
/// - `Nut`: Nutation/precession model type (default: [`DefaultNutationModel`]).
///
/// # Zero-Cost Abstraction
///
/// When using `DefaultEphemeris` (which is `Vsop87Ephemeris`, a zero-sized type),
/// `AstroContext` has no runtime overhead. The compiler fully inlines ephemeris
/// calls, making this equivalent to direct VSOP87 function calls.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::transform::context::AstroContext;
///
/// // Use defaults for typical applications
/// let ctx = AstroContext::new();
/// ```
#[derive(Debug, Clone)]
pub struct AstroContext<Eph = DefaultEphemeris, Eop = DefaultEop, Nut = DefaultNutationModel> {
    /// The ephemeris provider for planetary positions/velocities.
    eph: Eph,
    /// Earth Orientation Parameters (phantom for now).
    _eop: PhantomData<Eop>,
    /// Nutation/precession model (phantom for now).
    _nutation: PhantomData<Nut>,
}

impl Default for AstroContext<DefaultEphemeris, DefaultEop, DefaultNutationModel> {
    fn default() -> Self {
        Self {
            eph: Vsop87Ephemeris,
            _eop: PhantomData,
            _nutation: PhantomData,
        }
    }
}

impl AstroContext {
    /// Creates a new context with default configuration (VSOP87 ephemeris).
    ///
    /// This is equivalent to `AstroContext::default()` but provides a named
    /// constructor for clarity.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }
}

impl<Eph, Eop, Nut> AstroContext<Eph, Eop, Nut> {
    /// Returns a reference to the ephemeris provider.
    ///
    /// Use this to query planetary positions and velocities:
    ///
    /// ```rust
    /// use siderust::coordinates::transform::context::AstroContext;
    /// use siderust::coordinates::transform::ephemeris::{BodyId, BodyEphemeris};
    /// use siderust::astro::JulianDate;
    ///
    /// let ctx = AstroContext::new();
    /// let earth_pos = ctx.ephemeris().position_barycentric(BodyId::Earth, JulianDate::J2000);
    /// ```
    #[inline]
    pub fn ephemeris(&self) -> &Eph {
        &self.eph
    }

    /// Creates a context with a custom ephemeris provider.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let jpl_eph = JplDeEphemeris::load("de440.bsp")?;
    /// let ctx = AstroContext::with_ephemeris(jpl_eph);
    /// ```
    #[inline]
    pub fn with_ephemeris<E>(eph: E) -> AstroContext<E, Eop, Nut>
    where
        Eop: Default,
        Nut: Default,
    {
        AstroContext {
            eph,
            _eop: PhantomData,
            _nutation: PhantomData,
        }
    }
}

// Marker that this context uses default ephemeris
impl<Eop, Nut> AstroContext<DefaultEphemeris, Eop, Nut> {
    /// Returns true if this context uses the default ephemeris (VSOP87).
    #[inline]
    pub const fn uses_default_ephemeris(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::JulianDate;
    use crate::coordinates::transform::ephemeris::{BodyEphemeris, BodyId};

    #[test]
    fn test_context_creation() {
        let ctx = AstroContext::new();
        assert!(ctx.uses_default_ephemeris());
    }

    #[test]
    fn test_context_default() {
        let ctx: AstroContext = Default::default();
        assert!(ctx.uses_default_ephemeris());
    }

    #[test]
    fn test_ephemeris_accessor() {
        let ctx = AstroContext::new();
        let pos = ctx.ephemeris().position_barycentric(BodyId::Earth, JulianDate::J2000);

        // Earth should be roughly 1 AU from barycenter
        let dist = (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt();
        assert!((dist - 1.0).abs() < 0.02);
    }

    #[test]
    fn test_zero_size() {
        // AstroContext with ZST ephemeris should be zero-sized
        assert_eq!(
            std::mem::size_of::<AstroContext>(),
            0,
            "AstroContext with Vsop87Ephemeris should be zero-sized"
        );
    }
}
