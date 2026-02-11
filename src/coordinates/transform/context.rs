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
//! - **Zero-cost when unused**: Default type parameters allow simple usage.
//! - **Passed by reference**: Context is borrowed, not consumed, in transforms.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::transform::context::AstroContext;
//!
//! // Create a default context
//! let ctx = AstroContext::new();
//!
//! // Use with transforms:
//! // position.to_frame::<Ecliptic>(&jd, &ctx);
//! ```

use std::marker::PhantomData;

#[cfg(not(feature = "de440"))]
use crate::calculus::ephemeris::Vsop87Ephemeris;

/// Default ephemeris type.
///
/// - Without `de440` feature: [`Vsop87Ephemeris`] (VSOP87 + ELP2000-82B).
/// - With `de440` feature: `De440Ephemeris` (JPL DE440).
///
/// This type alias is used as the default `Eph` parameter in [`AstroContext`],
/// so all code using `AstroContext::default()` will automatically use the
/// selected backend.
#[cfg(not(feature = "de440"))]
pub type DefaultEphemeris = Vsop87Ephemeris;

#[cfg(feature = "de440")]
pub type DefaultEphemeris = crate::calculus::de440::De440Ephemeris;

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
/// - `Eph`: Ephemeris provider type (default: [`DefaultEphemeris`]).
/// - `Eop`: Earth Orientation Parameters type (default: [`DefaultEop`]).
/// - `Nut`: Nutation/precession model type (default: [`DefaultNutationModel`]).
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
    _ephemeris: PhantomData<Eph>,
    _eop: PhantomData<Eop>,
    _nutation: PhantomData<Nut>,
}

impl<Eph, Eop, Nut> Default for AstroContext<Eph, Eop, Nut> {
    fn default() -> Self {
        Self {
            _ephemeris: PhantomData,
            _eop: PhantomData,
            _nutation: PhantomData,
        }
    }
}

impl AstroContext {
    /// Creates a new context with default configuration.
    ///
    /// This is equivalent to `AstroContext::default()` but provides a named
    /// constructor for clarity.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }
}

impl<Eph, Eop, Nut> AstroContext<Eph, Eop, Nut> {
    /// Creates a context with custom type parameters.
    ///
    /// Use this when you need to specify custom ephemeris or model types.
    #[inline]
    pub fn with_types() -> Self {
        Self::default()
    }
}

// Marker that this context uses default ephemeris
impl<Eop, Nut> AstroContext<DefaultEphemeris, Eop, Nut> {
    /// Returns true if this context uses the default ephemeris.
    #[inline]
    pub const fn uses_default_ephemeris(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_context_with_types() {
        let ctx: AstroContext<DefaultEphemeris, DefaultEop, DefaultNutationModel> =
            AstroContext::with_types();
        assert!(ctx.uses_default_ephemeris());
    }
}
