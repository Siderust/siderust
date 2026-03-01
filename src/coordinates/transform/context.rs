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
//! // position.to_frame::<EclipticMeanJ2000>(&jd, &ctx);
//! ```

use std::marker::PhantomData;

use crate::astro::eop::{EopProvider, EopValues, IersEop};
use crate::time::JulianDate;

#[cfg(not(feature = "de440"))]
use crate::calculus::ephemeris::Vsop87Ephemeris;

/// Default ephemeris type.
///
/// - Without `de440` feature: [`Vsop87Ephemeris`] (VSOP87 + ELP2000-82B).
/// - With `de440` feature: `De440Ephemeris` (JPL DE440, compile-time).
/// - For DE441 or other large datasets: use [`RuntimeEphemeris`](crate::calculus::ephemeris::RuntimeEphemeris)
///   with a BSP file loaded at runtime via [`DataManager`](crate::data::DataManager).
///
/// This type alias is used as the default `Eph` parameter in [`AstroContext`],
/// so all code using `AstroContext::default()` will automatically use the
/// selected backend.
#[cfg(not(feature = "de440"))]
pub type DefaultEphemeris = Vsop87Ephemeris;

#[cfg(feature = "de440")]
pub type DefaultEphemeris = crate::calculus::ephemeris::De440Ephemeris;

/// Default Earth orientation model: [`IersEop`], backed by the
/// build-time embedded `finals2000A.all` table.
///
/// For zero-overhead use (no EOP corrections), substitute
/// [`NullEop`](crate::astro::eop::NullEop) as the `Eop` type parameter.
pub type DefaultEop = IersEop;

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
    /// Earth Orientation Parameters provider.
    eop: Eop,
    _nutation: PhantomData<Nut>,
}

impl<Eph, Eop: Default, Nut> Default for AstroContext<Eph, Eop, Nut> {
    fn default() -> Self {
        Self {
            _ephemeris: PhantomData,
            eop: Eop::default(),
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

impl<Eph, Eop: Default, Nut> AstroContext<Eph, Eop, Nut> {
    /// Creates a context with custom type parameters.
    ///
    /// Use this when you need to specify custom ephemeris or model types.
    #[inline]
    pub fn with_types() -> Self {
        Self::default()
    }
}

impl<Eph, Eop: EopProvider, Nut> AstroContext<Eph, Eop, Nut> {
    /// Look up EOP values for the given **UTC** Julian Date.
    ///
    /// # Time-scale contract
    /// `jd_utc` **must** be a UTC Julian Date.  Passing TT or UT1 values will
    /// corrupt the interpolated `dUT1`, `xp`, and `yp` values because the IERS
    /// tables are indexed by UTC civil date (see [`crate::astro::eop`] for
    /// details).
    ///
    /// Delegates to the context's [`EopProvider`].
    #[inline]
    pub fn eop_at(&self, jd_utc: JulianDate) -> EopValues {
        self.eop.eop_at(jd_utc)
    }

    /// Reference to the underlying EOP provider.
    #[inline]
    pub fn eop(&self) -> &Eop {
        &self.eop
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

// ═══════════════════════════════════════════════════════════════════════════
// DynAstroContext — runtime-selected ephemeris via DynEphemeris trait object
// ═══════════════════════════════════════════════════════════════════════════

use crate::calculus::ephemeris::DynEphemeris;

/// Astronomical context with a runtime-selected ephemeris backend.
///
/// This is the dynamic counterpart to [`AstroContext`]. Instead of selecting
/// an ephemeris backend at compile time via type parameters, it stores a
/// `Box<dyn DynEphemeris>` and dispatches via virtual calls.
///
/// Use this when:
/// - You load BSP files at runtime via [`RuntimeEphemeris`](crate::calculus::ephemeris::RuntimeEphemeris)
/// - You need to switch between backends without recompiling
/// - You want to override the default ephemeris at runtime
///
/// # Example
///
/// ```rust,ignore
/// use siderust::calculus::ephemeris::{RuntimeEphemeris, DynEphemeris};
/// use siderust::coordinates::transform::context::DynAstroContext;
///
/// let eph = RuntimeEphemeris::from_bsp("path/to/de441.bsp")?;
/// let ctx = DynAstroContext::with_ephemeris(Box::new(eph));
/// ```
pub struct DynAstroContext<Eop = DefaultEop> {
    ephemeris: Box<dyn DynEphemeris>,
    eop: Eop,
}

impl DynAstroContext<DefaultEop> {
    /// Create a dynamic context from a `DynEphemeris` implementor.
    pub fn with_ephemeris(eph: Box<dyn DynEphemeris>) -> Self {
        Self {
            ephemeris: eph,
            eop: DefaultEop::default(),
        }
    }
}

impl<Eop: Default> DynAstroContext<Eop> {
    /// Create a dynamic context with a custom EOP provider.
    pub fn with_ephemeris_and_eop(eph: Box<dyn DynEphemeris>, eop: Eop) -> Self {
        Self {
            ephemeris: eph,
            eop,
        }
    }
}

impl<Eop: EopProvider> DynAstroContext<Eop> {
    /// Reference to the runtime ephemeris backend.
    #[inline]
    pub fn ephemeris(&self) -> &dyn DynEphemeris {
        &*self.ephemeris
    }

    /// Look up EOP values for the given **UTC** Julian Date.
    #[inline]
    pub fn eop_at(&self, jd_utc: JulianDate) -> EopValues {
        self.eop.eop_at(jd_utc)
    }

    /// Reference to the underlying EOP provider.
    #[inline]
    pub fn eop(&self) -> &Eop {
        &self.eop
    }
}

impl<Eop> std::fmt::Debug for DynAstroContext<Eop>
where
    Eop: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DynAstroContext")
            .field("ephemeris", &"<dyn DynEphemeris>")
            .field("eop", &self.eop)
            .finish()
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
