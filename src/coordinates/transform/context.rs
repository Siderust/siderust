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
use crate::astro::nutation::NutationModel;
use crate::time::JulianDate;

#[cfg(not(any(feature = "de440", feature = "de441")))]
use crate::calculus::ephemeris::Vsop87Ephemeris;

/// Default ephemeris type.
///
/// - Without a DE feature: [`Vsop87Ephemeris`] (VSOP87 + ELP2000-82B).
/// - With `de441` feature (and real data): `De441Ephemeris` (JPL DE441, compile-time).
/// - With `de440` feature (and real data): `De440Ephemeris` (JPL DE440, compile-time).
/// - With a DE feature but matching `SIDERUST_JPL_STUB` set: falls back to
///   [`Vsop87Ephemeris`] so tests run without downloading the BSP.
/// - For other large datasets: use
///   [`RuntimeEphemeris`](crate::calculus::ephemeris::RuntimeEphemeris) with a
///   BSP file loaded at runtime via [`DataManager`](crate::data::DataManager).
///
/// This type alias is used as the default `Eph` parameter in [`AstroContext`],
/// so all code using `AstroContext::default()` will automatically use the
/// selected backend.
#[cfg(not(any(feature = "de440", feature = "de441")))]
pub type DefaultEphemeris = Vsop87Ephemeris;

#[cfg(all(feature = "de441", not(siderust_mock_de441)))]
pub type DefaultEphemeris = crate::calculus::ephemeris::De441Ephemeris;

#[cfg(all(
    feature = "de440",
    not(feature = "de441"),
    not(siderust_mock_de440)
))]
pub type DefaultEphemeris = crate::calculus::ephemeris::De440Ephemeris;

// Stub: DE feature is on but SIDERUST_JPL_STUB is set, fall back to VSOP87 so
// tests work without the BSP download. DE441 takes precedence when both DE
// features are enabled.
#[cfg(any(
    all(feature = "de441", siderust_mock_de441),
    all(
        feature = "de440",
        not(feature = "de441"),
        siderust_mock_de440
    )
))]
pub type DefaultEphemeris = crate::calculus::ephemeris::Vsop87Ephemeris;

/// Default Earth orientation model: [`IersEop`], backed by the
/// build-time embedded `finals2000A.all` table.
///
/// For zero-overhead use (no EOP corrections), substitute
/// [`NullEop`](crate::astro::eop::NullEop) as the `Eop` type parameter.
pub type DefaultEop = IersEop;

/// Default nutation/precession model marker.
///
/// Uses the full IAU 2006/2000A nutation model for SOFA-grade accuracy.
/// For the abridged 77-term model, substitute
/// [`Iau2000B`](crate::astro::nutation::Iau2000B).
pub type DefaultNutationModel = crate::astro::nutation::Iau2006A;

/// Runtime-selectable Earth orientation model presets.
///
/// These presets control nutation behavior while keeping the current IAU 2006
/// precession implementation used across the transform pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EarthOrientationModel {
    /// Full IAU 2000A nutation.
    Iau2000A,
    /// Abridged IAU 2000B nutation.
    Iau2000B,
    /// IAU 2006 precession-only profile (zero nutation offsets).
    Iau2006,
    /// High-precision IAU 2006A convention (IAU 2006 + corrected 2000A).
    Iau2006A,
}

/// Astronomical context for coordinate transformations.
///
/// This structure holds the runtime configuration for time-dependent
/// transformations:
/// - Ephemeris source for planetary positions
/// - Earth Orientation Parameters (EOP) for polar motion and UT1-UTC
///
/// Compile-time model selection is layered on top through
/// [`ModelContext`], keeping the base context focused on providers/data.
///
/// # Type Parameters
///
/// - `Eph`: Ephemeris provider type (default: [`DefaultEphemeris`]).
/// - `Eop`: Earth Orientation Parameters type (default: [`DefaultEop`]).
///
/// To override the default nutation model, call [`AstroContext::with_model`]
/// and pass the resulting [`ModelContext`] into the `_with` transform APIs.
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
pub struct AstroContext<Eph = DefaultEphemeris, Eop = DefaultEop> {
    _ephemeris: PhantomData<Eph>,
    /// Earth Orientation Parameters provider.
    eop: Eop,
}

impl<Eph, Eop: Default> Default for AstroContext<Eph, Eop> {
    fn default() -> Self {
        Self {
            _ephemeris: PhantomData,
            eop: Eop::default(),
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

impl<Eph, Eop: Default> AstroContext<Eph, Eop> {
    /// Creates a context with custom type parameters.
    ///
    /// Use this when you need to specify custom ephemeris or EOP types.
    #[inline]
    pub fn with_types() -> Self {
        Self::default()
    }
}

impl<Eph, Eop> AstroContext<Eph, Eop> {
    /// Binds a compile-time nutation model to this context.
    ///
    /// This returns a lightweight wrapper carrying only a reference to the
    /// base context plus a [`PhantomData`] marker for the selected model.
    #[inline]
    pub fn with_model<Nut: NutationModel>(&self) -> ModelContext<'_, Eph, Eop, Nut> {
        ModelContext {
            ctx: self,
            _nutation: PhantomData,
        }
    }
}

impl<Eph, Eop: EopProvider> AstroContext<Eph, Eop> {
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
impl<Eop> AstroContext<DefaultEphemeris, Eop> {
    /// Returns true if this context uses the default ephemeris.
    #[inline]
    pub const fn uses_default_ephemeris(&self) -> bool {
        true
    }
}

/// Context wrapper that binds a compile-time nutation model to an
/// [`AstroContext`] without duplicating the runtime provider state.
#[derive(Debug, Clone, Copy)]
pub struct ModelContext<'a, Eph = DefaultEphemeris, Eop = DefaultEop, Nut = DefaultNutationModel> {
    ctx: &'a AstroContext<Eph, Eop>,
    _nutation: PhantomData<Nut>,
}

impl<'a, Eph, Eop, Nut> ModelContext<'a, Eph, Eop, Nut> {
    /// Returns the underlying runtime context.
    #[inline]
    pub fn astro_context(&self) -> &'a AstroContext<Eph, Eop> {
        self.ctx
    }
}

/// Common interface for transformation contexts accepted by the ergonomic
/// `_with` APIs.
///
/// A plain [`AstroContext`] uses [`DefaultNutationModel`], while a
/// [`ModelContext`] overrides that compile-time model through its phantom
/// type parameter.
pub trait TransformContext {
    /// Ephemeris provider type carried by the context.
    type Eph;
    /// EOP provider type carried by the context.
    type Eop: EopProvider;
    /// Compile-time nutation model associated with this context value.
    type Nut: NutationModel;

    /// Returns the underlying runtime context.
    fn astro_context(&self) -> &AstroContext<Self::Eph, Self::Eop>;
}

impl<Eph, Eop: EopProvider> TransformContext for AstroContext<Eph, Eop> {
    type Eph = Eph;
    type Eop = Eop;
    type Nut = DefaultNutationModel;

    #[inline]
    fn astro_context(&self) -> &AstroContext<Self::Eph, Self::Eop> {
        self
    }
}

impl<'a, Eph, Eop: EopProvider, Nut: NutationModel> TransformContext
    for ModelContext<'a, Eph, Eop, Nut>
{
    type Eph = Eph;
    type Eop = Eop;
    type Nut = Nut;

    #[inline]
    fn astro_context(&self) -> &AstroContext<Self::Eph, Self::Eop> {
        self.ctx
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// DynAstroContext, runtime-selected ephemeris via DynEphemeris trait object
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
        let ctx: AstroContext<DefaultEphemeris, DefaultEop> = AstroContext::with_types();
        assert!(ctx.uses_default_ephemeris());
    }

    #[test]
    fn test_context_with_model() {
        let ctx = AstroContext::new();
        let model_ctx = ctx.with_model::<DefaultNutationModel>();
        let _ = model_ctx.astro_context();
    }

    #[test]
    fn test_eop_at_returns_values() {
        let ctx = AstroContext::new();
        let jd = JulianDate::J2000;
        let eop = ctx.eop_at(jd);
        // EOP values should be finite
        assert!(eop.dut1.is_finite());
    }

    #[test]
    fn test_eop_ref() {
        let ctx = AstroContext::new();
        let _eop = ctx.eop();
    }

    #[test]
    fn test_dyn_context_creation() {
        use crate::calculus::ephemeris::Vsop87Ephemeris;

        // Vsop87Ephemeris implements DynEphemeris via blanket impl
        let dyn_ctx = DynAstroContext::with_ephemeris(Box::new(Vsop87Ephemeris));
        let _eph = dyn_ctx.ephemeris();
        let _eop = dyn_ctx.eop_at(JulianDate::J2000);
        let _eop_ref = dyn_ctx.eop();

        // Debug should work
        let debug_str = format!("{:?}", dyn_ctx);
        assert!(debug_str.contains("DynAstroContext"));
    }

    #[test]
    fn test_dyn_context_with_eop() {
        use crate::calculus::ephemeris::Vsop87Ephemeris;

        let dyn_ctx = DynAstroContext::with_ephemeris_and_eop(
            Box::new(Vsop87Ephemeris),
            DefaultEop::default(),
        );
        let eop = dyn_ctx.eop_at(JulianDate::J2000);
        assert!(eop.dut1.is_finite());
    }
}
