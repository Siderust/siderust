// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical Transform Context
//!
//! ## Scientific scope
//!
//! Every time-dependent coordinate transformation in siderust needs access to
//! three external providers: an *ephemeris* that supplies planetary positions,
//! an *Earth Orientation Parameters* (EOP) table that supplies the observed
//! offset between UT1 and UTC and the pole coordinates, and a *nutation /
//! precession model* that determines which IAU series is used to orient the
//! terrestrial frame in the celestial reference system. These three concerns
//! are bundled in this module so that transform call sites carry no ambient
//! state and the same transform path can be reused with any combination of
//! backends.
//!
//! ## Technical scope
//!
//! - [`AstroContext<Eph, Eop>`] — the primary runtime provider container.
//!   - `Eph` selects the ephemeris backend (default: [`DefaultEphemeris`],
//!     VSOP87/ELP2000; DE440/DE441 selectable via Cargo features).
//!   - `Eop` selects the EOP source (default: [`DefaultEop`] = [`IersEop`],
//!     backed by the build-time `finals2000A.all` table; substitute
//!     [`NullEop`](crate::astro::eop::NullEop) for zero overhead).
//! - [`AstroContext::with_model::<Nut>()`] — attaches a compile-time nutation
//!   model marker, returning a zero-cost [`ModelContext`]. The nutation model
//!   is always a **phantom type**: use [`Iau2006A`](crate::astro::nutation::Iau2006A)
//!   (default, full 1365-term MHB2000 + P03 correction),
//!   [`Iau2000B`](crate::astro::nutation::Iau2000B) (77-term abridged),
//!   [`Iau2000A`](crate::astro::nutation::Iau2000A) (uncorrected 2000A), or
//!   [`Iau2006`](crate::astro::nutation::Iau2006) (precession-only). This is
//!   the canonical model-selection mechanism — no runtime enum dispatch.
//! - [`ModelContext<Eph, Eop, Nut>`] — lightweight borrow wrapper that satisfies
//!   [`TransformContext`] and exposes the chosen model to the transform
//!   providers via a single trait.
//! - [`DynAstroContext<Eop>`] — dynamic variant carrying a `Box<dyn DynEphemeris>`
//!   for cases where the BSP file is loaded at runtime.
//!
//! ## Model-selection example
//!
//! ```rust
//! use siderust::coordinates::transform::context::AstroContext;
//! use siderust::astro::nutation::Iau2000B;
//!
//! // Default context: IAU 2006A nutation, VSOP87 ephemeris, IERS EOP.
//! let ctx = AstroContext::new();
//!
//! // Abridged 77-term IAU 2000B nutation — zero runtime overhead.
//! let low_cost = ctx.with_model::<Iau2000B>();
//! // low_cost can now be passed to any _with transform variant.
//! ```
//!
//! ## References
//!
//! - IAU 2006 Resolution B1 (precession) and Resolution B1.6 (nutation).
//! - IERS Conventions (2010), §5.4 and §5.5.
//! - Wallace, P. T., & Capitaine, N. (2006). *A&A* 459, 981 (P03 correction).
//! - SOFA software collection.

use std::marker::PhantomData;

use crate::astro::eop::{EopError, EopProvider, EopValues, IersEop};
use crate::astro::nutation::NutationModel;
use crate::time::JulianDate;

use crate::ephemeris::Vsop87Ephemeris;

/// Default ephemeris type: [`Vsop87Ephemeris`] (VSOP87 + ELP2000-82B).
///
/// For larger JPL datasets, use [`RuntimeEphemeris`](crate::ephemeris::RuntimeEphemeris)
/// with a BSP file loaded at runtime.
pub type DefaultEphemeris = Vsop87Ephemeris;

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
    /// Fallibly look up EOP values for the given **UTC** Julian Date.
    #[inline]
    pub fn try_eop_at(&self, jd_utc: JulianDate) -> Result<EopValues, EopError> {
        self.eop.try_eop_at(jd_utc)
    }

    /// Fallibly look up EOP values for a **TT** observation epoch.
    ///
    /// This converts TT to UTC before querying the provider, preserving the
    /// EOP provider contract while keeping transform call sites ergonomic.
    #[inline]
    pub fn try_eop_at_tt(&self, jd_tt: JulianDate) -> Result<EopValues, EopError> {
        let jd_utc = crate::astro::earth_rotation::try_jd_utc_from_tt(jd_tt)?;
        self.try_eop_at(jd_utc)
    }

    /// Look up EOP values for a **TT** observation epoch.
    #[inline]
    pub fn eop_at_tt(&self, jd_tt: JulianDate) -> EopValues {
        match crate::astro::earth_rotation::try_jd_utc_from_tt(jd_tt) {
            Ok(jd_utc) => self.eop_at(jd_utc),
            Err(_) => self
                .eop
                .eop_at(crate::astro::earth_rotation::jd_utc_from_tt_delta_t(jd_tt)),
        }
    }

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

use crate::ephemeris::DynEphemeris;

/// Astronomical context with a runtime-selected ephemeris backend.
///
/// This is the dynamic counterpart to [`AstroContext`]. Instead of selecting
/// an ephemeris backend at compile time via type parameters, it stores a
/// `Box<dyn DynEphemeris>` and dispatches via virtual calls.
///
/// Use this when:
/// - You load BSP files at runtime via [`RuntimeEphemeris`](crate::ephemeris::RuntimeEphemeris)
/// - You need to switch between backends without recompiling
/// - You want to override the default ephemeris at runtime
///
/// # Example
///
/// ```rust,ignore
/// use siderust::ephemeris::{RuntimeEphemeris, DynEphemeris};
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

    /// Fallibly look up EOP values for the given **UTC** Julian Date.
    #[inline]
    pub fn try_eop_at(&self, jd_utc: JulianDate) -> Result<EopValues, EopError> {
        self.eop.try_eop_at(jd_utc)
    }

    /// Fallibly look up EOP values for a **TT** observation epoch.
    #[inline]
    pub fn try_eop_at_tt(&self, jd_tt: JulianDate) -> Result<EopValues, EopError> {
        let jd_utc = crate::astro::earth_rotation::try_jd_utc_from_tt(jd_tt)?;
        self.try_eop_at(jd_utc)
    }

    /// Look up EOP values for a **TT** observation epoch.
    #[inline]
    pub fn eop_at_tt(&self, jd_tt: JulianDate) -> EopValues {
        match crate::astro::earth_rotation::try_jd_utc_from_tt(jd_tt) {
            Ok(jd_utc) => self.eop_at(jd_utc),
            Err(_) => self
                .eop
                .eop_at(crate::astro::earth_rotation::jd_utc_from_tt_delta_t(jd_tt)),
        }
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
        let jd = crate::J2000;
        let eop = ctx.eop_at_tt(jd);
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
        use crate::ephemeris::Vsop87Ephemeris;

        // Vsop87Ephemeris implements DynEphemeris via blanket impl
        let dyn_ctx = DynAstroContext::with_ephemeris(Box::new(Vsop87Ephemeris));
        let _eph = dyn_ctx.ephemeris();
        let _eop = dyn_ctx.eop_at(crate::J2000);
        let _eop_ref = dyn_ctx.eop();

        // Debug should work
        let debug_str = format!("{:?}", dyn_ctx);
        assert!(debug_str.contains("DynAstroContext"));
    }

    #[test]
    fn test_dyn_context_with_eop() {
        use crate::ephemeris::Vsop87Ephemeris;

        let dyn_ctx = DynAstroContext::with_ephemeris_and_eop(
            Box::new(Vsop87Ephemeris),
            DefaultEop::default(),
        );
        let eop = dyn_ctx.eop_at(crate::J2000);
        assert!(eop.dut1.is_finite());
    }
}
