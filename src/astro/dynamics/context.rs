// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime context threaded through every dynamics call.
//!
//! ## Design
//!
//! [`DynamicsContext`] is passed *by reference* at every `acceleration()` /
//! `partials()` call site.  **Force models hold no providers**; they hold only
//! their tunable physical parameters (Cd, Cr, A/m, gravity truncation degree,
//! empirical coefficients, …).  The context is the single location where
//! providers are resolved.
//!
//! This separation means a force model can be constructed once and re-used
//! across different provider configurations (e.g. a Vsop87 context in tests, a
//! DE440 context in production) without any code duplication.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::astro::dynamics::context::{DynamicsContext, DynamicsContextBuilder};
//!
//! // Minimal context with no providers (useful for two-body propagation).
//! let ctx = DynamicsContext::empty();
//! assert!(ctx.require_ephemeris().is_err());
//! assert!(ctx.require_atmosphere().is_err());
//! assert!(ctx.require_gravity_field().is_err());
//!
//! // Builder with explicit providers.
//! let ctx = DynamicsContextBuilder::new()
//!     // .with_ephemeris(Arc::new(my_ephemeris))
//!     .build();
//! assert!(ctx.require_ephemeris().is_err());
//! ```

use std::sync::Arc;

use crate::astro::eop::EopValues;
use crate::ephemeris::DynEphemeris;
use crate::time::JulianDate;

use super::density::DensityProvider;
use super::errors::DynamicsError;
use super::gravity::GravityFieldProvider;

// =============================================================================
// Placeholder provider traits
// =============================================================================

/// Earth Orientation Parameters provider.
///
/// Interface for supplying Earth Orientation Parameters (EOP) to dynamics
/// computations.  The trait exists so [`DynamicsContext`] can carry an optional
/// `Arc<dyn EarthOrientationProvider>` without creating a hard dependency on
/// any particular EOP backend.
///
/// Downstream crates that supply real EOP data should implement this trait and
/// inject the implementation via [`DynamicsContextBuilder::with_eop`].
///
/// # Default implementations
///
/// All methods provide default implementations that return `None`, so existing
/// marker implementations continue to compile unchanged.  Override only the
/// methods relevant to your provider.
pub trait EarthOrientationProvider: Send + Sync {
    /// EOP values at the given UTC Julian Date.
    ///
    /// Returns `None` if data is unavailable for the requested epoch.
    ///
    /// # Time-scale contract
    /// `jd_utc` **must** be a UTC Julian Date. Passing a TT or UT1 Julian
    /// Date will silently return incorrect (or no) EOP values.
    fn eop_at(&self, _jd_utc: JulianDate) -> Option<EopValues> {
        None
    }
}

/// Solar activity (F10.7 flux, Ap) provider.
///
/// Interface for supplying space-weather indices to atmospheric density models
/// such as NRLMSISE-00 or JB2008.  The trait exists so [`DynamicsContext`] can
/// carry an optional `Arc<dyn SolarActivityProvider>` without creating a hard
/// dependency on any particular space-weather backend.
///
/// Downstream crates that supply real solar-activity indices should implement
/// this trait and inject the implementation via
/// [`DynamicsContextBuilder::with_solar_activity`].
///
/// # Default implementations
///
/// All methods provide default implementations that return `None`, so existing
/// marker implementations continue to compile unchanged.  Override only the
/// methods relevant to your provider.
pub trait SolarActivityProvider: Send + Sync {
    /// F10.7 solar flux index (solar flux units, SFU) at the given UTC Julian Date.
    ///
    /// Returns `None` if data is unavailable for the requested epoch.
    fn f107_sfu(&self, _jd_utc: JulianDate) -> Option<f64> {
        None
    }

    /// Daily Ap geomagnetic index at the given UTC Julian Date.
    ///
    /// Returns `None` if data is unavailable for the requested epoch.
    fn ap_daily(&self, _jd_utc: JulianDate) -> Option<f64> {
        None
    }
}

// =============================================================================
// Conventions
// =============================================================================

/// IAU precession model selection.
///
/// Controls which precession/nutation model is used when transforming between
/// inertial and Earth-fixed frames.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PrecessionModel {
    /// IAU 2006/2000A precession-nutation model (recommended for modern work).
    #[default]
    IAU2006,

    /// IAU 2000 precession-nutation model.
    IAU2000,
}

/// Hint for the default time scale used in dynamics computations.
///
/// This is a placeholder for future work; concrete time-scale semantics are
/// managed by [`tempoch`](https://docs.rs/tempoch) typed instants.  The hint
/// allows a context to advertise its preference without enforcing it at the
/// type level today.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeScaleHint {
    /// Terrestrial Time (TT) — the preferred scale for dynamics.
    #[default]
    TT,

    /// Geocentric Coordinate Time (TCG).
    TCG,

    /// Barycentric Dynamical Time (TDB).
    TDB,
}

/// Convention metadata bundled into a [`DynamicsContext`].
///
/// Collects choices that affect numerical results but are not tied to a
/// single force model (e.g. which IAU precession model to use, what time
/// scale computations default to).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Conventions {
    /// Precession/nutation model to use for frame transforms.
    pub iau_precession_model: PrecessionModel,

    /// Default time scale for dynamics computations.
    pub time_scale_default: TimeScaleHint,
}

impl Default for Conventions {
    /// Returns IAU 2006 precession with TT as the default time scale.
    fn default() -> Self {
        Self {
            iau_precession_model: PrecessionModel::IAU2006,
            time_scale_default: TimeScaleHint::TT,
        }
    }
}

// =============================================================================
// DynamicsContext
// =============================================================================

/// Runtime context passed by reference to every force-model call.
///
/// All provider slots are optional [`Arc`]-wrapped trait objects.  Use
/// [`DynamicsContext::empty()`] when no providers are needed (e.g. pure
/// two-body propagation) or [`DynamicsContext::builder()`] to construct a
/// fully populated context.
///
/// Convenience accessors ([`require_ephemeris`][Self::require_ephemeris],
/// [`require_atmosphere`][Self::require_atmosphere],
/// [`require_gravity_field`][Self::require_gravity_field]) return a typed
/// [`DynamicsError`] when the requested provider is absent, allowing force
/// models to propagate the error cleanly without `unwrap`.
pub struct DynamicsContext {
    /// Ephemeris provider for third-body and SRP computations.
    pub ephemeris: Option<Arc<dyn DynEphemeris + Send + Sync>>,

    /// Atmospheric density provider for drag force models.
    pub atmosphere: Option<Arc<dyn DensityProvider + Send + Sync>>,

    /// Geopotential field provider for high-fidelity gravity models.
    pub gravity: Option<Arc<dyn GravityFieldProvider + Send + Sync>>,

    /// Earth Orientation Parameters provider.
    pub eop: Option<Arc<dyn EarthOrientationProvider + Send + Sync>>,

    /// Solar activity provider for atmosphere-density and SRP modelling.
    pub solar_activity: Option<Arc<dyn SolarActivityProvider + Send + Sync>>,

    /// Convention choices (precession model, default time scale, …).
    pub conventions: Conventions,
}

impl DynamicsContext {
    /// Construct an empty context with no providers and default conventions.
    ///
    /// Useful for two-body propagation or unit tests that do not need real
    /// provider data.
    pub fn empty() -> Self {
        Self {
            ephemeris: None,
            atmosphere: None,
            gravity: None,
            eop: None,
            solar_activity: None,
            conventions: Conventions::default(),
        }
    }

    /// Start building a [`DynamicsContext`] with an explicit provider set.
    pub fn builder() -> DynamicsContextBuilder {
        DynamicsContextBuilder::new()
    }

    // ---- Convenience accessors ----

    /// Return a reference to the ephemeris provider, or an
    /// [`DynamicsError::EphemerisUnavailable`] error if none is configured.
    pub fn require_ephemeris(&self) -> Result<&Arc<dyn DynEphemeris + Send + Sync>, DynamicsError> {
        self.ephemeris
            .as_ref()
            .ok_or(DynamicsError::EphemerisUnavailable {
                body: "(any)",
                source: None,
            })
    }

    /// Return a reference to the atmosphere density provider, or an
    /// [`DynamicsError::AtmosphereProviderError`]-wrapped error if none is
    /// configured.
    ///
    /// The error is wrapped in a minimal `io::Error` to keep the
    /// `Box<dyn Error + Send + Sync>` type requirement satisfied.
    pub fn require_atmosphere(
        &self,
    ) -> Result<&Arc<dyn DensityProvider + Send + Sync>, DynamicsError> {
        self.atmosphere.as_ref().ok_or_else(|| {
            DynamicsError::AtmosphereProviderError(Box::new(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "no atmosphere provider in DynamicsContext",
            )))
        })
    }

    /// Return a reference to the gravity field provider, or
    /// [`DynamicsError::GravityFieldUnavailable`] if none is configured.
    pub fn require_gravity_field(
        &self,
    ) -> Result<&Arc<dyn GravityFieldProvider + Send + Sync>, DynamicsError> {
        self.gravity
            .as_ref()
            .ok_or(DynamicsError::GravityFieldUnavailable)
    }
}

// =============================================================================
// DynamicsContextBuilder
// =============================================================================

/// Builder for [`DynamicsContext`].
///
/// Construct via [`DynamicsContext::builder()`] or [`DynamicsContextBuilder::new()`].
/// Call the `with_*` methods to set providers, then call [`build`][Self::build]
/// to obtain the finished context.
#[derive(Default)]
pub struct DynamicsContextBuilder {
    ephemeris: Option<Arc<dyn DynEphemeris + Send + Sync>>,
    atmosphere: Option<Arc<dyn DensityProvider + Send + Sync>>,
    gravity: Option<Arc<dyn GravityFieldProvider + Send + Sync>>,
    eop: Option<Arc<dyn EarthOrientationProvider + Send + Sync>>,
    solar_activity: Option<Arc<dyn SolarActivityProvider + Send + Sync>>,
    conventions: Conventions,
}

impl DynamicsContextBuilder {
    /// Create a builder with no providers and default conventions.
    pub fn new() -> Self {
        Self {
            ephemeris: None,
            atmosphere: None,
            gravity: None,
            eop: None,
            solar_activity: None,
            conventions: Conventions::default(),
        }
    }

    /// Set the ephemeris provider.
    pub fn with_ephemeris(mut self, e: Arc<dyn DynEphemeris + Send + Sync>) -> Self {
        self.ephemeris = Some(e);
        self
    }

    /// Set the atmosphere density provider.
    pub fn with_atmosphere(mut self, a: Arc<dyn DensityProvider + Send + Sync>) -> Self {
        self.atmosphere = Some(a);
        self
    }

    /// Set the gravity field provider.
    pub fn with_gravity(mut self, g: Arc<dyn GravityFieldProvider + Send + Sync>) -> Self {
        self.gravity = Some(g);
        self
    }

    /// Set the Earth Orientation Parameters provider.
    pub fn with_eop(mut self, e: Arc<dyn EarthOrientationProvider + Send + Sync>) -> Self {
        self.eop = Some(e);
        self
    }

    /// Set the solar activity provider.
    pub fn with_solar_activity(mut self, s: Arc<dyn SolarActivityProvider + Send + Sync>) -> Self {
        self.solar_activity = Some(s);
        self
    }

    /// Override the convention metadata.
    pub fn with_conventions(mut self, c: Conventions) -> Self {
        self.conventions = c;
        self
    }

    /// Consume the builder and produce a [`DynamicsContext`].
    pub fn build(self) -> DynamicsContext {
        DynamicsContext {
            ephemeris: self.ephemeris,
            atmosphere: self.atmosphere,
            gravity: self.gravity,
            eop: self.eop,
            solar_activity: self.solar_activity,
            conventions: self.conventions,
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::density::ExponentialAtmosphere;
    use crate::astro::dynamics::gravity::TwoBodyEarth;

    // ---- empty() ----

    #[test]
    fn empty_has_no_providers() {
        let ctx = DynamicsContext::empty();
        assert!(ctx.ephemeris.is_none());
        assert!(ctx.atmosphere.is_none());
        assert!(ctx.gravity.is_none());
        assert!(ctx.eop.is_none());
        assert!(ctx.solar_activity.is_none());
    }

    #[test]
    fn empty_uses_default_conventions() {
        let ctx = DynamicsContext::empty();
        assert_eq!(
            ctx.conventions.iau_precession_model,
            PrecessionModel::IAU2006
        );
        assert_eq!(ctx.conventions.time_scale_default, TimeScaleHint::TT);
    }

    // ---- builder roundtrip ----

    #[test]
    fn builder_with_gravity_roundtrip() {
        let g: Arc<dyn GravityFieldProvider + Send + Sync> = Arc::new(TwoBodyEarth);
        let ctx = DynamicsContextBuilder::new().with_gravity(g).build();
        assert!(ctx.gravity.is_some());
        assert!(ctx.ephemeris.is_none());
        assert!(ctx.atmosphere.is_none());
    }

    #[test]
    fn builder_with_atmosphere_roundtrip() {
        let a: Arc<dyn DensityProvider + Send + Sync> = Arc::new(ExponentialAtmosphere::LEO_500KM);
        let ctx = DynamicsContextBuilder::new().with_atmosphere(a).build();
        assert!(ctx.atmosphere.is_some());
        assert!(ctx.ephemeris.is_none());
        assert!(ctx.gravity.is_none());
    }

    #[test]
    fn builder_with_conventions() {
        let c = Conventions {
            iau_precession_model: PrecessionModel::IAU2000,
            time_scale_default: TimeScaleHint::TDB,
        };
        let ctx = DynamicsContextBuilder::new().with_conventions(c).build();
        assert_eq!(
            ctx.conventions.iau_precession_model,
            PrecessionModel::IAU2000
        );
        assert_eq!(ctx.conventions.time_scale_default, TimeScaleHint::TDB);
    }

    #[test]
    fn context_builder_shorthand() {
        let ctx = DynamicsContext::builder().build();
        assert!(ctx.ephemeris.is_none());
    }

    // ---- require_* accessors ----

    #[test]
    fn require_ephemeris_returns_error_when_absent() {
        let ctx = DynamicsContext::empty();
        let result = ctx.require_ephemeris();
        assert!(result.is_err());
        assert!(matches!(
            result.err().unwrap(),
            DynamicsError::EphemerisUnavailable { .. }
        ));
    }

    #[test]
    fn require_atmosphere_returns_error_when_absent() {
        let ctx = DynamicsContext::empty();
        let result = ctx.require_atmosphere();
        assert!(result.is_err());
        assert!(matches!(
            result.err().unwrap(),
            DynamicsError::AtmosphereProviderError(_)
        ));
    }

    #[test]
    fn require_gravity_field_returns_error_when_absent() {
        let ctx = DynamicsContext::empty();
        let result = ctx.require_gravity_field();
        assert!(result.is_err());
        assert!(matches!(
            result.err().unwrap(),
            DynamicsError::GravityFieldUnavailable
        ));
    }

    #[test]
    fn require_gravity_field_returns_ok_when_present() {
        let g: Arc<dyn GravityFieldProvider + Send + Sync> = Arc::new(TwoBodyEarth);
        let ctx = DynamicsContextBuilder::new().with_gravity(g).build();
        assert!(ctx.require_gravity_field().is_ok());
    }

    #[test]
    fn require_atmosphere_returns_ok_when_present() {
        let a: Arc<dyn DensityProvider + Send + Sync> = Arc::new(ExponentialAtmosphere::LEO_500KM);
        let ctx = DynamicsContextBuilder::new().with_atmosphere(a).build();
        assert!(ctx.require_atmosphere().is_ok());
    }
}
