// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Core observation trait, type aliases, and residual types.

use crate::astro::dynamics::OrbitState;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::time::JulianDate;

use super::error::PodObservationsError;
use super::provider_bundle::ProviderBundle;

// ─── Type aliases ────────────────────────────────────────────────────────────

/// Geocentric inertial Cartesian orbit state, position + velocity in GCRS.
///
/// Equivalent to `OrbitState<Geocentric, GCRS>` from `siderust`.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::CartesianState;
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// let state: CartesianState = CartesianState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
/// );
/// assert!((state.position.x().value() - 7_000.0).abs() < 1e-10);
/// ```
pub type CartesianState = OrbitState<Geocentric, GCRS>;

// ─── Sensor taxonomy ─────────────────────────────────────────────────────────

/// Observable type tag.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::ObsType;
/// assert_ne!(ObsType::GnssPseudorange, ObsType::SlrNormalPoint);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ObsType {
    /// GNSS pseudorange (code-based one-way range, metres).
    GnssPseudorange,
    /// GNSS carrier-phase range (phase-based, metres + integer ambiguity).
    GnssCarrierPhase,
    /// SLR two-way normal-point range (metres).
    SlrNormalPoint,
    /// Inter-satellite range (metres), e.g. LISA arm lengths.
    InterSatRange,
    /// Doppler/range-rate (m/s).
    Doppler,
}

// ─── Residual types ──────────────────────────────────────────────────────────

/// Carrier-phase residual carrying both the metric residual and the fractional
/// cycle count.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::PhaseResidual;
///
/// let r = PhaseResidual { residual_m: 0.5, cycles: 0.5 / 0.1903 };
/// assert!((r.residual_m).abs() < 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseResidual {
    /// Residual (observed − modelled) in metres.
    pub residual_m: f64,
    /// Residual expressed in carrier cycles (residual_m / wavelength_m).
    pub cycles: f64,
}

/// Type-erased residual for use in heterogeneous [`crate::pod::observation::ObservationBatch`]
/// collections.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::{ObsResidual, PhaseResidual};
///
/// let r = ObsResidual::Scalar(3.0);
/// let p = ObsResidual::Phase(PhaseResidual { residual_m: 0.1, cycles: 0.5 });
/// assert!(matches!(r, ObsResidual::Scalar(_)));
/// assert!(matches!(p, ObsResidual::Phase(_)));
/// ```
#[derive(Debug, Clone)]
pub enum ObsResidual {
    /// Scalar residual (pseudorange, SLR, Doppler) in metres.
    Scalar(f64),
    /// Carrier-phase residual.
    Phase(PhaseResidual),
}

impl From<f64> for ObsResidual {
    fn from(v: f64) -> Self {
        ObsResidual::Scalar(v)
    }
}

impl From<PhaseResidual> for ObsResidual {
    fn from(p: PhaseResidual) -> Self {
        ObsResidual::Phase(p)
    }
}

// ─── Observation trait ───────────────────────────────────────────────────────

/// Typed observation model that connects a tracked spacecraft state to one
/// observable.
///
/// Implementors embed all measurement-specific data (observed value, noise
/// model, atmospheric parameters, satellite ephemeris, …) and expose a single
/// entry point that computes the O−C residual given the current estimated
/// state and an optional provider bundle for dynamic look-ups (clock biases,
/// station positions, …).
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::{CartesianState, Observation, ObsType};
/// use siderust::pod::observation::provider_bundle::NullProviderBundle;
/// use siderust::pod::observation::PodObservationsError;
/// use siderust::astro::dynamics::{Position, Velocity};
/// use siderust::coordinates::frames::GCRS;
/// use siderust::time::JulianDate;
///
/// struct ConstantObs;
/// impl Observation for ConstantObs {
///     type Residual = f64;
///     fn modeled_value(
///         &self,
///         _state: &CartesianState,
///         _providers: &dyn siderust::pod::observation::provider_bundle::ProviderBundle,
///     ) -> Result<f64, PodObservationsError> {
///         Ok(0.0)
///     }
///     fn obs_type(&self) -> ObsType { ObsType::GnssPseudorange }
///     fn epoch(&self) -> JulianDate { JulianDate::new(2_451_545.0) }
///     fn sigma(&self) -> siderust::qtty::Meters { siderust::qtty::Meters::new(1.0) }
/// }
/// let state: CartesianState = CartesianState::new(
///     JulianDate::new(2_451_545.0).to_j2000s(),
///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
///     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
/// );
/// let r = ConstantObs.modeled_value(&state, &NullProviderBundle).unwrap();
/// assert_eq!(r, 0.0);
/// ```
pub trait Observation: Send + Sync {
    /// Residual type (e.g. `f64` for range, [`PhaseResidual`] for carrier-phase).
    type Residual;

    /// Compute the O−C residual at the given estimated state.
    ///
    /// The `providers` argument is consulted for dynamic auxiliary data such as
    /// satellite clock biases, station positions, and receiver clock biases.
    fn modeled_value(
        &self,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Result<Self::Residual, PodObservationsError>;

    /// Observable type (sensor taxonomy).
    fn obs_type(&self) -> ObsType;

    /// Epoch at which the observation was collected (TT Julian date).
    fn epoch(&self) -> JulianDate;

    /// Assumed measurement standard deviation (same units as the residual).
    fn sigma(&self) -> qtty::Meter;
}

// ─── Object-safe erasure ─────────────────────────────────────────────────────

/// Object-safe version of [`Observation`] used inside
/// [`crate::pod::observation::ObservationBatch`].
///
/// This trait is blanket-implemented for every `T: Observation` whose
/// `Residual` type can be converted into [`ObsResidual`].
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::obs_trait::{AnyObservation, ObsResidual};
/// // See ObservationBatch for a usage example.
/// let _: &dyn AnyObservation;
/// ```
pub trait AnyObservation: Send + Sync {
    /// Compute the type-erased O−C residual.
    fn any_modeled_value(
        &self,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Result<ObsResidual, PodObservationsError>;

    /// Observable type.
    fn obs_type(&self) -> ObsType;

    /// Epoch.
    fn epoch(&self) -> JulianDate;

    /// Standard deviation.
    fn sigma(&self) -> qtty::Meter;
}

impl<T> AnyObservation for T
where
    T: Observation + 'static,
    ObsResidual: From<T::Residual>,
{
    fn any_modeled_value(
        &self,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Result<ObsResidual, PodObservationsError> {
        self.modeled_value(state, providers).map(ObsResidual::from)
    }

    fn obs_type(&self) -> ObsType {
        Observation::obs_type(self)
    }

    fn epoch(&self) -> JulianDate {
        Observation::epoch(self)
    }

    fn sigma(&self) -> qtty::Meter {
        Observation::sigma(self)
    }
}
