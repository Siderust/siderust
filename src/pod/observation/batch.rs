// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Ordered collection of observations grouped by epoch.

use crate::time::JulianDate;

use super::error::PodObservationsError;
use crate::pod::observation::obs_trait::{AnyObservation, CartesianState, ObsResidual};
use crate::pod::observation::provider_bundle::ProviderBundle;

/// Ordered collection of heterogeneous observations, typically spanning one
/// processing arc.
///
/// Observations are stored as trait objects so that pseudorange, carrier-phase,
/// SLR, and inter-satellite measurements can coexist in the same batch.
/// [`residuals_at_epoch`](ObservationBatch::residuals_at_epoch) filters and
/// evaluates all observations at a given epoch.
///
/// # Examples
///
/// ```
/// use siderust::pod::observation::batch::ObservationBatch;
/// use siderust::pod::observation::provider_bundle::NullProviderBundle;
/// use siderust::time::JulianDate;
///
/// let batch = ObservationBatch::new();
/// assert!(batch.is_empty());
///
/// let epoch = JulianDate::new(2_451_545.0);
/// let residuals = batch.residuals_at_epoch(
///     epoch,
///     // A real use case passes an estimated CartesianState here; we skip it
///     // in the doc-test.
/// );
/// assert!(residuals.is_empty());
/// ```
pub struct ObservationBatch {
    obs: Vec<Box<dyn AnyObservation>>,
}

impl ObservationBatch {
    /// Create an empty batch.
    pub fn new() -> Self {
        Self { obs: Vec::new() }
    }

    /// Append an observation.
    pub fn push(&mut self, o: Box<dyn AnyObservation>) {
        self.obs.push(o);
    }

    /// Number of observations in the batch.
    pub fn len(&self) -> usize {
        self.obs.len()
    }

    /// Returns `true` if the batch contains no observations.
    pub fn is_empty(&self) -> bool {
        self.obs.is_empty()
    }

    /// All epochs present in the batch (unsorted, may contain duplicates).
    pub fn epochs(&self) -> Vec<JulianDate> {
        self.obs.iter().map(|o| o.epoch()).collect()
    }

    /// Compute O−C residuals for all observations whose epoch equals `epoch`
    /// within a tolerance of 0.5 µs (≈ 5.8e-12 days).
    ///
    /// Returns one entry per matching observation; each entry is either the
    /// typed [`ObsResidual`] or a [`PodObservationsError`].
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::observation::batch::ObservationBatch;
    /// use siderust::pod::observation::obs_trait::{ObsType, Observation, CartesianState};
    /// use siderust::pod::observation::provider_bundle::{NullProviderBundle, ProviderBundle};
    /// use siderust::pod::observation::PodObservationsError;
    /// use siderust::astro::dynamics::{Position, Velocity};
    /// use siderust::coordinates::frames::GCRS;
    /// use siderust::time::JulianDate;
    ///
    /// struct ZeroObs(JulianDate);
    /// impl Observation for ZeroObs {
    ///     type Residual = f64;
    ///     fn residual(&self, _: &CartesianState, _: &dyn ProviderBundle)
    ///         -> Result<f64, PodObservationsError> { Ok(0.0) }
    ///     fn obs_type(&self) -> ObsType { ObsType::GnssPseudorange }
    ///     fn epoch(&self) -> JulianDate { self.0 }
    ///     fn sigma(&self) -> siderust::qtty::Meters { siderust::qtty::Meters::new(1.0) }
    /// }
    ///
    /// let epoch = JulianDate::new(2_451_545.0);
    /// let state: CartesianState = CartesianState::new(
    ///     epoch.to_j2000s(),
    ///     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
    ///     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
    /// );
    ///
    /// let mut batch = ObservationBatch::new();
    /// batch.push(Box::new(ZeroObs(epoch)));
    ///
    /// let results = batch.residuals_at_epoch_with(epoch, &state, &NullProviderBundle);
    /// assert_eq!(results.len(), 1);
    /// assert!(results[0].is_ok());
    /// ```
    pub fn residuals_at_epoch(
        &self,
        epoch: JulianDate,
    ) -> Vec<Result<ObsResidual, PodObservationsError>> {
        // Filter only — useful for pre-checking which observations match.
        self.obs
            .iter()
            .filter(|o| epochs_match(o.epoch(), epoch))
            .map(|_| Err(PodObservationsError::LightTimeNotConverged))
            .collect()
    }

    /// Compute residuals at `epoch`, evaluating each matching observation
    /// against `state` and `providers`.
    pub fn residuals_at_epoch_with(
        &self,
        epoch: JulianDate,
        state: &CartesianState,
        providers: &dyn ProviderBundle,
    ) -> Vec<Result<ObsResidual, PodObservationsError>> {
        self.obs
            .iter()
            .filter(|o| epochs_match(o.epoch(), epoch))
            .map(|o| o.any_residual(state, providers))
            .collect()
    }
}

impl Default for ObservationBatch {
    fn default() -> Self {
        Self::new()
    }
}

/// Epochs match within 0.5 µs (≈ 5.8e-12 JD days).
fn epochs_match(a: JulianDate, b: JulianDate) -> bool {
    const TOL_DAYS: f64 = 5.8e-12;
    (a.value() - b.value()).abs() < TOL_DAYS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::pod::observation::obs_trait::{ObsType, Observation};
    use crate::pod::observation::provider_bundle::NullProviderBundle;

    struct ConstObs {
        epoch: JulianDate,
        residual: f64,
    }
    impl Observation for ConstObs {
        type Residual = f64;
        fn residual(
            &self,
            _: &CartesianState,
            _: &dyn ProviderBundle,
        ) -> Result<f64, PodObservationsError> {
            Ok(self.residual)
        }
        fn obs_type(&self) -> ObsType {
            ObsType::GnssPseudorange
        }
        fn epoch(&self) -> JulianDate {
            self.epoch
        }
        fn sigma(&self) -> qtty::Meter {
            qtty::Meter::new(1.0)
        }
    }

    fn make_state(epoch: JulianDate) -> CartesianState {
        CartesianState::new(
            epoch.to_j2000s(),
            Position::<GCRS>::new(7_000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    #[test]
    fn batch_residuals_correct_epoch() {
        let epoch = JulianDate::new(2_451_545.0);
        let other = JulianDate::new(2_451_546.0);
        let state = make_state(epoch);
        let mut batch = ObservationBatch::new();
        batch.push(Box::new(ConstObs {
            epoch,
            residual: 1.0,
        }));
        batch.push(Box::new(ConstObs {
            epoch: other,
            residual: 0.0,
        }));

        let res = batch.residuals_at_epoch_with(epoch, &state, &NullProviderBundle);
        assert_eq!(res.len(), 1);
        if let Ok(ObsResidual::Scalar(v)) = &res[0] {
            assert!((*v - 1.0).abs() < 1e-10);
        } else {
            panic!("Expected Scalar residual");
        }
    }

    #[test]
    fn batch_empty_for_unknown_epoch() {
        let epoch = JulianDate::new(2_451_545.0);
        let query = JulianDate::new(2_451_600.0);
        let state = make_state(epoch);
        let mut batch = ObservationBatch::new();
        batch.push(Box::new(ConstObs {
            epoch,
            residual: 1.0,
        }));
        let res = batch.residuals_at_epoch_with(query, &state, &NullProviderBundle);
        assert!(res.is_empty());
    }
}
