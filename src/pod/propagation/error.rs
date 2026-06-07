// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error enum for `siderust-dynamics`.
//!
//! Wraps the two failure surfaces this crate exposes:
//!
//! * Upstream [`crate::astro::dynamics::errors::DynamicsError`] for any
//!   force-model / integrator / variational-equation failure.
//! * [`crate::pod::force::thrust::ManeuverError`] for finite-burn / thrust-arc input
//!   validation.
//!
//! Both upstream errors are funnelled here via [`From`] so consumer code
//! can use `?` against a single [`DynamicsError`] return type.

use thiserror::Error;

use principia::PrincipiaError;

use crate::pod::force::thrust::ManeuverError;

/// Unified error type for `siderust-dynamics`.
///
/// # Example
///
/// ```ignore
/// use siderust::pod::propagation::DynamicsError;
/// use siderust::pod::force::ManeuverError;
///
/// fn check(isp: f64) -> Result<(), DynamicsError> {
///     siderust::pod::force::mass_flow_rate(qtty::force::Newtons::new(1.0), isp)?;
///     Ok(())
/// }
/// assert!(matches!(check(0.0), Err(DynamicsError::Maneuver(ManeuverError::NonPositiveIsp(_)))));
/// ```
#[derive(Debug, Error)]
pub enum DynamicsError {
    /// Error originating from a force model, integrator, or variational
    /// propagator in upstream `siderust`.
    #[error(transparent)]
    Upstream(#[from] crate::astro::dynamics::errors::DynamicsError),

    /// Error originating from finite-burn / thrust-arc input validation.
    #[error(transparent)]
    Maneuver(#[from] ManeuverError),
}

impl From<PrincipiaError> for DynamicsError {
    fn from(err: PrincipiaError) -> Self {
        Self::Upstream(err.into())
    }
}

impl From<principia::PropagationError> for DynamicsError {
    fn from(err: principia::PropagationError) -> Self {
        let upstream: crate::astro::dynamics::errors::DynamicsError = err.into();
        Self::Upstream(upstream)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pod::force::{mass_flow_rate, ManeuverError};
    use principia::{PrincipiaError, PropagationError};
    use qtty::{force::Newtons, Second};

    #[test]
    fn maneuver_error_converts() {
        let err = DynamicsError::from(ManeuverError::NonPositiveIsp(0.0));
        assert!(matches!(
            err,
            DynamicsError::Maneuver(ManeuverError::NonPositiveIsp(_))
        ));
    }

    #[test]
    fn mass_flow_rate_failure_maps_to_maneuver_variant() {
        let err =
            DynamicsError::from(mass_flow_rate(Newtons::new(1.0), Second::new(0.0)).unwrap_err());
        assert!(matches!(err, DynamicsError::Maneuver(_)));
    }

    #[test]
    fn principia_error_converts_to_upstream() {
        let err = DynamicsError::from(PrincipiaError::DegenerateGeometry { reason: "test" });
        assert!(matches!(err, DynamicsError::Upstream(_)));
    }

    #[test]
    fn propagation_error_converts_to_upstream() {
        let err = DynamicsError::from(PropagationError::MaxStepsExceeded { max_steps: 1 });
        assert!(matches!(err, DynamicsError::Upstream(_)));
    }
}
