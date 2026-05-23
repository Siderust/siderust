// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error enum for `siderust-dynamics`.
//!
//! Wraps the two failure surfaces this crate exposes:
//!
//! * Upstream [`siderust::astro::dynamics::errors::DynamicsError`] for any
//!   force-model / integrator / variational-equation failure.
//! * [`siderust::pod::dynamics::thrust::ManeuverError`] for finite-burn / thrust-arc input
//!   validation.
//!
//! Both upstream errors are funnelled here via [`From`] so consumer code
//! can use `?` against a single [`DynamicsError`] return type.

use thiserror::Error;

use principia::PrincipiaError;

use super::thrust::ManeuverError;

/// Unified error type for `siderust-dynamics`.
///
/// # Example
///
/// ```ignore
/// use siderust::pod::dynamics::{DynamicsError, ManeuverError};
///
/// fn check(isp: f64) -> Result<(), DynamicsError> {
///     siderust::pod::dynamics::mass_flow_rate(qtty::force::Newtons::new(1.0), isp)?;
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
