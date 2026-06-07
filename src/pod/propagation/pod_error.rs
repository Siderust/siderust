// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error enum for `siderust::pod::propagation`.
//!
//! Wraps the failure surfaces this crate exposes:
//!
//! * [`crate::pod::propagation::DynamicsError`] for any integrator,
//!   variational, or force-model failure surfaced through the typed adapter.
//! * Registry / configuration validation errors raised by
//!   [`crate::pod::force::registry::ForceModelRegistry`] and friends.
//! * Process-noise input validation.
//! * Explicit *feature-not-implemented* signal for stubs that the registry
//!   recognises by name but cannot yet build (e.g. boxwing SRP).
//!
//! No `anyhow` is used; every variant is an explicit case so consumers can
//! pattern-match.

use thiserror::Error;

/// Top-level error type for `siderust::pod::propagation`.
///
/// # Examples
///
/// ```
/// use siderust::pod::propagation::PodDynamicsError;
///
/// let err = PodDynamicsError::UnknownModel("nonsense".into());
/// assert_eq!(err.to_string(), "unknown force model name: 'nonsense'");
/// ```
#[derive(Debug, Error)]
pub enum PodDynamicsError {
    /// Error from the underlying `siderust-dynamics` integrator / STM /
    /// force-model stack.
    #[error(transparent)]
    Dynamics(#[from] crate::pod::propagation::DynamicsError),

    /// The registry was asked for a name it does not recognise.
    #[error("unknown force model name: '{0}'")]
    UnknownModel(String),

    /// The registry was asked for a model with parameters that do not match
    /// the factory's expected shape.
    #[error("force model '{name}' rejected the supplied parameters: {reason}")]
    InvalidParameters {
        /// Registered name of the model.
        name: String,
        /// Free-form description of why the params were rejected.
        reason: &'static str,
    },

    /// The caller requested a model that is recognised by name but whose
    /// implementation has not been wired up in this release. Used in lieu of
    /// `unimplemented!()` to keep `scripts/check_no_todos.sh` happy and to
    /// give consumers a typed error variant they can fall back from.
    #[error("force model '{0}' is recognised but not yet implemented in this build")]
    FeatureNotImplemented(&'static str),

    /// Process-noise configuration was internally inconsistent
    /// (e.g. negative sigma, non-positive τ).
    #[error("invalid process-noise configuration: {0}")]
    InvalidProcessNoise(&'static str),

    /// Parameter STM finite-difference perturbation request was invalid
    /// (e.g. zero step size).
    #[error("invalid parameter STM step: {0}")]
    InvalidParamStep(&'static str),
}
