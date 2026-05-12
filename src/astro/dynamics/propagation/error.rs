// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Errors specific to the high-level propagation driver.
//!
//! ## Scope
//!
//! Provides [`PropagationError`], which wraps integrator failures, step-control
//! issues, event evaluation errors, and budget violations into a unified error type.

use crate::astro::dynamics::errors::DynamicsError;
use std::fmt;

#[derive(Debug)]
pub enum PropagationError {
    /// Wrap an underlying integrator step error.
    StepControl(DynamicsError),
    /// The controller would have shrunk the step below the configured minimum.
    StepBelowMinimum { h_requested: f64, h_min: f64 },
    /// The driver exceeded its configured step budget.
    MaxStepsExceeded { max_steps: u32 },
    /// An event detector returned an error while evaluating its switching function.
    EventEvaluation {
        name: &'static str,
        source: DynamicsError,
    },
}

impl fmt::Display for PropagationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StepControl(e) => write!(f, "integrator step control error: {e}"),
            Self::StepBelowMinimum { h_requested, h_min } => write!(
                f,
                "requested step {h_requested:e} s falls below configured minimum {h_min:e} s"
            ),
            Self::MaxStepsExceeded { max_steps } => {
                write!(f, "propagation exceeded max_steps={max_steps}")
            }
            Self::EventEvaluation { name, source } => {
                write!(f, "event '{name}' evaluation failed: {source}")
            }
        }
    }
}

impl std::error::Error for PropagationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::StepControl(e) => Some(e),
            Self::EventEvaluation { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl From<PropagationError> for DynamicsError {
    fn from(e: PropagationError) -> Self {
        DynamicsError::Provider(Box::new(e))
    }
}
