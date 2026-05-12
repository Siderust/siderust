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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::errors::DynamicsError;

    fn dummy_err() -> DynamicsError {
        DynamicsError::DegenerateGeometry { reason: "test" }
    }

    #[test]
    fn step_control_display() {
        let e = PropagationError::StepControl(dummy_err());
        let s = format!("{e}");
        assert!(s.contains("integrator step control"), "got: {s}");
    }

    #[test]
    fn step_below_minimum_display() {
        let e = PropagationError::StepBelowMinimum { h_requested: 1e-10, h_min: 1e-6 };
        let s = format!("{e}");
        assert!(s.contains("minimum"), "got: {s}");
    }

    #[test]
    fn max_steps_exceeded_display() {
        let e = PropagationError::MaxStepsExceeded { max_steps: 1000 };
        let s = format!("{e}");
        assert!(s.contains("1000"), "got: {s}");
    }

    #[test]
    fn event_evaluation_display() {
        let e = PropagationError::EventEvaluation {
            name: "altitude",
            source: dummy_err(),
        };
        let s = format!("{e}");
        assert!(s.contains("altitude"), "got: {s}");
    }

    #[test]
    fn step_control_has_source() {
        use std::error::Error;
        let e = PropagationError::StepControl(dummy_err());
        assert!(e.source().is_some());
    }

    #[test]
    fn event_evaluation_has_source() {
        use std::error::Error;
        let e = PropagationError::EventEvaluation {
            name: "test",
            source: dummy_err(),
        };
        assert!(e.source().is_some());
    }

    #[test]
    fn step_below_minimum_no_source() {
        use std::error::Error;
        let e = PropagationError::StepBelowMinimum { h_requested: 1e-10, h_min: 1e-6 };
        assert!(e.source().is_none());
    }

    #[test]
    fn max_steps_no_source() {
        use std::error::Error;
        let e = PropagationError::MaxStepsExceeded { max_steps: 100 };
        assert!(e.source().is_none());
    }

    #[test]
    fn from_propagation_error_for_dynamics_error() {
        let e = PropagationError::MaxStepsExceeded { max_steps: 42 };
        let d: DynamicsError = DynamicsError::from(e);
        let s = format!("{d:?}");
        assert!(!s.is_empty());
    }
}
