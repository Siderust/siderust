// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Errors for surface-observer construction and transforms.

use core::fmt;

/// Error returned by surface-site construction or body-surface transforms.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceError {
    /// Semi-major axis, radius, or height is non-finite or non-positive where required.
    InvalidDimension {
        /// Which parameter failed validation.
        parameter: &'static str,
        /// Offending value.
        value: f64,
    },
    /// Flattening is outside the open interval `(0, 1)`.
    InvalidFlattening {
        /// Offending flattening.
        value: f64,
    },
    /// Latitude is outside `[-90°, +90°]` after normalisation attempts failed.
    InvalidLatitude {
        /// Offending latitude in degrees.
        degrees: f64,
    },
    /// The selected body has no documented solid/reference surface model.
    UnsupportedSurface {
        /// Body name.
        body: &'static str,
        /// Short reason.
        reason: &'static str,
    },
    /// Degenerate geometry (zero range, undefined azimuth at pole, etc.).
    DegenerateGeometry {
        /// Short reason.
        reason: &'static str,
    },
    /// Non-finite input component.
    NonFiniteInput {
        /// Which field failed.
        field: &'static str,
    },
}

impl fmt::Display for SurfaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimension { parameter, value } => {
                write!(f, "invalid surface dimension {parameter}={value}")
            }
            Self::InvalidFlattening { value } => {
                write!(f, "invalid ellipsoid flattening f={value}")
            }
            Self::InvalidLatitude { degrees } => {
                write!(f, "invalid latitude {degrees}°")
            }
            Self::UnsupportedSurface { body, reason } => {
                write!(f, "unsupported surface model for {body}: {reason}")
            }
            Self::DegenerateGeometry { reason } => {
                write!(f, "degenerate surface geometry: {reason}")
            }
            Self::NonFiniteInput { field } => {
                write!(f, "non-finite input in {field}")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for SurfaceError {}
