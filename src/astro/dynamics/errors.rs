// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error types for the astronomy-specific dynamics layer.

use std::fmt;

use principia::{PrincipiaError, PropagationError};

/// Errors produced by force models, runtime-context accessors, and astronomy-specific wrappers.
#[derive(Debug)]
pub enum DynamicsError {
    /// An ephemeris provider could not return the requested body's position.
    EphemerisUnavailable {
        /// Human-readable body name (e.g. `"Sun"`, `"Moon"`).
        body: &'static str,
        /// Underlying provider error, if one was returned.
        source: Option<Box<dyn std::error::Error + Send + Sync>>,
    },
    /// Earth Orientation Parameters (EOP) are required but not available.
    EOPUnavailable {
        /// Underlying provider error, if one was returned.
        source: Option<Box<dyn std::error::Error + Send + Sync>>,
    },
    /// A geopotential coefficient at degree/order `(n, m)` is not available.
    GravityCoefficientUnavailable {
        /// Spherical-harmonic degree `n`.
        degree: u16,
        /// Spherical-harmonic order `m`.
        order: u16,
    },
    /// The spacecraft is below the planetary surface (altitude < 0).
    AltitudeBelowSurface {
        /// Computed altitude in kilometres (may be negative).
        altitude_km: f64,
    },
    /// A geometric computation degenerated.
    DegenerateGeometry {
        /// Short human-readable explanation of what went wrong.
        reason: &'static str,
    },
    /// An integration step size or count is invalid.
    InvalidStepRequest {
        /// Short human-readable explanation of what constraint was violated.
        reason: &'static str,
    },
    /// An atmosphere density provider returned an error.
    AtmosphereProviderError(Box<dyn std::error::Error + Send + Sync>),
    /// An opaque provider error not covered by the more specific variants.
    Provider(Box<dyn std::error::Error + Send + Sync>),
    /// A gravity field is required but no provider was set in the context.
    GravityFieldUnavailable,
    /// The requested degree/order exceeds what the gravity field provider supports.
    GeopotentialDegreeOutOfRange {
        /// Degree requested by the force model.
        requested: usize,
        /// Maximum degree the provider supports.
        max: usize,
    },
}

impl DynamicsError {
    /// Collapse this astronomy-layer error into the generic `principia` error surface.
    pub(crate) fn into_principia(self) -> PrincipiaError {
        match self {
            Self::EphemerisUnavailable { .. } => {
                PrincipiaError::ContextDataUnavailable { what: "ephemeris" }
            }
            Self::EOPUnavailable { .. } => PrincipiaError::ContextDataUnavailable { what: "eop" },
            Self::GravityCoefficientUnavailable { degree, order } => {
                PrincipiaError::GravityCoefficientUnavailable { degree, order }
            }
            Self::AltitudeBelowSurface { .. } => PrincipiaError::DegenerateGeometry {
                reason: "altitude below surface",
            },
            Self::DegenerateGeometry { reason } => PrincipiaError::DegenerateGeometry { reason },
            Self::InvalidStepRequest { reason } => PrincipiaError::InvalidStepRequest { reason },
            Self::AtmosphereProviderError(_) => {
                PrincipiaError::ContextDataUnavailable { what: "atmosphere" }
            }
            Self::Provider(_) => PrincipiaError::ContextDataUnavailable { what: "provider" },
            Self::GravityFieldUnavailable => PrincipiaError::ContextDataUnavailable {
                what: "gravity field",
            },
            Self::GeopotentialDegreeOutOfRange { requested, max } => {
                PrincipiaError::GeopotentialDegreeOutOfRange { requested, max }
            }
        }
    }
}

impl From<PrincipiaError> for DynamicsError {
    fn from(value: PrincipiaError) -> Self {
        match value {
            PrincipiaError::DegenerateGeometry { reason } => Self::DegenerateGeometry { reason },
            PrincipiaError::InvalidStepRequest { reason }
            | PrincipiaError::StepControlFailed { reason }
            | PrincipiaError::StepBelowMinimum { reason }
            | PrincipiaError::PropagationFailed { reason } => Self::InvalidStepRequest { reason },
            PrincipiaError::GravityCoefficientUnavailable { degree, order } => {
                Self::GravityCoefficientUnavailable { degree, order }
            }
            PrincipiaError::GeopotentialDegreeOutOfRange { requested, max } => {
                Self::GeopotentialDegreeOutOfRange { requested, max }
            }
            PrincipiaError::PartialsUnavailable { model } => {
                Self::Provider(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    format!("analytic partials not available for model '{model}'"),
                )))
            }
            PrincipiaError::ContextDataUnavailable { what } => match what {
                "ephemeris" => Self::EphemerisUnavailable {
                    body: "(any)",
                    source: None,
                },
                "eop" => Self::EOPUnavailable { source: None },
                "gravity" | "gravity field" => Self::GravityFieldUnavailable,
                "atmosphere" => Self::AtmosphereProviderError(Box::new(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "no atmosphere provider in DynamicsContext",
                ))),
                other => Self::Provider(Box::new(std::io::Error::other(format!(
                    "context data unavailable: {other}"
                )))),
            },
            _ => Self::Provider(Box::new(std::io::Error::other("unmapped principia error"))),
        }
    }
}

impl From<PropagationError> for DynamicsError {
    fn from(value: PropagationError) -> Self {
        match value {
            PropagationError::StepControl(source) => source.into(),
            PropagationError::StepBelowMinimum { .. } => Self::InvalidStepRequest {
                reason: "requested step below configured minimum",
            },
            PropagationError::MaxStepsExceeded { .. } => Self::InvalidStepRequest {
                reason: "propagator max_steps exceeded",
            },
            PropagationError::EventEvaluation { source, .. } => source.into(),
        }
    }
}

impl fmt::Display for DynamicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EphemerisUnavailable { body, .. } => {
                write!(f, "ephemeris unavailable for body '{body}'")
            }
            Self::EOPUnavailable { .. } => write!(f, "Earth Orientation Parameters unavailable"),
            Self::GravityCoefficientUnavailable { degree, order } => {
                write!(
                    f,
                    "gravity coefficient C_{degree},{order} not available in current model"
                )
            }
            Self::AltitudeBelowSurface { altitude_km } => {
                write!(f, "altitude {altitude_km:.3} km is below the surface")
            }
            Self::DegenerateGeometry { reason } => write!(f, "degenerate geometry: {reason}"),
            Self::InvalidStepRequest { reason } => write!(f, "invalid step request: {reason}"),
            Self::AtmosphereProviderError(e) => write!(f, "atmosphere provider error: {e}"),
            Self::Provider(e) => write!(f, "provider error: {e}"),
            Self::GravityFieldUnavailable => {
                write!(f, "no gravity field provider set in dynamics context")
            }
            Self::GeopotentialDegreeOutOfRange { requested, max } => write!(
                f,
                "requested geopotential degree {requested} exceeds provider maximum {max}"
            ),
        }
    }
}

impl std::error::Error for DynamicsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::EphemerisUnavailable { source, .. } => source
                .as_ref()
                .map(|e| e.as_ref() as &(dyn std::error::Error + 'static)),
            Self::EOPUnavailable { source } => source
                .as_ref()
                .map(|e| e.as_ref() as &(dyn std::error::Error + 'static)),
            Self::AtmosphereProviderError(e) | Self::Provider(e) => Some(e.as_ref()),
            _ => None,
        }
    }
}

/// Errors produced when constructing a local orbital frame from a degenerate state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalFrameError {
    /// The position and velocity vectors are parallel (or anti-parallel).
    PositionAndVelocityParallel,
    /// The position vector has zero magnitude.
    ZeroPositionMagnitude,
    /// The velocity vector has zero magnitude.
    ZeroVelocityMagnitude,
}

impl fmt::Display for LocalFrameError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PositionAndVelocityParallel => write!(
                f,
                "cannot build local orbital frame: position and velocity are parallel"
            ),
            Self::ZeroPositionMagnitude => write!(
                f,
                "cannot build local orbital frame: position vector has zero magnitude"
            ),
            Self::ZeroVelocityMagnitude => write!(
                f,
                "cannot build local orbital frame: velocity vector has zero magnitude"
            ),
        }
    }
}

impl std::error::Error for LocalFrameError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn principia_geopotential_error_maps_back() {
        let err: DynamicsError = PrincipiaError::GeopotentialDegreeOutOfRange {
            requested: 10,
            max: 4,
        }
        .into();
        assert!(matches!(
            err,
            DynamicsError::GeopotentialDegreeOutOfRange {
                requested: 10,
                max: 4
            }
        ));
    }

    #[test]
    fn astronomy_ephemeris_error_maps_to_context_unavailable() {
        let err = DynamicsError::EphemerisUnavailable {
            body: "Sun",
            source: None,
        };
        assert!(matches!(
            err.into_principia(),
            PrincipiaError::ContextDataUnavailable { what: "ephemeris" }
        ));
    }
}
