// SPDX-License-Identifier: AGPL-3.0-only
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
            PropagationError::InvalidConfiguration(source) => source.into(),
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

    // Display impls
    #[test]
    fn display_ephemeris_unavailable() {
        let e = DynamicsError::EphemerisUnavailable {
            body: "Mars",
            source: None,
        };
        assert!(e.to_string().contains("Mars"));
    }

    #[test]
    fn display_eop_unavailable() {
        let e = DynamicsError::EOPUnavailable { source: None };
        assert!(e.to_string().contains("Earth Orientation"));
    }

    #[test]
    fn display_gravity_coeff_unavailable() {
        let e = DynamicsError::GravityCoefficientUnavailable {
            degree: 3,
            order: 2,
        };
        assert!(e.to_string().contains("C_3,2"));
    }

    #[test]
    fn display_altitude_below_surface() {
        let e = DynamicsError::AltitudeBelowSurface { altitude_km: -5.0 };
        assert!(e.to_string().contains("-5.000"));
    }

    #[test]
    fn display_degenerate_geometry() {
        let e = DynamicsError::DegenerateGeometry {
            reason: "collinear",
        };
        assert!(e.to_string().contains("collinear"));
    }

    #[test]
    fn display_invalid_step_request() {
        let e = DynamicsError::InvalidStepRequest {
            reason: "too large",
        };
        assert!(e.to_string().contains("too large"));
    }

    #[test]
    fn display_atmosphere_provider_error() {
        let e = DynamicsError::AtmosphereProviderError(Box::new(std::io::Error::other("no atm")));
        assert!(e.to_string().contains("atmosphere"));
    }

    #[test]
    fn display_provider_error() {
        let e = DynamicsError::Provider(Box::new(std::io::Error::other("fail")));
        assert!(e.to_string().contains("provider"));
    }

    #[test]
    fn display_gravity_field_unavailable() {
        let e = DynamicsError::GravityFieldUnavailable;
        assert!(e.to_string().contains("gravity field"));
    }

    #[test]
    fn display_geopotential_degree_out_of_range() {
        let e = DynamicsError::GeopotentialDegreeOutOfRange {
            requested: 8,
            max: 4,
        };
        let s = e.to_string();
        assert!(s.contains('8') && s.contains('4'));
    }

    // std::error::Error::source()
    #[test]
    fn error_source_with_source() {
        let inner: Box<dyn std::error::Error + Send + Sync> =
            Box::new(std::io::Error::other("fail"));
        let e = DynamicsError::EOPUnavailable {
            source: Some(inner),
        };
        assert!(std::error::Error::source(&e).is_some());
    }

    #[test]
    fn error_source_without_source() {
        let e = DynamicsError::GravityFieldUnavailable;
        assert!(std::error::Error::source(&e).is_none());
    }

    #[test]
    fn error_source_atmosphere_has_source() {
        let e = DynamicsError::AtmosphereProviderError(Box::new(std::io::Error::other("atm")));
        assert!(std::error::Error::source(&e).is_some());
    }

    // From<PropagationError>
    #[test]
    fn from_propagation_step_below_minimum() {
        let e: DynamicsError = PropagationError::StepBelowMinimum {
            h_requested: 1e-10,
            h_min: 1e-9,
        }
        .into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn from_propagation_max_steps() {
        let e: DynamicsError = PropagationError::MaxStepsExceeded { max_steps: 10_000 }.into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn from_propagation_step_control() {
        let inner = PrincipiaError::StepControlFailed { reason: "diverged" };
        let e: DynamicsError = PropagationError::StepControl(inner).into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn from_propagation_invalid_config() {
        let inner = PrincipiaError::DegenerateGeometry { reason: "test" };
        let e: DynamicsError = PropagationError::InvalidConfiguration(inner).into();
        assert!(matches!(e, DynamicsError::DegenerateGeometry { .. }));
    }

    #[test]
    fn from_propagation_event_evaluation() {
        let src = PrincipiaError::DegenerateGeometry {
            reason: "degenerate",
        };
        let e: DynamicsError = PropagationError::EventEvaluation {
            name: "test_event",
            source: src,
        }
        .into();
        assert!(matches!(e, DynamicsError::DegenerateGeometry { .. }));
    }

    // From<PrincipiaError>
    #[test]
    fn principia_degenerate_geometry_maps() {
        let e: DynamicsError = PrincipiaError::DegenerateGeometry { reason: "zero cp" }.into();
        assert!(matches!(e, DynamicsError::DegenerateGeometry { .. }));
    }

    #[test]
    fn principia_step_control_failed_maps() {
        let e: DynamicsError = PrincipiaError::StepControlFailed {
            reason: "step ctrl",
        }
        .into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn principia_step_below_minimum_maps() {
        let e: DynamicsError = PrincipiaError::StepBelowMinimum {
            reason: "too small",
        }
        .into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn principia_propagation_failed_maps() {
        let e: DynamicsError = PrincipiaError::PropagationFailed { reason: "blowup" }.into();
        assert!(matches!(e, DynamicsError::InvalidStepRequest { .. }));
    }

    #[test]
    fn principia_gravity_coeff_maps() {
        let e: DynamicsError = PrincipiaError::GravityCoefficientUnavailable {
            degree: 5,
            order: 3,
        }
        .into();
        assert!(matches!(
            e,
            DynamicsError::GravityCoefficientUnavailable {
                degree: 5,
                order: 3
            }
        ));
    }

    #[test]
    fn principia_partials_unavailable_maps() {
        let e: DynamicsError = PrincipiaError::PartialsUnavailable { model: "drag" }.into();
        assert!(matches!(e, DynamicsError::Provider(_)));
    }

    #[test]
    fn principia_context_eop_maps() {
        let e: DynamicsError = PrincipiaError::ContextDataUnavailable { what: "eop" }.into();
        assert!(matches!(e, DynamicsError::EOPUnavailable { .. }));
    }

    #[test]
    fn principia_context_gravity_maps() {
        let e: DynamicsError = PrincipiaError::ContextDataUnavailable { what: "gravity" }.into();
        assert!(matches!(e, DynamicsError::GravityFieldUnavailable));
    }

    #[test]
    fn principia_context_gravity_field_maps() {
        let e: DynamicsError = PrincipiaError::ContextDataUnavailable {
            what: "gravity field",
        }
        .into();
        assert!(matches!(e, DynamicsError::GravityFieldUnavailable));
    }

    #[test]
    fn principia_context_atmosphere_maps() {
        let e: DynamicsError = PrincipiaError::ContextDataUnavailable { what: "atmosphere" }.into();
        assert!(matches!(e, DynamicsError::AtmosphereProviderError(_)));
    }

    #[test]
    fn principia_context_unknown_maps() {
        let e: DynamicsError = PrincipiaError::ContextDataUnavailable {
            what: "unknown_field",
        }
        .into();
        assert!(matches!(e, DynamicsError::Provider(_)));
    }

    #[test]
    fn principia_catchall_maps_to_provider() {
        let e: DynamicsError = PrincipiaError::InvalidPropagationConfig { reason: "bad" }.into();
        assert!(matches!(e, DynamicsError::Provider(_)));
    }

    // into_principia()
    #[test]
    fn into_principia_altitude_below_surface() {
        let e = DynamicsError::AltitudeBelowSurface { altitude_km: -1.0 };
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::DegenerateGeometry { .. }
        ));
    }

    #[test]
    fn into_principia_gravity_field_unavailable() {
        let e = DynamicsError::GravityFieldUnavailable;
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::ContextDataUnavailable {
                what: "gravity field"
            }
        ));
    }

    #[test]
    fn into_principia_atmosphere_provider() {
        let e = DynamicsError::AtmosphereProviderError(Box::new(std::io::Error::other("no atm")));
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::ContextDataUnavailable { what: "atmosphere" }
        ));
    }

    #[test]
    fn into_principia_provider() {
        let e = DynamicsError::Provider(Box::new(std::io::Error::other("fail")));
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::ContextDataUnavailable { what: "provider" }
        ));
    }

    #[test]
    fn into_principia_degenerate_geometry() {
        let e = DynamicsError::DegenerateGeometry { reason: "zero vec" };
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::DegenerateGeometry { .. }
        ));
    }

    #[test]
    fn into_principia_invalid_step() {
        let e = DynamicsError::InvalidStepRequest { reason: "overflow" };
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::InvalidStepRequest { .. }
        ));
    }

    #[test]
    fn into_principia_geopotential_out_of_range() {
        let e = DynamicsError::GeopotentialDegreeOutOfRange {
            requested: 8,
            max: 4,
        };
        assert!(matches!(
            e.into_principia(),
            PrincipiaError::GeopotentialDegreeOutOfRange {
                requested: 8,
                max: 4
            }
        ));
    }

    // LocalFrameError Display
    #[test]
    fn local_frame_error_display_parallel() {
        let e = LocalFrameError::PositionAndVelocityParallel;
        assert!(e.to_string().contains("parallel"));
    }

    #[test]
    fn local_frame_error_display_zero_pos() {
        let e = LocalFrameError::ZeroPositionMagnitude;
        assert!(e.to_string().contains("position"));
    }

    #[test]
    fn local_frame_error_display_zero_vel() {
        let e = LocalFrameError::ZeroVelocityMagnitude;
        assert!(e.to_string().contains("velocity"));
    }

    #[test]
    fn local_frame_error_is_std_error() {
        let e = LocalFrameError::ZeroPositionMagnitude;
        assert!(std::error::Error::source(&e).is_none());
    }
}
