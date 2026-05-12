// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error types for the dynamics module.
//!
//! [`DynamicsError`] covers every failure mode that a force model, integrator,
//! or provider accessor can encounter at runtime.  [`LocalFrameError`] covers
//! the degenerate-geometry cases that arise when constructing an orbital local
//! frame from a degenerate state.
//!
//! ## Design notes
//!
//! Both types implement the standard [`std::error::Error`] trait and are
//! `Send + Sync` so they can cross thread and async-task boundaries without
//! wrapping.  Variants that originate from a lower-level provider carry an
//! optional `source` chain so callers can inspect the root cause.

use std::fmt;

// =============================================================================
// DynamicsError
// =============================================================================

/// Errors produced by force models, integrators, and dynamics providers.
///
/// Variants are non-exhaustive in spirit — new ones may be added as the
/// dynamics kernel grows.  Match on the variants you care about and use a
/// wildcard arm for forward compatibility.
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

    /// A geopotential coefficient at degree/order `(n, m)` is not available
    /// from the current [`GravityFieldProvider`](super::gravity::GravityFieldProvider).
    GravityCoefficientUnavailable {
        /// Spherical-harmonic degree `n`.
        degree: u16,
        /// Spherical-harmonic order `m`.
        order: u16,
    },

    /// The spacecraft is below the planetary surface (altitude < 0).
    ///
    /// Force models that assume a valid atmospheric or geopotential altitude
    /// return this variant rather than extrapolating into nonsensical regimes.
    AltitudeBelowSurface {
        /// Computed altitude in kilometres (may be negative).
        altitude_km: f64,
    },

    /// A geometric computation degenerated (e.g. zero cross-product, singular
    /// rotation matrix).
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

    /// The requested degree/order exceeds what the gravity field provider
    /// supports.
    GeopotentialDegreeOutOfRange {
        /// Degree requested by the force model.
        requested: usize,
        /// Maximum degree the provider supports.
        max: usize,
    },
}

impl fmt::Display for DynamicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EphemerisUnavailable { body, .. } => {
                write!(f, "ephemeris unavailable for body '{body}'")
            }
            Self::EOPUnavailable { .. } => {
                write!(f, "Earth Orientation Parameters unavailable")
            }
            Self::GravityCoefficientUnavailable { degree, order } => {
                write!(
                    f,
                    "gravity coefficient C_{degree},{order} not available in current model"
                )
            }
            Self::AltitudeBelowSurface { altitude_km } => {
                write!(f, "altitude {altitude_km:.3} km is below the surface")
            }
            Self::DegenerateGeometry { reason } => {
                write!(f, "degenerate geometry: {reason}")
            }
            Self::InvalidStepRequest { reason } => {
                write!(f, "invalid step request: {reason}")
            }
            Self::AtmosphereProviderError(e) => {
                write!(f, "atmosphere provider error: {e}")
            }
            Self::Provider(e) => {
                write!(f, "provider error: {e}")
            }
            Self::GravityFieldUnavailable => {
                write!(f, "no gravity field provider set in dynamics context")
            }
            Self::GeopotentialDegreeOutOfRange { requested, max } => {
                write!(
                    f,
                    "requested geopotential degree {requested} exceeds provider maximum {max}"
                )
            }
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
            Self::AtmosphereProviderError(e) => Some(e.as_ref()),
            Self::Provider(e) => Some(e.as_ref()),
            _ => None,
        }
    }
}

// =============================================================================
// LocalFrameError
// =============================================================================

/// Errors produced when constructing a local orbital frame from a degenerate
/// state.
///
/// Returned by the fallible `from_state_checked` constructors on
/// [`LocalOrbitalFrame`](super::frames::LocalOrbitalFrame) once they are
/// converted to return `Result` (tracked in a later todo).  Also used by any
/// code that detects degenerate geometry before delegating to the frame type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalFrameError {
    /// The position and velocity vectors are parallel (or anti-parallel), so
    /// the orbit-normal direction `r × v` is zero and the frame is undefined.
    PositionAndVelocityParallel,

    /// The position vector has zero magnitude; the radial direction is
    /// undefined.
    ZeroPositionMagnitude,

    /// The velocity vector has zero magnitude; the along-track / VNC
    /// directions are undefined.
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

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    // ---- DynamicsError Display ----

    #[test]
    fn ephemeris_unavailable_display() {
        let e = DynamicsError::EphemerisUnavailable {
            body: "Sun",
            source: None,
        };
        assert_eq!(e.to_string(), "ephemeris unavailable for body 'Sun'");
    }

    #[test]
    fn eop_unavailable_display() {
        let e = DynamicsError::EOPUnavailable { source: None };
        assert_eq!(e.to_string(), "Earth Orientation Parameters unavailable");
    }

    #[test]
    fn gravity_coefficient_display() {
        let e = DynamicsError::GravityCoefficientUnavailable {
            degree: 8,
            order: 3,
        };
        assert_eq!(
            e.to_string(),
            "gravity coefficient C_8,3 not available in current model"
        );
    }

    #[test]
    fn altitude_below_surface_display() {
        let e = DynamicsError::AltitudeBelowSurface { altitude_km: -12.5 };
        assert!(e.to_string().contains("-12.500"));
    }

    #[test]
    fn degenerate_geometry_display() {
        let e = DynamicsError::DegenerateGeometry {
            reason: "zero cross-product",
        };
        assert!(e.to_string().contains("zero cross-product"));
    }

    #[test]
    fn invalid_step_request_display() {
        let e = DynamicsError::InvalidStepRequest {
            reason: "step size must be positive",
        };
        assert!(e.to_string().contains("step size must be positive"));
    }

    // ---- DynamicsError source() forwarding ----

    #[test]
    fn ephemeris_unavailable_source_none() {
        let e = DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: None,
        };
        assert!(e.source().is_none());
    }

    #[test]
    fn ephemeris_unavailable_source_some() {
        let inner: Box<dyn Error + Send + Sync> =
            Box::new(std::io::Error::new(std::io::ErrorKind::NotFound, "no data"));
        let e = DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: Some(inner),
        };
        assert!(e.source().is_some());
    }

    #[test]
    fn eop_unavailable_source_forwarded() {
        let inner: Box<dyn Error + Send + Sync> =
            Box::new(std::io::Error::new(std::io::ErrorKind::TimedOut, "timeout"));
        let e = DynamicsError::EOPUnavailable {
            source: Some(inner),
        };
        assert!(e.source().is_some());
    }

    #[test]
    fn provider_error_source_forwarded() {
        let inner: Box<dyn Error + Send + Sync> =
            Box::new(std::io::Error::new(std::io::ErrorKind::Other, "boom"));
        let e = DynamicsError::Provider(inner);
        assert!(e.source().is_some());
        assert!(e.to_string().contains("boom"));
    }

    #[test]
    fn atmosphere_provider_error_source_forwarded() {
        let inner: Box<dyn Error + Send + Sync> =
            Box::new(std::io::Error::new(std::io::ErrorKind::Other, "atm error"));
        let e = DynamicsError::AtmosphereProviderError(inner);
        assert!(e.source().is_some());
    }

    // ---- LocalFrameError Display ----

    #[test]
    fn local_frame_error_parallel_display() {
        let e = LocalFrameError::PositionAndVelocityParallel;
        assert!(e.to_string().contains("parallel"));
    }

    #[test]
    fn local_frame_error_zero_position_display() {
        let e = LocalFrameError::ZeroPositionMagnitude;
        assert!(e.to_string().contains("position vector has zero magnitude"));
    }

    #[test]
    fn local_frame_error_zero_velocity_display() {
        let e = LocalFrameError::ZeroVelocityMagnitude;
        assert!(e.to_string().contains("velocity vector has zero magnitude"));
    }

    #[test]
    fn local_frame_error_is_error_trait() {
        let e: &dyn Error = &LocalFrameError::PositionAndVelocityParallel;
        assert!(e.source().is_none());
    }
}
