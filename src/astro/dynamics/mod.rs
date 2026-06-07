// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomy-specific spacecraft dynamics layer built on top of `principia`.
//!
//! ## Scientific scope
//!
//! Provides force models, runtime context, atmospheric density providers, and
//! a high-level [`Propagator`] facade for geocentric orbit propagation in the
//! GCRS frame on the TT time scale.  Numerical integration kernels, generic
//! integrators, covariance transport, and variational equations live in the
//! `principia` crate; this module adds the astronomy semantics (body
//! constants, ephemeris access, Earth orientation, geodetic altitude) that
//! `principia` deliberately omits.
//!
//! ## Technical scope
//!
//! Entry points for most callers:
//! * [`Propagator`] / [`Dop853Propagator`] / [`Dopri5Propagator`] /
//!   [`Rk4Propagator`] — typed orbit propagators.
//! * [`DynamicsContext`] / [`DynamicsContextBuilder`] — provider injection.
//! * [`forces`] — force-model implementations (two-body, J2, drag, SRP, …).
//! * [`gravity`] — geopotential helpers.
//! * [`density`] — atmospheric density models.
//! * [`state`] — typed orbit state, position, velocity aliases.
//!
//! Types returned by the propagator but defined in `principia`
//! ([`EventDetector`], [`PropagationResult`], [`StateTransitionMatrix`]) are
//! re-exported here so callers of this crate-level API do not need to add a
//! direct `principia` dependency.
//!
//! ## References
//!
//! * Montenbruck & Gill, *Satellite Orbits*, Springer, 2000.
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed.
//! * IERS Conventions 2010, IERS Technical Note 36.

pub mod context;
pub mod covariance;
pub mod density;
pub mod errors;
pub mod forces;
pub mod gravity;
pub mod propagation;
pub mod state;
pub mod units;

pub mod frames {
    //! Astronomy-specific local orbital frame aliases.
    pub use principia::frames::LocalTrajectoryFrame;
    pub use principia::frames::{
        lvlh_from_raw_km_km_s as lvlh_from_state, rtn_from_raw_km_km_s as rtn_from_state,
        vnc_from_raw_km_km_s as vnc_from_state, LVLH, RTN, VNC,
    };

    /// Local orbital frame anchored to the GCRS inertial frame.
    pub type LocalOrbitalFrame<M> = LocalTrajectoryFrame<crate::coordinates::frames::GCRS, M>;
}

// ── Astronomy-specific sub-modules ──────────────────────────────────────────

pub use context::{
    Conventions, DynamicsContext, DynamicsContextBuilder, EarthOrientationProvider,
    PrecessionModel, SolarActivityProvider, TimeScaleHint,
};
pub use density::{
    geodetic_altitude, AtmosphereProvider, ConstantDensity, DensityProvider, ExponentialAtmosphere,
    Nrlmsise00LiteApprox,
};
pub use errors::{DynamicsError, LocalFrameError};
pub use forces::{
    CannonballSrp, CentralBodyRelativity1Pn, Conical, Cylindrical, DragForce, EclipseModel,
    EmpiricalAcceleration, ExponentialDrag, Geopotential, NoEclipse, ShadowModel, SunPerturbation,
    ThirdBody, ThirdBodyProvider, TwoBody, AU_IN_KM, DEGENERATE_RADIUS_KM, EARTH_J2, GM_EARTH,
    GM_MOON, GM_SUN, J2, OMEGA_EARTH_RAD_S, P0, R_EARTH,
};
pub use gravity::{
    spherical_harmonic_acceleration, GravityConstants, GravityFieldProvider, LowDegreeEarth,
    TwoBodyEarth,
};
pub use propagation::{
    Dop853Propagator, Dopri5Propagator, EventDetector, EventOccurrence, FixedRk4Adapter,
    PropagationError, PropagationResult, Propagator, PropagatorConfig, RadialThresholdEvent,
    Rk4Propagator,
};
pub use state::{
    Acceleration, AccelerationUnit, Force, OrbitState, Position, SpacecraftProperties,
    SpacecraftState, Velocity, VelocityUnit,
};
pub use units::GravitationalParameter;

// ── `principia` types used in the Propagator return API ─────────────────────
//
// Re-exported so callers of `siderust::astro::dynamics` do not need a direct
// `principia` dependency to work with propagation results.
pub use covariance::StateCovariance;
pub use principia::StateTransitionMatrix;
