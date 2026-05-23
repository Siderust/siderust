// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI for spacecraft dynamics: contexts, orbit states, and propagators.
//!
//! ## Lifecycle
//!
//! ### Dynamics context
//!
//! 1. Create an empty context: `siderust_dynamics_context_new(out)`.
//! 2. Optionally attach providers:
//!    - `siderust_dynamics_context_with_ephemeris(ctx, eph_handle)` — shares a
//!      previously loaded `SiderustRuntimeEphemeris`.
//!    - `siderust_dynamics_context_with_atmosphere(ctx, vtable)` — injects a
//!      C-callback density provider.
//!    - `siderust_dynamics_context_with_gravity_field(ctx, vtable)` — injects a
//!      C-callback gravity field.
//! 3. Free when done: `siderust_dynamics_context_free(ctx)`.
//!
//! ### Orbit state
//!
//! - `siderust_orbit_state_new(jd, x, y, z, vx, vy, vz, out)` — construct from
//!   raw doubles (position km, velocity km/s, epoch as Julian Date).
//! - `siderust_orbit_state_position(handle, out_x, out_y, out_z)`
//! - `siderust_orbit_state_velocity(handle, out_vx, out_vy, out_vz)`
//! - `siderust_orbit_state_epoch_jd(handle, out_jd)`
//! - `siderust_orbit_state_free(handle)`
//!
//! ### Two-body propagator
//!
//! 1. Create: `siderust_propagator_two_body_earth_new(out)` (DOP853, Earth GM).
//!    Or `siderust_propagator_two_body_new(gm_km3_s2, out)` for a custom GM.
//! 2. Propagate: `siderust_propagator_propagate(handle, state, dt_s, ctx, out)`.
//!    Returns a [`SiderustDynamicsStatus`]; on `Ok` writes the final state
//!    through `out`.
//! 3. Free: `siderust_propagator_free(handle)`.
//!
//! ## Status codes
//!
//! Functions that can produce dynamics-specific failures return
//! [`SiderustDynamicsStatus`] instead of the generic [`SiderustStatus`].
//! The `Ok` variant (0) indicates success; all other values map to a specific
//! [`DynamicsError`] variant.

use std::os::raw::c_void;
use std::sync::Arc;

use principia::integrators::dop853_propagate;
use principia::PrincipiaError;
use qtty::tolerances::IntegratorTolerances;
use siderust::astro::dynamics::context::DynamicsContext;
use siderust::astro::dynamics::density::DensityProvider;
use siderust::astro::dynamics::errors::DynamicsError;
use siderust::astro::dynamics::forces::TwoBody;
use siderust::astro::dynamics::gravity::GravityFieldProvider;
use siderust::astro::dynamics::units::{GravitationalParameter, GM_EARTH};
use siderust::astro::dynamics::{OrbitState, Position, Velocity};
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::frames::GCRS;
use siderust::ephemeris::DynEphemeris;
use siderust::qtty::{KilogramsPerCubicMeter, Kilometers, Second};
use siderust::time::JulianDate;

use crate::error::SiderustStatus;
use crate::runtime_ephemeris::SiderustRuntimeEphemeris;

// =============================================================================
// Panic-catching guard for dynamics functions
// =============================================================================

macro_rules! dyn_guard {
    ($body:block) => {{
        let result = ::std::panic::catch_unwind(::std::panic::AssertUnwindSafe(|| $body));
        match result {
            Ok(s) => s,
            Err(_) => SiderustDynamicsStatus::InternalPanic,
        }
    }};
}

// =============================================================================
// SiderustDynamicsStatus — mirrors DynamicsError
// =============================================================================

/// Status codes returned by dynamics propagation functions.
///
/// Maps to the variants of [`siderust::astro::dynamics::errors::DynamicsError`].
///
/// | Value | Meaning |
/// |-------|---------|
/// | 0 | Success |
/// | 1 | Required pointer was null |
/// | 2 | No ephemeris provider configured |
/// | 3 | EOP data unavailable |
/// | 4 | A gravity spherical-harmonic coefficient is unavailable |
/// | 5 | Altitude is below the body surface |
/// | 6 | Degenerate geometry (zero position or velocity, r ∥ v) |
/// | 7 | Invalid integrator step request |
/// | 8 | Atmosphere density provider returned an error |
/// | 9 | An opaque provider error occurred |
/// | 10 | No gravity field provider configured in the context |
/// | 11 | Requested geopotential degree/order exceeds the provider limit |
/// | 99 | A Rust panic was caught at the FFI boundary |
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustDynamicsStatus {
    /// Propagation succeeded.
    Ok = 0,
    /// A required pointer argument was null.
    NullPointer = 1,
    /// No ephemeris provider configured in the dynamics context.
    EphemerisUnavailable = 2,
    /// Earth Orientation Parameters unavailable for the requested epoch.
    EopUnavailable = 3,
    /// A gravity spherical-harmonic coefficient is unavailable.
    GravityCoefficientUnavailable = 4,
    /// Spacecraft altitude is below the body surface (negative altitude).
    AltitudeBelowSurface = 5,
    /// Degenerate geometry (zero position, zero velocity, or r ∥ v).
    DegenerateGeometry = 6,
    /// Integrator step request violated a constraint (e.g. step size is zero).
    InvalidStepRequest = 7,
    /// The atmosphere density provider returned an error.
    AtmosphereProviderError = 8,
    /// An opaque provider error not covered by the more specific variants.
    ProviderError = 9,
    /// No gravity field provider configured in the dynamics context.
    GravityFieldUnavailable = 10,
    /// Requested geopotential degree or order exceeds the provider's maximum.
    GeopotentialDegreeOutOfRange = 11,
    /// A Rust panic was caught at the FFI boundary.
    InternalPanic = 99,
}

impl SiderustDynamicsStatus {
    fn from_dynamics_error(e: &DynamicsError) -> Self {
        match e {
            DynamicsError::EphemerisUnavailable { .. } => Self::EphemerisUnavailable,
            DynamicsError::EOPUnavailable { .. } => Self::EopUnavailable,
            DynamicsError::GravityCoefficientUnavailable { .. } => {
                Self::GravityCoefficientUnavailable
            }
            DynamicsError::AltitudeBelowSurface { .. } => Self::AltitudeBelowSurface,
            DynamicsError::DegenerateGeometry { .. } => Self::DegenerateGeometry,
            DynamicsError::InvalidStepRequest { .. } => Self::InvalidStepRequest,
            DynamicsError::AtmosphereProviderError(_) => Self::AtmosphereProviderError,
            DynamicsError::Provider(_) => Self::ProviderError,
            DynamicsError::GravityFieldUnavailable => Self::GravityFieldUnavailable,
            DynamicsError::GeopotentialDegreeOutOfRange { .. } => {
                Self::GeopotentialDegreeOutOfRange
            }
        }
    }
}

// =============================================================================
// C vtable types for pluggable providers
// =============================================================================

/// C-callback vtable for an atmosphere density provider.
///
/// Pass a pointer to this struct to
/// `siderust_dynamics_context_with_atmosphere`. The function pointer
/// `density` will be called during propagation with the spacecraft's geodetic
/// altitude (km); it must return the air density in kg/m³, or a **negative**
/// value to signal an error (e.g. below-surface altitude).
///
/// `user_data` is passed unchanged to every call; use it to carry a pointer to
/// your density-model state.  You are responsible for its lifetime: keep the
/// data alive until the context is freed.
#[repr(C)]
pub struct SiderustAtmosphereVtable {
    /// Compute atmospheric density (kg/m³) at `altitude_km` above the surface.
    ///
    /// Return a negative value to signal that the altitude is below the surface
    /// (the propagator will map this to `SIDERUST_DYNAMICS_STATUS_ALTITUDE_BELOW_SURFACE`).
    pub density: unsafe extern "C" fn(altitude_km: f64, user_data: *mut c_void) -> f64,
    /// Caller-controlled context forwarded to every `density` call.
    pub user_data: *mut c_void,
}

// SAFETY: The C caller is responsible for ensuring the user_data is thread-safe.
unsafe impl Send for SiderustAtmosphereVtable {}
unsafe impl Sync for SiderustAtmosphereVtable {}

/// C-callback vtable for a spherical-harmonic gravity field provider.
///
/// Pass a pointer to this struct to
/// `siderust_dynamics_context_with_gravity_field`.  The function pointers
/// `c_normalized` and `s_normalized` are called during propagation to retrieve
/// fully-normalised Stokes coefficients C̄_{n,m} and S̄_{n,m}.
///
/// `user_data` is forwarded to every coefficient call.
#[repr(C)]
pub struct SiderustGravityVtable {
    /// Gravitational parameter GM (km³/s²).
    pub gm_km3_s2: f64,
    /// Equatorial reference radius (km).
    pub reference_radius_km: f64,
    /// Maximum spherical-harmonic degree available.
    pub max_degree: u32,
    /// Maximum spherical-harmonic order available.
    pub max_order: u32,
    /// Return the fully-normalised cosine coefficient C̄_{n,m}.
    pub c_normalized: unsafe extern "C" fn(n: u32, m: u32, user_data: *mut c_void) -> f64,
    /// Return the fully-normalised sine coefficient S̄_{n,m}.
    pub s_normalized: unsafe extern "C" fn(n: u32, m: u32, user_data: *mut c_void) -> f64,
    /// Caller-controlled context forwarded to every coefficient call.
    pub user_data: *mut c_void,
}

// SAFETY: The C caller is responsible for ensuring the user_data is thread-safe.
unsafe impl Send for SiderustGravityVtable {}
unsafe impl Sync for SiderustGravityVtable {}

// =============================================================================
// Internal provider wrappers
// =============================================================================

struct FfiDensityProvider {
    vtable: SiderustAtmosphereVtable,
}

impl DensityProvider for FfiDensityProvider {
    fn name(&self) -> &'static str {
        "ffi_atmosphere"
    }

    fn density(&self, altitude: Kilometers) -> Result<KilogramsPerCubicMeter, DynamicsError> {
        // SAFETY: The C caller guarantees the fn pointer and user_data are valid.
        let rho = unsafe { (self.vtable.density)(altitude.value(), self.vtable.user_data) };
        if rho < 0.0 {
            Err(DynamicsError::AltitudeBelowSurface {
                altitude_km: altitude.value(),
            })
        } else {
            Ok(KilogramsPerCubicMeter::new(rho))
        }
    }
}

// SAFETY: vtable is Send+Sync (see above).
unsafe impl Send for FfiDensityProvider {}
unsafe impl Sync for FfiDensityProvider {}

struct FfiGravityProvider {
    vtable: SiderustGravityVtable,
}

impl GravityFieldProvider for FfiGravityProvider {
    fn mu(&self) -> GravitationalParameter {
        GravitationalParameter::new(self.vtable.gm_km3_s2)
    }

    fn reference_radius(&self) -> Kilometers {
        Kilometers::new(self.vtable.reference_radius_km)
    }

    fn max_degree(&self) -> usize {
        self.vtable.max_degree as usize
    }

    fn max_order(&self) -> usize {
        self.vtable.max_order as usize
    }

    fn c_normalized(&self, n: usize, m: usize) -> Result<f64, PrincipiaError> {
        // SAFETY: The C caller guarantees the fn pointer and user_data are valid.
        Ok(unsafe { (self.vtable.c_normalized)(n as u32, m as u32, self.vtable.user_data) })
    }

    fn s_normalized(&self, n: usize, m: usize) -> Result<f64, PrincipiaError> {
        // SAFETY: The C caller guarantees the fn pointer and user_data are valid.
        Ok(unsafe { (self.vtable.s_normalized)(n as u32, m as u32, self.vtable.user_data) })
    }
}

// SAFETY: vtable is Send+Sync (see above).
unsafe impl Send for FfiGravityProvider {}
unsafe impl Sync for FfiGravityProvider {}

// Helper: arc a RuntimeEphemeris clone so it can be shared into DynamicsContext.

// =============================================================================
// SiderustDynamicsContext — opaque handle
// =============================================================================

/// Opaque handle to a [`DynamicsContext`] used by propagation entry points.
///
/// Created via `siderust_dynamics_context_new` and freed via
/// `siderust_dynamics_context_free`.  Providers are attached via the
/// `siderust_dynamics_context_with_*` family.
pub struct SiderustDynamicsContext {
    pub(crate) inner: DynamicsContext,
}

/// Create a dynamics context with no providers (suitable for pure two-body
/// propagation).
///
/// On success, `*out` receives a pointer to a newly heap-allocated context.
/// Free the context with `siderust_dynamics_context_free` when done.
#[no_mangle]
pub extern "C" fn siderust_dynamics_context_new(
    out: *mut *mut SiderustDynamicsContext,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let ctx = Box::new(SiderustDynamicsContext {
            inner: DynamicsContext::empty(),
        });
        unsafe { *out = Box::into_raw(ctx) };
        SiderustStatus::Ok
    }}
}

/// Free a dynamics context previously created by `siderust_dynamics_context_new`.
///
/// # Safety
///
/// `handle` must be either null or a live pointer produced by this crate.
#[no_mangle]
pub unsafe extern "C" fn siderust_dynamics_context_free(handle: *mut SiderustDynamicsContext) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Attach a runtime ephemeris provider to the dynamics context.
///
/// Clones the underlying `RuntimeEphemeris` so both the context and the
/// original handle remain independently valid.  `eph` must be a live handle
/// previously returned by `siderust_runtime_ephemeris_load_bsp` or
/// `siderust_runtime_ephemeris_load_bytes`.
#[no_mangle]
pub extern "C" fn siderust_dynamics_context_with_ephemeris(
    ctx: *mut SiderustDynamicsContext,
    eph: *const SiderustRuntimeEphemeris,
) -> SiderustStatus {
    ffi_guard! {{
        if ctx.is_null() || eph.is_null() {
            return SiderustStatus::NullPointer;
        }
        // Clone the inner RuntimeEphemeris and box it as an Arc<dyn DynEphemeris>.
        let eph_arc: Arc<dyn DynEphemeris + Send + Sync> =
            Arc::new(unsafe { (*eph).inner.clone() });
        unsafe { (*ctx).inner.ephemeris = Some(eph_arc) };
        SiderustStatus::Ok
    }}
}

/// Attach an atmosphere density provider (C-callback vtable) to the context.
///
/// The `vtable` is copied by value; the caller's `user_data` pointer must
/// remain valid until the context is freed.
///
/// # Safety
///
/// The `density` function pointer in `vtable` must remain valid for the
/// lifetime of the context.  `user_data` must be safe to pass across thread
/// boundaries (the propagator may call it from a different thread).
#[no_mangle]
pub extern "C" fn siderust_dynamics_context_with_atmosphere(
    ctx: *mut SiderustDynamicsContext,
    vtable: *const SiderustAtmosphereVtable,
) -> SiderustStatus {
    ffi_guard! {{
        if ctx.is_null() || vtable.is_null() {
            return SiderustStatus::NullPointer;
        }
        // Copy the vtable so the caller does not need to keep it alive.
        let vt = unsafe { std::ptr::read(vtable) };
        let provider: Arc<dyn DensityProvider + Send + Sync> =
            Arc::new(FfiDensityProvider { vtable: vt });
        unsafe { (*ctx).inner.atmosphere = Some(provider) };
        SiderustStatus::Ok
    }}
}

/// Attach a spherical-harmonic gravity field provider (C-callback vtable) to
/// the context.
///
/// The `vtable` is copied by value.  The function pointers and `user_data`
/// must remain valid until the context is freed.
///
/// # Safety
///
/// The `c_normalized` / `s_normalized` function pointers must remain valid for
/// the lifetime of the context.  `user_data` must be thread-safe.
#[no_mangle]
pub extern "C" fn siderust_dynamics_context_with_gravity_field(
    ctx: *mut SiderustDynamicsContext,
    vtable: *const SiderustGravityVtable,
) -> SiderustStatus {
    ffi_guard! {{
        if ctx.is_null() || vtable.is_null() {
            return SiderustStatus::NullPointer;
        }
        let vt = unsafe { std::ptr::read(vtable) };
        let provider: Arc<dyn GravityFieldProvider + Send + Sync> =
            Arc::new(FfiGravityProvider { vtable: vt });
        unsafe { (*ctx).inner.gravity = Some(provider) };
        SiderustStatus::Ok
    }}
}

// =============================================================================
// SiderustOrbitState — opaque handle
// =============================================================================

/// Opaque handle to an [`OrbitState`] in the Geocentric/GCRS inertial frame.
///
/// Units:
/// - Position: km
/// - Velocity: km/s
/// - Epoch: Julian Date (TT)
pub struct SiderustOrbitState {
    pub(crate) inner: OrbitState<Geocentric, GCRS>,
}

/// Construct an orbit state from raw components.
///
/// - `epoch_jd` — epoch as a Julian Date (TT).
/// - `x`, `y`, `z` — position components in km (Geocentric/GCRS).
/// - `vx`, `vy`, `vz` — velocity components in km/s (GCRS).
///
/// On success, `*out` receives a pointer to a new handle.
/// Free with `siderust_orbit_state_free`.
#[no_mangle]
pub extern "C" fn siderust_orbit_state_new(
    epoch_jd: f64,
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    out: *mut *mut SiderustOrbitState,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let state = OrbitState::new(
            JulianDate::new(epoch_jd).to_j2000s(),
            Position::<GCRS>::new(x, y, z),
            Velocity::<GCRS>::new(vx, vy, vz),
        );
        let handle = Box::new(SiderustOrbitState { inner: state });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Free an orbit state handle previously returned by this crate.
///
/// # Safety
///
/// `handle` must be either null or a live pointer produced by this crate.
#[no_mangle]
pub unsafe extern "C" fn siderust_orbit_state_free(handle: *mut SiderustOrbitState) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Read the position components (km, Geocentric/GCRS) from an orbit state.
#[no_mangle]
pub extern "C" fn siderust_orbit_state_position(
    handle: *const SiderustOrbitState,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out_x.is_null() || out_y.is_null() || out_z.is_null() {
            return SiderustStatus::NullPointer;
        }
        let s = unsafe { &(*handle).inner };
        unsafe {
            *out_x = s.position.x().value();
            *out_y = s.position.y().value();
            *out_z = s.position.z().value();
        }
        SiderustStatus::Ok
    }}
}

/// Read the velocity components (km/s, GCRS) from an orbit state.
#[no_mangle]
pub extern "C" fn siderust_orbit_state_velocity(
    handle: *const SiderustOrbitState,
    out_vx: *mut f64,
    out_vy: *mut f64,
    out_vz: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out_vx.is_null() || out_vy.is_null() || out_vz.is_null() {
            return SiderustStatus::NullPointer;
        }
        let s = unsafe { &(*handle).inner };
        unsafe {
            *out_vx = s.velocity.x().value();
            *out_vy = s.velocity.y().value();
            *out_vz = s.velocity.z().value();
        }
        SiderustStatus::Ok
    }}
}

/// Read the epoch (Julian Date, TT) from an orbit state.
#[no_mangle]
pub extern "C" fn siderust_orbit_state_epoch_jd(
    handle: *const SiderustOrbitState,
    out_jd: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out_jd.is_null() {
            return SiderustStatus::NullPointer;
        }
        let epoch_jd = unsafe { (*handle).inner.epoch.to::<siderust::JD>().raw().value() };
        unsafe { *out_jd = epoch_jd };
        SiderustStatus::Ok
    }}
}

// =============================================================================
// SiderustPropagator — opaque handle (TwoBody, DOP853)
// =============================================================================

/// Opaque handle to a DOP853 two-body propagator.
///
/// Stores the [`TwoBody`] force model; propagation is performed directly via
/// the DOP853 adaptive integrator.
/// Created via `siderust_propagator_two_body_earth_new` or
/// `siderust_propagator_two_body_new`.
pub struct SiderustPropagator {
    pub(crate) force: TwoBody,
}

/// Create a DOP853 two-body propagator using Earth's gravitational parameter
/// (EGM2008 GM = 398 600.441 8 km³/s²).
///
/// On success, `*out` receives a pointer to the new handle.
/// Free with `siderust_propagator_free`.
#[no_mangle]
pub extern "C" fn siderust_propagator_two_body_earth_new(
    out: *mut *mut SiderustPropagator,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let handle = Box::new(SiderustPropagator { force: TwoBody { mu: GM_EARTH } });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Create a DOP853 two-body propagator with a custom gravitational parameter.
///
/// - `gm_km3_s2` — gravitational parameter GM in km³/s² (e.g. 398 600.441 8
///   for Earth, 1.327 124 4e11 for the Sun).
///
/// On success, `*out` receives a pointer to the new handle.
/// Free with `siderust_propagator_free`.
#[no_mangle]
pub extern "C" fn siderust_propagator_two_body_new(
    gm_km3_s2: f64,
    out: *mut *mut SiderustPropagator,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let handle = Box::new(SiderustPropagator {
            force: TwoBody { mu: GravitationalParameter::new(gm_km3_s2) },
        });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Free a propagator handle previously returned by this crate.
///
/// # Safety
///
/// `handle` must be either null or a live pointer produced by this crate.
#[no_mangle]
pub unsafe extern "C" fn siderust_propagator_free(handle: *mut SiderustPropagator) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Propagate an orbit state forward (or backward) by `dt_s` seconds.
///
/// - `handle` — propagator created by `siderust_propagator_two_body_*_new`.
/// - `state_in` — initial orbit state (Geocentric/GCRS).
/// - `dt_s` — propagation interval in seconds (negative for backward propagation).
/// - `ctx` — dynamics context (may be null; null is treated as an empty context,
///   equivalent to `siderust_dynamics_context_new` with no providers attached).
/// - `out_state` — on success, `*out_state` is set to a newly allocated orbit
///   state representing the final state.  Free with `siderust_orbit_state_free`.
///
/// Returns [`SiderustDynamicsStatus`]:
/// - `OK` (0) on success.
/// - A specific non-zero code when propagation fails (see [`SiderustDynamicsStatus`]).
#[no_mangle]
pub extern "C" fn siderust_propagator_propagate(
    handle: *const SiderustPropagator,
    state_in: *const SiderustOrbitState,
    dt_s: f64,
    ctx: *const SiderustDynamicsContext,
    out_state: *mut *mut SiderustOrbitState,
) -> SiderustDynamicsStatus {
    dyn_guard! {{
        if handle.is_null() || state_in.is_null() || out_state.is_null() {
            return SiderustDynamicsStatus::NullPointer;
        }

        let propagator = unsafe { &(*handle) };
        let s0 = unsafe { (*state_in).inner };
        let dt = Second::new(dt_s);

        // Support null ctx as "empty context" so callers doing pure two-body
        // propagation do not need to allocate a context handle at all.
        let empty;
        let context: &DynamicsContext = if ctx.is_null() {
            empty = DynamicsContext::empty();
            &empty
        } else {
            unsafe { &(*ctx).inner }
        };

        // The propagator holds its own empty context internally; we need to
        // call the lower-level integration path that accepts an external context.
        // Use the force model directly via the integrator helpers.
        let result = propagate_with_context(propagator, s0, dt, context);

        match result {
            Ok(s_final) => {
                let handle_out = Box::new(SiderustOrbitState { inner: s_final });
                unsafe { *out_state = Box::into_raw(handle_out) };
                SiderustDynamicsStatus::Ok
            }
            Err(e) => SiderustDynamicsStatus::from_dynamics_error(&e),
        }
    }}
}

// =============================================================================
// Internal: propagate using an external DynamicsContext
// =============================================================================

/// Calls the DOP853 integrator directly, passing the caller-supplied context.
fn propagate_with_context(
    p: &SiderustPropagator,
    state: OrbitState<Geocentric, GCRS>,
    dt: Second,
    ctx: &DynamicsContext,
) -> Result<OrbitState<Geocentric, GCRS>, DynamicsError> {
    let tols = IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9);
    Ok(dop853_propagate(&p.force, state, dt, tols, ctx)?)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dynamics_status_ok_is_zero() {
        assert_eq!(SiderustDynamicsStatus::Ok as i32, 0);
    }

    #[test]
    fn dynamics_status_maps_ephemeris_unavailable() {
        let e = DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: None,
        };
        assert_eq!(
            SiderustDynamicsStatus::from_dynamics_error(&e),
            SiderustDynamicsStatus::EphemerisUnavailable
        );
    }

    #[test]
    fn dynamics_status_maps_gravity_field_unavailable() {
        let e = DynamicsError::GravityFieldUnavailable;
        assert_eq!(
            SiderustDynamicsStatus::from_dynamics_error(&e),
            SiderustDynamicsStatus::GravityFieldUnavailable
        );
    }
}
