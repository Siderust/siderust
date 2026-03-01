// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Ephemeris Provider Abstraction
//!
//! This module defines the [`Ephemeris`] trait (compile-time, static dispatch)
//! and the [`DynEphemeris`] trait (runtime, dynamic dispatch) for computing
//! fundamental solar-system body positions and velocities.
//!
//! ## Design
//!
//! - **Compile-time selection**: The [`Ephemeris`] trait uses associated
//!   functions (no `&self`) — backends are zero-sized marker types; all
//!   dispatch is monomorphized away.
//! - **Runtime selection**: The [`DynEphemeris`] trait uses `&self` methods,
//!   enabling trait objects (`Box<dyn DynEphemeris>`) for runtime-loaded data.
//!   A blanket impl bridges: every `Ephemeris` implementor automatically
//!   implements `DynEphemeris`.
//! - **Feature-gated**: VSOP87 is always available. DE440 requires
//!   the `de440` Cargo feature. Runtime loading requires `runtime-data`.
//! - **Body-centric user API preserved**: `Earth::vsop87e(jd)` etc. remain
//!   unchanged. The traits are used internally by the coordinate transform
//!   pipeline and `AstroContext`.
//!
//! ## Available Backends
//!
//! | Backend              | Feature        | Source                    |
//! |----------------------|----------------|---------------------------|
//! | [`Vsop87Ephemeris`]  | (always)       | VSOP87 + ELP2000-82B      |
//! | `De440Ephemeris`     | `de440`        | JPL DE440 (compile-time)   |
//! | [`RuntimeEphemeris`] | (always)       | Any BSP file (runtime)     |
//!
//! For DE441 (and other large datasets), use [`RuntimeEphemeris`] with
//! [`DataManager`](crate::data::DataManager) to download and load at runtime.

#[cfg(feature = "de440")]
mod de440_backend;
mod runtime_backend;
mod vsop87_backend;

#[cfg(feature = "de440")]
pub use de440_backend::De440Ephemeris;
pub use runtime_backend::RuntimeEphemeris;
pub use vsop87_backend::Vsop87Ephemeris;

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, Day, Kilometer};

/// Velocity measured in AU/day.
pub type AuPerDay = qtty::Per<AstronomicalUnit, Day>;

/// Trait abstracting over ephemeris backends for fundamental body states
/// (**compile-time, static dispatch**).
///
/// All methods are associated functions (no `&self` receiver) because backends
/// are zero-sized marker types — the coefficient data lives in static arrays
/// generated at compile time.
///
/// The five methods cover every ephemeris call site in the coordinate transform
/// pipeline, aberration, and observer-state modules.
///
/// See also [`DynEphemeris`] for the runtime/dynamic-dispatch counterpart.
pub trait Ephemeris {
    /// Sun's position in barycentric ecliptic coordinates (AU).
    fn sun_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's position in barycentric ecliptic coordinates (AU).
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's position in heliocentric ecliptic coordinates (AU).
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's velocity in barycentric ecliptic coordinates (AU/day).
    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay>;

    /// Moon's position in geocentric ecliptic coordinates (km).
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, EclipticMeanJ2000, Kilometer>;
}

/// Trait abstracting over ephemeris backends for fundamental body states
/// (**runtime, dynamic dispatch**).
///
/// This is the instance-based counterpart to [`Ephemeris`]. All methods take
/// `&self`, making this trait **object-safe** — you can use `Box<dyn DynEphemeris>`
/// or `&dyn DynEphemeris` for runtime-selected backends.
///
/// ## Blanket implementation
///
/// Every type implementing [`Ephemeris`] automatically implements `DynEphemeris`
/// (the `&self` receiver is ignored since `Ephemeris` uses associated functions).
/// This means `Vsop87Ephemeris`, `De440Ephemeris`, etc. can all be used as
/// `&dyn DynEphemeris` with zero additional code.
///
/// ## Example
///
/// ```rust,ignore
/// use siderust::calculus::ephemeris::{DynEphemeris, Vsop87Ephemeris, RuntimeEphemeris};
///
/// // Compile-time backend used dynamically:
/// let eph: &dyn DynEphemeris = &Vsop87Ephemeris;
/// let sun = eph.sun_barycentric(jd);
///
/// // Runtime-loaded backend:
/// let rt = RuntimeEphemeris::from_bsp("path/to/de441.bsp")?;
/// let sun = rt.sun_barycentric(jd);
/// ```
pub trait DynEphemeris {
    /// Sun's position in barycentric ecliptic coordinates (AU).
    fn sun_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's position in barycentric ecliptic coordinates (AU).
    fn earth_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's position in heliocentric ecliptic coordinates (AU).
    fn earth_heliocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Earth's velocity in barycentric ecliptic coordinates (AU/day).
    fn earth_barycentric_velocity(&self, jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay>;

    /// Moon's position in geocentric ecliptic coordinates (km).
    fn moon_geocentric(&self, jd: JulianDate)
        -> Position<Geocentric, EclipticMeanJ2000, Kilometer>;
}

/// Blanket implementation: every static `Ephemeris` backend automatically
/// implements `DynEphemeris`.
impl<T: Ephemeris> DynEphemeris for T {
    #[inline]
    fn sun_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        <T as Ephemeris>::sun_barycentric(jd)
    }

    #[inline]
    fn earth_barycentric(
        &self,
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        <T as Ephemeris>::earth_barycentric(jd)
    }

    #[inline]
    fn earth_heliocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        <T as Ephemeris>::earth_heliocentric(jd)
    }

    #[inline]
    fn earth_barycentric_velocity(&self, jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        <T as Ephemeris>::earth_barycentric_velocity(jd)
    }

    #[inline]
    fn moon_geocentric(
        &self,
        jd: JulianDate,
    ) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        <T as Ephemeris>::moon_geocentric(jd)
    }
}
