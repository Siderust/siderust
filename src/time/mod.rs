// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Time Module
//!
//! This module provides time-related types and abstractions for astronomical calculations.
//!
//! # Core types
//!
//! - [`Time<S>`] — generic instant parameterised by a [`TimeScale`] marker.
//! - [`TimeScale`] — trait that defines a time scale (epoch offset + conversions).
//! - [`JulianDate`] — type alias for `Time<JD>`.
//! - [`JulianEphemerisDay`] — type alias for `Time<JDE>`.
//! - [`ModifiedJulianDate`] — type alias for `Time<MJD>`.
//! - [`Period<S>`] — a time interval parameterised by a [`TimeScale`] marker.
//! - [`Interval<T>`] — a generic interval over any [`TimeInstant`].
//! - [`TimeInstant`] — trait for points in time usable with [`Interval`].
//!
//! # Time scales
//!
//! The following markers implement [`TimeScale`]:
//!
//! | Marker | Scale |
//! |--------|-------|
//! | [`JD`] | Julian Date |
//! | [`JDE`] | Julian Ephemeris Day |
//! | [`MJD`] | Modified Julian Date |
//! | [`TDB`] | Barycentric Dynamical Time |
//! | [`TT`] | Terrestrial Time |
//! | [`TAI`] | International Atomic Time |
//! | [`GPS`] | GPS Time |
//! | [`UnixTime`] | Unix / POSIX time |
//! | [`UT`] | Universal Time (Earth rotation) |
//!
//! # ΔT (Delta T)
//!
//! The difference **ΔT = TT − UT** is applied automatically by the
//! [`UT`] time scale.  Use `Time::<UT>::new(jd_ut)` for UT-based values,
//! or construct any scale via `from_utc()` which routes through `UT` internally.
//! The raw ΔT value (in seconds) is available via [`Time::<UT>::delta_t()`](Time::delta_t).

mod delta_t;
pub(crate) mod instant;
mod julian_date_ext;
mod period;
pub(crate) mod scales;

// ── Re-exports ────────────────────────────────────────────────────────────

pub use instant::{Time, TimeInstant, TimeScale};
pub use period::{complement_within, intersect_periods, Interval, Period, UtcPeriod};
pub use scales::{UnixTime, GPS, JD, JDE, MJD, TAI, TDB, TT, UT};

// ── Backward-compatible type aliases ──────────────────────────────────────

/// Julian Date — continuous count of days since the Julian Period.
///
/// This is a type alias for [`Time<JD>`].  All historical call-sites
/// (`JulianDate::new(...)`, `JulianDate::J2000`, `.julian_centuries()`, …)
/// continue to work without modification.
pub type JulianDate = Time<JD>;

/// Julian Ephemeris Day — dynamical Julian day used by many ephemeris formulas.
///
/// This is a type alias for [`Time<JDE>`].
pub type JulianEphemerisDay = Time<JDE>;

/// Modified Julian Date — `JD − 2 400 000.5`.
///
/// This is a type alias for [`Time<MJD>`].
pub type ModifiedJulianDate = Time<MJD>;

/// Universal Time — Earth-rotation civil time scale.
///
/// This is a type alias for [`Time<UT>`].
pub type UniversalTime = Time<UT>;
