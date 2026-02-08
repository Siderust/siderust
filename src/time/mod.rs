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
//! - [`ModifiedJulianDate`] — type alias for `Time<MJD>`.
//! - [`Period<T>`] — a time interval between two [`TimeInstant`]s.
//! - [`TimeInstant`] — trait for points in time usable with [`Period`].
//!
//! # Time scales
//!
//! The following markers implement [`TimeScale`]:
//!
//! | Marker | Scale |
//! |--------|-------|
//! | [`JD`] | Julian Date |
//! | [`MJD`] | Modified Julian Date |
//! | [`TDB`] | Barycentric Dynamical Time |
//! | [`TT`] | Terrestrial Time |
//! | [`TAI`] | International Atomic Time |
//! | [`GPS`] | GPS Time |
//! | [`UnixTime`] | Unix / POSIX time |

pub(crate) mod instant;
mod julian_date_ext;
mod period;
pub(crate) mod scales;

// ── Re-exports ────────────────────────────────────────────────────────────

pub use instant::{Time, TimeInstant, TimeScale};
pub use period::{complement_within, intersect_periods, Period};
pub use scales::{UnixTime, GPS, JD, MJD, TAI, TDB, TT};

// ── Backward-compatible type aliases ──────────────────────────────────────

/// Julian Date — continuous count of days since the Julian Period.
///
/// This is a type alias for [`Time<JD>`].  All historical call-sites
/// (`JulianDate::new(...)`, `JulianDate::J2000`, `.julian_centuries()`, …)
/// continue to work without modification.
pub type JulianDate = Time<JD>;

/// Modified Julian Date — `JD − 2 400 000.5`.
///
/// This is a type alias for [`Time<MJD>`].
pub type ModifiedJulianDate = Time<MJD>;
