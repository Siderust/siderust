// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Time scales and instants
//!
//! ## Scientific scope
//!
//! Astronomical computations live on top of a small zoo of time scales —
//! TT (Terrestrial Time), TDB and TCB (the dynamical scales used by
//! barycentric ephemerides), TCG (geocentric coordinate time), TAI
//! (atomic time), UTC (civil time, with leap seconds), and UT1 (the
//! Earth-rotation scale). Mixing them silently produces errors of
//! seconds to tens of seconds — large enough to invalidate any
//! sub-arcsecond positional result. Following the IAU SOFA conventions
//! and the IERS Conventions (2010), every astronomical instant in this
//! crate carries its time scale at the type level so that conversions
//! between scales must be explicit.
//!
//! The reference epoch used throughout is **J2000.0 = JD 2 451 545.0
//! TT** (= 2000 January 1, 12:00:00 TT), as defined by IAU 2006
//! Resolution B1.
//!
//! ## Technical scope
//!
//! All time handling is delegated to the `tempoch` crate. This module
//! re-exports its full public surface (`Scale`, `EncodedTime`, `Time`,
//! the `TT`/`TDB`/`TCB`/`TCG`/`TAI`/`UTC`/`UT1` scale markers,
//! `Interval`, `delta_t_seconds`, …) and defines the v1 astronomy-facing
//! defaults as TT-based encoded dates plus the named same-scale helpers
//! (`to_jd`, `to_mjd`, `to_j2000_seconds`, `shifted_by`, `duration_since`)
//! provided by tempoch core:
//!
//! - [`JulianDate`] = `tempoch::JulianDate<TT>` = `EncodedTime<TT, JD>`.
//! - [`ModifiedJulianDate`] = `tempoch::ModifiedJulianDate<TT>` =
//!   `EncodedTime<TT, MJD>`.
//!
//! Advanced scale-aware work still uses `tempoch::JulianDate<S>` and
//! `tempoch::ModifiedJulianDate<S>` directly. The [`J2000`] constant exposes
//! the J2000.0 TT epoch as a [`JulianDate`]. [`UT`] is a backward-compatible
//! alias for [`UT1`]; [`Period<T>`] is a backward-compatible alias for
//! [`Interval<T>`].
//!
//! ## References
//!
//! - IAU SOFA Board (2021). *IAU SOFA Software Collection*. International
//!   Astronomical Union. <http://www.iausofa.org/>.
//! - Petit, G., Luzum, B. (Eds.) (2010). *IERS Conventions (2010)*. IERS
//!   Technical Note 36. Verlag des Bundesamts für Kartographie und
//!   Geodäsie, Frankfurt am Main.
//! - Seidelmann, P. K. (Ed.) (1992). *Explanatory Supplement to the
//!   Astronomical Almanac*. University Science Books. ISBN 0-935702-68-7.
//! - International Astronomical Union (2006). Resolution B3 ("Re-definition
//!   of Barycentric Dynamical Time, TDB"). IAU Transactions XXVIB.

pub use tempoch::{
    complement_within, constats, delta_t_seconds, delta_t_seconds_extrapolated, eop,
    ContinuousScale, ConversionError, ConversionTarget, CoordinateScale, EncodedTime,
    FormatForScale, InfallibleConversionTarget, InfallibleFormatForScale, Interval,
    InvalidIntervalError, InvalidPeriodError, PeriodListError, Scale, ScaleKind, Time, TimeContext,
    TimeDataError, TimeInstant, JD, JULIAN_YEAR_DAYS, MJD, TAI, TCB, TCG, TDB, TT, UT1, UTC,
};

/// Backward-compatible alias: old siderust code used `UT` for the UT1 axis.
pub type UT = UT1;

/// Backward-compatible generic period alias over an instant type `T`.
pub type Period<T> = Interval<T>;

/// TT-scale Julian Date.  Equivalent to `tempoch::JulianDate<TT>`.
pub type JulianDate = tempoch::JulianDate<TT>;

/// TT-scale Modified Julian Date.  Equivalent to `tempoch::ModifiedJulianDate<TT>`.
pub type ModifiedJulianDate = tempoch::ModifiedJulianDate<TT>;

/// J2000.0 epoch as a [`JulianDate`] (TT, JD 2 451 545.0).
pub const J2000: JulianDate = JulianDate::J2000;

/// Checked TT Julian Date constructor from a typed day quantity.
#[inline]
pub fn try_jd(raw: qtty::Day) -> Result<JulianDate, ConversionError> {
    JulianDate::try_new(raw)
}

/// Checked TT Julian Date constructor from a scalar JD value.
#[inline]
pub fn try_jd_f64(value: f64) -> Result<JulianDate, ConversionError> {
    try_jd(qtty::Day::new(value))
}

/// Panicking TT Julian Date constructor for internal finite constants and tests.
#[track_caller]
#[inline]
pub fn jd(raw: qtty::Day) -> JulianDate {
    try_jd(raw).expect("JulianDate must be finite")
}

/// Checked TT Modified Julian Date constructor from a typed day quantity.
#[inline]
pub fn try_mjd(raw: qtty::Day) -> Result<ModifiedJulianDate, ConversionError> {
    ModifiedJulianDate::try_new(raw)
}

/// Checked TT Modified Julian Date constructor from a scalar MJD value.
#[inline]
pub fn try_mjd_f64(value: f64) -> Result<ModifiedJulianDate, ConversionError> {
    try_mjd(qtty::Day::new(value))
}

/// Panicking TT Modified Julian Date constructor for internal finite constants and tests.
#[track_caller]
#[inline]
pub fn mjd(raw: qtty::Day) -> ModifiedJulianDate {
    try_mjd(raw).expect("ModifiedJulianDate must be finite")
}

/// Helper: build a `ModifiedJulianDate` (TT scale) directly from a `chrono::DateTime<Utc>`.
///
/// This mirrors the common pattern `Time::<UTC>::from_chrono(dt).into()` but provides
/// a concise, single-call helper for callers that only need an MJD.
#[inline]
pub fn modified_julian_date_from_chrono(dt: chrono::DateTime<UTC>) -> ModifiedJulianDate {
    Time::<UTC>::from_chrono(dt).into()
}

pub use crate::calculus::math_core::intervals::intersect as intersect_periods;
