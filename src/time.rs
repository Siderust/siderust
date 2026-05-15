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
//! `JulianDate`, `ModifiedJulianDate`, the `TT`/`TDB`/`TCB`/`TCG`/
//! `TAI`/`UTC`/`UT1` scale markers, `Interval`, `delta_t_seconds`, …)
//! and adds two backward-compatible non-generic aliases that pin the
//! scale to TT for legacy call sites:
//!
//! - [`JulianDate`] = `tempoch::JulianDate<TT>` = `EncodedTime<TT, JD>`.
//! - [`ModifiedJulianDate`] = `tempoch::ModifiedJulianDate<TT>` =
//!   `EncodedTime<TT, MJD>`.
//!
//! Generic re-exports `JulianDateG`/`ModifiedJulianDateG` keep the
//! `JulianDate<TT>` / `JulianDate<TDB>` form compiling. The
//! [`J2000`] constant exposes the J2000.0 TT epoch as a [`JulianDate`].
//! [`UT`] is a backward-compatible alias for [`UT1`]; [`Period<T>`] is a
//! backward-compatible alias for [`Interval<T>`].
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

/// Re-export `JulianDate` and `ModifiedJulianDate` under their generic names
/// so code that explicitly writes `JulianDate<TT>` or `JulianDate<TDB>` still compiles.
pub use tempoch::{JulianDate as JulianDateG, ModifiedJulianDate as ModifiedJulianDateG};

/// Backward-compatible alias: old siderust code used `UT` for the UT1 axis.
pub type UT = UT1;

/// Backward-compatible generic period alias over an instant type `T`.
pub type Period<T> = Interval<T>;

/// TT-scale Julian Date.  Equivalent to `tempoch::JulianDate<TT>`.
pub type JulianDate = tempoch::JulianDate<TT>;

/// TT-scale Modified Julian Date.  Equivalent to `tempoch::ModifiedJulianDate<TT>`.
pub type ModifiedJulianDate = tempoch::ModifiedJulianDate<TT>;

/// J2000.0 epoch as a [`JulianDate`] (TT, JD 2 451 545.0).
pub const J2000: JulianDate =
    JulianDate::from_raw_unchecked(qtty::Quantity::<qtty::unit::Day>::new(2_451_545.0));

pub use crate::calculus::math_core::intervals::intersect as intersect_periods;
