// SPDX-License-Identifier: AGPL-3.0-only
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
//! This module is a **thin delegation facade** over `tempoch`. It re-exports
//! `tempoch`'s full public surface (scales, `Time`, `TimeContext`, format tags,
//! `Interval`, ΔT/EOP accessors) and defines two siderust-specific shortcuts:
//!
//! - [`JulianDate`] = `tempoch::JulianDate<TT>` (TT-default Julian Date)
//! - [`ModifiedJulianDate`] = `tempoch::ModifiedJulianDate<TT>`
//! - [`J2000`] — the J2000.0 TT epoch constant
//!
//! All time-scale ownership (ΔT tables, EOP, UTC/TAI leap seconds, time-data
//! freshness) belongs to `tempoch`.  Callers that need advanced scale-aware
//! work can use `tempoch` types directly; those types are identical to what
//! `siderust::time::*` re-exports.
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

pub use tempoch::*;

/// TT-scale Julian Date.  Equivalent to `tempoch::JulianDate<TT>`.
pub type JulianDate = tempoch::JulianDate<TT>;

/// TT-scale Modified Julian Date.  Equivalent to `tempoch::ModifiedJulianDate<TT>`.
pub type ModifiedJulianDate = tempoch::ModifiedJulianDate<TT>;

/// J2000.0 TT as a [`JulianDate`] (`JD 2 451 545.0 TT`); **`const`** everywhere.
pub const J2000: JulianDate = tempoch::JulianDate::<TT>::JD_EPOCH_J2000_0;

pub use crate::event::search::intervals::intersect as intersect_periods;
