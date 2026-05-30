// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IERS Earth Orientation Parameters
//!
//! Re-exports the IERS EOP data and lookup API from [`tempoch`].
//!
//! The EOP table (daily `finals2000A.all` values for polar motion `xp`, `yp`,
//! `UT1 − UTC`, and celestial-pole offsets `dX`, `dY`) is generated and
//! maintained by `tempoch`.  Siderust no longer embeds its own copy.
//!
//! ## Usage
//!
//! ```no_run
//! use siderust::astro::iers_data::{builtin_eop_at, builtin_eop_covers};
//! use siderust::qtty::Days;
//!
//! // MJD 51544.5 ≈ J2000.0 in UTC
//! let mjd = Days::new(51544.5);
//! if builtin_eop_covers(mjd) {
//!     let eop = builtin_eop_at(mjd).unwrap();
//!     let _ = eop.ut1_minus_utc;
//! }
//! ```
//!
//! For the full typed EOP provider, use [`crate::astro::eop::IersEop`].

pub use tempoch::eop::{
    builtin_eop_at, builtin_eop_covers, eop_end, eop_observed_end, eop_start, EopValues,
};
