// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Stellar Altitude Periods
//!
//! Efficient algorithms for finding time intervals where a **fixed star**
//! (static ICRS direction) is above, below, or within a given altitude range
//! as seen from an observer on Earth.
//!
//! ## Key Insight
//!
//! A star's altitude varies sinusoidally with Earth's rotation — all stars
//! share the same period (one sidereal day) and differ only in amplitude,
//! offset, and phase.  This module exploits that structure for **analytical
//! bracket discovery**, avoiding the expensive uniform scan used for bodies
//! with non‑trivial orbital motion (Sun, Moon).
//!
//! ## Usage
//!
//! ```ignore
//! use siderust::calculus::stellar::*;
//! use qtty::*;
//!
//! let ra  = Degrees::new(101.287);   // Sirius RA  (J2000)
//! let dec = Degrees::new(-16.716);   // Sirius Dec (J2000)
//!
//! let periods = find_star_above_periods(ra, dec, site, window, Degrees::new(0.0));
//! ```
//!
//! ## See Also
//!
//! - [`crate::calculus::altitude`] — unified API (dispatches here for
//!   [`AltitudeTarget::FixedEquatorial`](crate::calculus::altitude::AltitudeTarget::FixedEquatorial))
//! - [`crate::calculus::solar`] — analogous module for the Sun
//! - [`crate::calculus::lunar`] — analogous module for the Moon

pub(crate) mod star_equations;

pub mod altitude_periods;

pub use altitude_periods::*;
