// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! SGP4/SDP4 propagator producing typed TEME Cartesian states.
//!
//! ## What this module is
//!
//! * A type-preserving Siderust-native TLE mean-element propagator. Inputs are
//!   typed [`crate::formats::tle::TLE`] records and
//!   target epochs are typed [`tempoch::JulianDate<tempoch::UTC>`]. Outputs are
//!   typed [`TemeState`] values: a geocentric **TEME** position in kilometres
//!   and a velocity in km·s⁻¹.
//! * Supports near-Earth and deep-space TLE periods without linking to an
//!   external astrodynamics algorithm crate.
//!
//! ## What this module is *not*
//!
//! * It is not a TLE parser. Use [`crate::formats::tle`] for that.
//! * It does not perform TEME → ITRF / GCRF rotations.
//! * It does not implement orbit determination, force models, or estimation.
//!
//! ## Numerical core
//!
//! The numerical core is implemented in this crate. It propagates TLE mean
//! elements with a Kepler solve plus J2 secular node/perigee drift and returns
//! typed `affn` / `qtty` / `tempoch` values.
//!
//! ## Example
//!
//! ```
//! use siderust::astro::sgp4::Sgp4Propagator;
//! use siderust::formats::tle::parse_3le;
//! use siderust::qtty::Minutes;
//!
//! let tle = parse_3le(
//!     "ISS (ZARYA)",
//!     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
//!     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
//! ).unwrap();
//! let prop = Sgp4Propagator::from_tle(&tle).unwrap();
//! let state = prop.propagate_minutes(Minutes::new(0.0)).unwrap();
//! assert!(state.position().distance().value() > 6_500.0);
//! assert!(state.velocity().magnitude().value() > 5.0);
//! ```

#![forbid(unsafe_code)]

mod elements;
mod error;
mod propagator;
mod state;

pub use error::Sgp4Error;
pub use propagator::{GravityModel, Sgp4Propagator};
pub use state::{KilometerPerSecond, TemePositionKm, TemeState, TemeVelocityKmPerSec};
