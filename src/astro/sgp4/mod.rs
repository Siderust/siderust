// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! TLE mean-element propagator producing typed TEME Cartesian states.
//!
//! Full Vallado SGP4 (near-Earth) and SDP4 (deep-space) propagation, including
//! secular/periodic corrections, BSTAR drag, lunar/solar terms and resonances.
//! The numerical core is native Rust, following Vallado et al.,
//! [AIAA 2006-6753](https://celestrak.org/publications/AIAA/2006-6753/).
//! It uses the existing affn/qtty/tempoch stack, with no SGP4 dependency.
//!
//! [`GravityModel::Wgs72`] selects WGS-72 and AFSPC initialization **and**
//! propagation conventions. Vallado SGP4-VER vectors are tested through the
//! typed API at absolute component tolerances of 1e-6 km and 1e-9 km/s.
//! [`GravityModel::Wgs72Iau`] explicitly keeps WGS-72 while using IAU sidereal
//! time and epoch conventions. [`GravityModel::Wgs84`] uses WGS-84 with IAU.
//! These alternatives are not expected to reproduce the WGS-72/AFSPC vectors.
//!
//! Inputs are typed TLE **mean** elements, not osculating Kepler elements.
//! Outputs are geocentric TEME kilometres and kilometres per second; use the
//! independent-vector transformation infrastructure for other frames (and EOP data
//! where required). This replaces the previous secular J2 approximation;
//! numerical results intentionally change, including at epoch.
//!
//! Time follows 1,440 nominal minutes per civil UTC day, rather than elapsed SI
//! seconds. Tempoch performs civil UTC validation and conversion. Leap-second
//! labels cannot be represented on this grid and are rejected. Calendar
//! conversion preserves nanoseconds and avoids large scalar JD subtraction.
//! A scalar JD supplied by a caller has already lost precision (about 40 µs
//! per ULP near modern epochs); this cannot be recovered. Offsets must be finite
//! and fit the calendar conversion's signed 64-bit nanosecond range
//! (approximately ±292 years); actual useful TLE prediction spans are far shorter.
//!
//! Invalid elements, unrepresentable epochs and divergent propagation return
//! [`Sgp4Error`]. As in the reference model, sub-Earth positions are possible:
//! this API is not a reentry detector. TLE derivatives are metadata; SGP4 uses
//! mean motion and BSTAR rather than independently integrating those derivatives.
//!
//! ## Example
//!
//! ```
//! use siderust::astro::sgp4::TlePropagator;
//! use siderust::formats::tle::parse_3le;
//! use siderust::qtty::Minutes;
//!
//! let tle = parse_3le(
//!     "ISS (ZARYA)",
//!     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
//!     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
//! ).unwrap();
//! let prop = TlePropagator::from_tle(&tle).unwrap();
//! let state = prop.propagate_minutes(Minutes::new(0.0)).unwrap();
//! assert!(state.position().distance().value() > 6_500.0);
//! assert!(state.velocity().magnitude().value() > 5.0);
//! ```
//!
//! Full coordinate-velocity conversion between time-dependent frames additionally
//! needs `Rdot r` and, for translating origins, origin velocity. The frame
//! rotation providers and `BlockDiagRotation6` alone do not provide that map.

#![forbid(unsafe_code)]

mod deep_space;
mod elements;
mod error;
mod model;
mod propagator;
mod state;
mod time;

pub use error::Sgp4Error;
pub use propagator::{GravityModel, TlePropagator};
pub use state::{KilometerPerSecond, TemePositionKm, TemeState, TemeVelocityKmPerSec};

/// Backward-compatible name for the TLE mean-element propagator.
///
/// New code may use [`TlePropagator`] when the TLE input contract is the
/// important part of the API, or `Sgp4Propagator` when matching existing SGP4
/// example code.
pub type Sgp4Propagator = TlePropagator;
