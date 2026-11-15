// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Calculus Module
//!
//! Lunar geometry, phase, brightness, and event searches. The standard
//! altitude and azimuth APIs use geometric sightlines from an Earth
//! observer. The default ephemeris determines position accuracy, and the
//! Earth-orientation context supplies the rotation to terrestrial axes.
//!
//! Altitude searches cache geocentric Cartesian positions and precompute
//! the site's Earth-fixed position and local basis. They evaluate Earth
//! orientation at each query. The crossing search uses Chebyshev candidate
//! discovery, refinement, and a local scan fallback.

//! ## References
//! - Chapront‑Touzé, M. & Chapront, J. (1983). "The lunar ephemeris ELP 2000".
//!   *A&A*, 124, 50–62.
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann‑Bell.
//! - Folkner, W. M. et al. (2014). "The Planetary and Lunar Ephemeris DE 430
//!   and DE 431". *IPN Progress Report* 42‑196.
//! - Park, R. S. et al. (2021). "The JPL Planetary and Lunar Ephemerides
//!   DE 440 and DE 441". *AJ*, 161, 105.
//! - Jones, A. (2013). Lunar reflectance / phase attenuation model
//!   (used by [`photometry`]).

mod altitude_periods;
pub(crate) mod azimuth;
pub mod meeus_ch47;
pub(crate) mod moon_cache;
mod moon_equations;
pub mod phase;
pub mod photometry;
pub mod position_angle;

pub(crate) use altitude_periods::*;
pub(crate) use azimuth::*;
