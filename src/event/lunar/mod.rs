// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Calculus Module
//!
//! ## Scientific scope
//!
//! Apparent geometry and event calculus of the Moon as seen from a
//! topocentric Earth observer. Covers position, altitude, azimuth,
//! illuminated fraction and phase events, and standard photometric
//! attenuation models. Position accuracy is bounded by the chosen engine:
//! ELP2000‑82B / Meeus chapter 47 truncations for the analytical path,
//! JPL DE440/441 Chebyshev coefficients for the high‑precision path.
//! Atmospheric refraction is not applied.
//!
//! ## Technical scope
//!
//! - `altitude_periods`: Moon altitude band/threshold periods.
//! - `azimuth`: scalar azimuth function used by the unified azimuth API.
//! - `meeus_ch47`: closed‑form Moon position from *Astronomical Algorithms*.
//! - `moon_cache`: Chebyshev caches for fast batch evaluation of
//!   geocentric Moon position, nutation rotation, and altitude.
//! - `moon_equations`: high‑level `Moon::*` API on top of the engines.
//! - `phase`: phase angle, illuminated fraction, phase events
//!   (new/first‑quarter/full/last‑quarter), illumination period finders.
//! - `photometry`: empirical lunar reflectance / phase attenuation models.
//!
//! All period‑finding delegates to `crate::event::search::intervals`
//! which provides scan + Brent refinement + interval assembly. This module
//! supplies the Moon‑altitude closure and JD↔MJD conversions.
//!
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

pub(crate) use altitude_periods::*;
pub(crate) use azimuth::*;
