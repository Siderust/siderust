// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # JPL Planetary and Lunar Ephemerides (DE4xx)
//!
//! ## Scientific scope
//!
//! The JPL DE (*Development Ephemeris*) series are numerical integrations of
//! the equations of motion for the Solar System bodies produced by the Jet
//! Propulsion Laboratory.  Each release (DE440, DE441, …) provides highly
//! accurate Chebyshev polynomial representations of barycentric and
//! heliocentric positions and velocities.
//!
//! ## Technical scope
//!
//! - [`eval`] — Chebyshev polynomial evaluation and [`eval::DynSegmentDescriptor`],
//!   which stores a JD interval and its coefficient block, and evaluates position /
//!   velocity for one body.
//! - [`bodies`] — generic body-chain resolution that derives Earth, Sun, and Moon
//!   from the natively integrated barycentric states (Earth–Moon Barycenter +
//!   Moon offset).
//!
//! DE4xx files are loaded at runtime via [`crate::ephemeris::RuntimeEphemeris`].
//!
//! ## References
//!
//! - Standish, E. M. (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405".
//!   *JPL Interoffice Memorandum* 312.F-98-048.
//! - Folkner, W. M., Williams, J. G., Boggs, D. H., Park, R. S., & Kuchynka, P.
//!   (2014). "The Planetary and Lunar Ephemerides DE430 and DE431".
//!   *IPN Progress Report* 42-196, 1–81.
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440
//!   and DE441". *The Astronomical Journal* 161, 105.
//!   <https://doi.org/10.3847/1538-3881/abd414>

pub(crate) mod bodies;
pub(crate) mod eval;
