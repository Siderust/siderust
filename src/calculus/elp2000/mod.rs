// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ELP2000 Lunar Theory Module
//!
//! This module implements the ELP2000-82B theory, a high-precision analytical model for the motion of the Moon.
//!
//! ## What is ELP2000?
//!
//! ELP2000 (Éphémérides Lunaires Parisiennes 2000) is a semi-analytical lunar theory developed by M. Chapront-Touzé and J. Chapront at the Bureau des Longitudes, Paris. It provides a mathematical model for the Moon's orbital motion, including its position and velocity, with high accuracy over long time spans.
//!
//! The "2000" in ELP2000 refers to the reference epoch J2000.0 (Julian Date 2451545.0, January 1, 2000, 12:00 TT).
//!
//! ## Where and Why is ELP2000 Used?
//!
//! ELP2000 is widely used in astronomy, space science, and navigation for:
//! - Computing precise lunar ephemerides (positions and velocities).
//! - Predicting lunar eclipses, occultations, and other astronomical events.
//! - Spacecraft navigation and mission planning involving the Moon.
//! - Historical and future lunar position calculations.
//!
//! ## Mathematical, Astronomical, and Algorithmic Background
//!
//! ELP2000 is based on the analytical expansion of the Moon's motion using Poisson series, taking into account:
//! - Perturbations from the Sun, Earth, and planets.
//! - Relativistic corrections.
//! - Tidal effects and secular variations.
//!
//! The theory expresses the Moon's longitude, latitude, and distance as sums of periodic terms, each with specific amplitudes, frequencies, and phases. The coefficients are derived from observational data and numerical integrations.
//!
//! The implementation typically involves:
//! - Loading large sets of coefficients (the "ELP data").
//! - Evaluating trigonometric series for a given time (usually in Julian centuries since J2000.0).
//! - Applying corrections for nutation, precession, and other effects as needed.
//!
//! ## References
//!
//! - Chapront-Touzé, M., & Chapront, J. (1983). "The lunar ephemeris ELP 2000". Astronomy and Astrophysics, 124, 50-62.
//! - Chapront-Touzé, M., & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris adequate for historical times". Astronomy and Astrophysics, 190, 342-352.
//! - [Bureau des Longitudes ELP2000 Resources](https://www.imcce.fr/Equipes/ASD/insola/elp/index.html)
//! - Meeus, J. (1998). "Astronomical Algorithms", 2nd Edition. Willmann-Bell, Inc.
//!
//! ## Module Structure
//!
//! - `elp_structs`: Data structures for ELP2000 coefficients and series.
//! - `elp_data`: The actual numerical coefficients and tables for the theory.
//! - `elp2000`: Algorithms for evaluating the lunar position and velocity using ELP2000.
//!
//! This module provides the foundation for high-precision lunar calculations in Siderust.

mod elp_series;
mod elp_structs;
