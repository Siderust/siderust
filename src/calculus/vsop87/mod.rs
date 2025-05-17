//! # VSOP87 Planetary Theory Module
//!
//! This module provides high-precision planetary positions using the VSOP87 theory.
//!
//! ## What is VSOP87?
//!
//! VSOP87 (Variations Séculaires des Orbites Planétaires 1987) is a semi-analytical planetary theory developed by Pierre Bretagnon and Gérard Francou at the Bureau des Longitudes, Paris. It models the motion of the major planets in the Solar System using trigonometric series expansions, providing accurate heliocentric and barycentric positions over long time spans.
//!
//! The theory is available in several versions (VSOP87A, VSOP87B, VSOP87C, VSOP87D, VSOP87E), each tailored for different coordinate systems and applications. This module focuses on the most widely used forms: VSOP87A (heliocentric ecliptic rectangular coordinates, equinox J2000.0) and VSOP87E (barycentric ecliptic rectangular coordinates, equinox J2000.0).
//!
//! ## Where and Why is VSOP87 Used?
//!
//! VSOP87 is a standard in astronomy and space science for:
//! - Computing planetary ephemerides (positions and velocities).
//! - Predicting planetary conjunctions, oppositions, and transits.
//! - Spacecraft navigation and mission planning.
//! - Astronomical software and almanacs (e.g., Jean Meeus' "Astronomical Algorithms").
//!
//! ## Mathematical, Astronomical, and Algorithmic Background
//!
//! VSOP87 expresses the position of each planet as a sum of periodic terms (trigonometric series) in time, with coefficients derived from numerical integrations and observations. The expansions account for:
//! - Mutual gravitational perturbations between planets.
//! - Secular variations and long-term trends.
//! - High-order periodic effects.
//!
//! The general form for each coordinate (X, Y, Z) is:
//! ```text
//! X = Σ (A * cos(B + C * t))
//! ```
//! where `A`, `B`, and `C` are coefficients, and `t` is time in Julian millennia from J2000.0.
//!
//! The implementation involves:
//! - Loading large tables of coefficients for each planet and coordinate.
//! - Evaluating the series for a given Julian date.
//! - Supporting both heliocentric and barycentric reference frames.
//!
//! ## References
//!
//! - Bretagnon, P., & Francou, G. (1988). "VSOP87: Theory of the motions of the five major planets". Astronomy and Astrophysics, 202, 309-315.
//! - [IMCCE VSOP87 Resources](https://www.imcce.fr/inpop/ephemerides/vsop87/)
//! - Meeus, J. (1998). "Astronomical Algorithms", 2nd Edition. Willmann-Bell, Inc.
//!
//! ## Module Structure
//!
//! - `vsop87_trait`: Trait abstraction for VSOP87 implementations.
//! - `vsop87a`: Implementation for VSOP87A (heliocentric).
//! - `vsop87e`: Implementation for VSOP87E (barycentric).
//! - `vsop87`: Core algorithms and shared logic.
//!
//! This module enables high-precision planetary calculations for Siderust and related astronomical applications.

mod vsop87_impl;

#[allow(clippy::all)]
mod vsop87a;
#[allow(clippy::all)]
mod vsop87e;

mod vsop87_trait;

pub use vsop87_trait::VSOP87;

use vsop87_impl::Vsop87;
use vsop87_impl::compute_vsop87;
