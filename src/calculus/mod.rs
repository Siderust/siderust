// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Celestial Mechanics Calculus Module
//!
//! This module provides the core algorithms and data structures for high-precision
//! celestial mechanics calculations in siderust. It is the computational backbone
//! for planetary, lunar, and solar ephemerides, as well as for general orbital mechanics.
//!
//! ## Overview
//!
//! The `calculus` module is organized into several submodules, each responsible for
//! a specific aspect of astronomical computation. The main functionalities include:
//!
//! - **VSOP87 planetary theory**: High-precision positions for the major planets.
//! - **ELP2000 lunar theory**: Accurate lunar positions using the ELP2000-82B model.
//! - **Solar calculations**: Solar coordinates and related phenomena.
//! - **Keplerian orbits**: General elliptical orbit propagation and Kepler's equation solvers.
//! - **Extremas finding**: Algorithms to find maxima/minima (e.g., culmination, rise/set).
//!
//! ## Submodules
//!
//! ### `vsop87`
//! Implements the VSOP87 theory for planetary positions. This includes:
//! - Data tables for VSOP87A/E (and potentially other variants).
//! - Efficient evaluation of planetary heliocentric coordinates (X, Y, Z) as a function of time.
//! - Traits and structures for extensibility and abstraction over different VSOP versions.
//!
//! ### `elp2000`
//! Implements the ELP2000-82B theory for the Moon's position. Features:
//! - Series expansions for lunar longitude, latitude, and distance.
//! - Precession corrections and planetary perturbations.
//! - Functions to compute geocentric ecliptic coordinates of the Moon for any Julian date.
//!
//! ### `solar`
//! Contains algorithms for solar coordinates and related calculations, such as:
//! - Apparent and mean solar longitude/latitude.
//! - Solar equations of time, declination, and related phenomena.
//! - Utility functions for solar ephemerides.
//!
//! ### `kepler_equations`
//! Provides general-purpose routines for elliptical orbit propagation, including:
//! - Solvers for Kepler's equation (Newton-Raphson and bisection methods).
//! - Calculation of true anomaly, eccentric anomaly, and orbital positions.
//! - Conversion between orbital elements and Cartesian coordinates.
//! - Orbital period and mean motion calculations.
//!
//! ### `events`
//! Algorithms to catch astro events such as:
//! - Finding upper/lower culmination (transit) times for celestial objects.
//! - Root-finding and refinement routines for event timing.
//! - Support for both static and dynamic targets (e.g., stars, planets).
//!
//! ## Usage
//!
//! This module is intended for internal use by higher-level ephemeris and observer modules,
//! but its public API can be used directly for custom calculations. For example:
//!
//! ```rust
//! use siderust::calculus::kepler_equations::solve_keplers_equation;
//! use qtty::*;
//!
//! let m = 1.0 * RAD;
//! let e = 0.0167;
//! let eccentric_anomaly = solve_keplers_equation(m, e);
//! ```
//!
//! ## Re-exports
//!
//! The module re-exports the most important functions and types from its submodules for convenience:
//! - `kepler_equations::*`
//! - `solar::*`
//! - `find_extremas::*`
//!
//! ## References
//!
//! - VSOP87: P. Bretagnon & G. Francou, "Planetary theories in rectangular and spherical variables. VSOP87 solutions", A&A 1988
//! - ELP2000: M. Chapront-Touzé & J. Chapront, "ELP2000-82B: A semi-analytical lunar ephemeris", A&A 1983
//! - Duffett-Smith & Zwart, "Practical Astronomy with your Calculator or Spreadsheet", 4th ed.
//!
//! ## File Structure
//!
//! - `vsop87/`         — VSOP87 planetary theory implementation
//! - `elp2000/`        — ELP2000 lunar theory implementation
//! - `solar/`          — Solar coordinate calculations
//! - `kepler_equations/`— General orbital mechanics
//! - `find_extremas/`  — Extremum-finding algorithms
//!
//! ## See Also
//!
//! - [`coordinates`](../coordinates/index.html): Coordinate systems and transformations
//! - [`units`](../units/index.html): Physical and astronomical units
//! - [`astro`](../astro/index.html): Higher-level astronomical models
//!
//! ---
//! _This module is designed for extensibility and scientific rigor, supporting both
//! high-precision ephemerides and general-purpose celestial mechanics._
//!

pub mod elp2000;
pub mod events;
pub mod kepler_equations;
pub mod pluto;
pub mod root_finding;
pub mod solar;
pub mod vsop87;

pub use pluto::Pluto;
pub use vsop87::VSOP87;
