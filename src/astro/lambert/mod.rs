// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lambert solver — typed 0/N-revolution Lambert boundary-value solver.
//!
//! ## Overview
//!
//! Solves Lambert's two-point boundary-value problem (find the conic transfer
//! connecting two position vectors in a prescribed time of flight) using
//! Izzo's 2014 reformulation with third-order Householder iteration. Both the
//! single-revolution branch (`N = 0`) and the multi-revolution branches
//! (`N ≥ 1`, two solutions per `N`) are implemented.
//!
//! ## Public API
//!
//! - [`lambert`] / [`lambert_n_rev`] — **typed** entry-points operating on
//!   [`affn::cartesian::Position`] / [`affn::cartesian::Velocity`] and
//!   [`qtty::Second`]. These are the documented public interface.
//! - [`solve_lambert`] / [`solve_lambert_n_rev`] — low-level numeric backend
//!   on plain `[f64; 3]` km / km/s arrays. Useful for FFI and for inner loops
//!   in mission-design search where allocating typed wrappers per evaluation
//!   is not desired.
//! - [`LambertError`] — single, unified error enum, re-exported here.
//!
//! ## Examples
//!
//! Earth-centred prograde transfer (Vallado *Fundamentals of Astrodynamics*
//! Ex. 7-5, 4th ed., p. 467):
//!
//! ```
//! use affn::cartesian::Position;
//! use affn::centers::ReferenceCenter;
//! use affn::frames::ICRS;
//! use qtty::dynamics::GravitationalParameter;
//! use qtty::length::Kilometer;
//! use qtty::Second;
//! use siderust::astro::lambert::{lambert, LambertBranch};
//!
//! let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
//! let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
//! let tof = Second::new(4_560.0);
//! let mu = GravitationalParameter::new(398_600.4418);
//!
//! let sol = lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
//! assert!((sol.v1.x().value() - 2.058913).abs() < 1e-3);
//! ```
//!
//! ## Out of scope
//!
//! * Mission-design optimisation, porkchop search, transfer-window
//!   enumeration → future `siderust-mission-design`.
//! * Frame / center bookkeeping at the algorithm boundary beyond what
//!   the typed entry-points provide.
//!
//! ## References
//!
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*. Celest. Mech. Dyn.
//!   Astron., 121(1):1–15.
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and
//!   Methods of Astrodynamics* (Rev. ed.). AIAA.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications* (4th ed.). Microcosm Press.

#![forbid(unsafe_code)]

mod error;
mod izzo;
mod typed;

pub use error::LambertError;
pub use izzo::{
    solve_lambert, solve_lambert_n_rev, LambertBranch, LambertDiagnostics, LambertSolution,
    NRevBranch,
};
pub use typed::{lambert, lambert_n_rev, TypedLambertSolution};
