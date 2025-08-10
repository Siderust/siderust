//! Numerical kernels underpinning high‑precision astronomy.
//!
//! `calculus` gathers algorithms for orbit propagation and analytic ephemerides
//! used throughout the crate.  The routines aim for reproducibility and are
//! expressed in terms of strongly typed quantities from [`units`](../units/index.html).
//!
//! ## Components
//! - [`vsop87`]: evaluation of the VSOP87 planetary theory.
//! - [`elp2000`]: implementation of the ELP2000‑82B lunar solution.
//! - [`solar`]: solar coordinates and related utilities.
//! - [`kepler_equations`]: Kepler solvers and state‑vector conversions.
//! - [`events`]: searches for transits, culminations and other extrema.
//!
//! ## Example
//!
//! ```rust
//! use siderust::calculus::kepler_equations::solve_keplers_equation;
//! use siderust::units::Radians;
//!
//! let m = Radians::new(1.0);
//! let e = 0.0167;
//! let eccentric_anomaly = solve_keplers_equation(m, e);
//! ```
//!

pub mod solar;
pub mod events;
pub mod vsop87;
pub mod elp2000;

pub mod kepler_equations;

pub use vsop87::VSOP87;

