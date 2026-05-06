// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Math Core — Astronomy-Agnostic Numerical Engine
//!
//! ## Scientific scope
//!
//! Provides the complete set of scalar numerical primitives needed for
//! astronomical event detection: root finding, extremum search, interval
//! assembly, and bracket generation.  All algorithms are classical numerical
//! analysis methods with well-characterised convergence and error bounds.
//! The module is intentionally astronomy-agnostic — no coordinate frames,
//! ephemeris libraries, or body types appear here.
//!
//! ## Technical scope
//!
//! Every routine operates on a scalar closure `Fn(ModifiedJulianDate) ->
//! Quantity<V>` and uses typed [`qtty`] quantities for tolerances and
//! physical arguments.
//!
//! | Sub-module | Purpose |
//! |------------|---------|
//! | [`root_finding`] | Brent hybrid and bisection solvers for f(t) = 0 |
//! | [`extrema`] | Golden-section minimiser / maximiser; classify max vs min |
//! | [`intervals`] | Assemble "in-range" intervals from roots of f(t) − h |
//! | [`bracketing`] | Seed / bracket generation policies (fixed step, adaptive) |
//!
//! ## References
//!
//! - Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*.
//!   Prentice-Hall.  (Brent hybrid root-solver and golden-section search.)
//! - Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P.
//!   (2007). *Numerical Recipes in C++*, 3rd ed. Cambridge University Press.

pub mod bracketing;
pub mod extrema;
pub mod intervals;
pub mod root_finding;
