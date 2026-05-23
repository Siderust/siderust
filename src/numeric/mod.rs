// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Numeric
//!
//! ## Scientific scope
//! Reusable numerical primitives for astronomical computations.
//!
//! ## Technical scope
//! Root-finding, extremum detection, interval algebra, and bracketing policies.
//!
//! ## References
//! - Brent, R.P. (1973). *Algorithms for Minimization without Derivatives*. Prentice-Hall.

pub mod bracketing;
pub mod extrema;
pub mod intervals;
pub mod root_finding;
