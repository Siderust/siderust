// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orbit-state covariance matrix re-exported from `principia`.
//!
//! [`StateCovariance`] is a 6×6 symmetric positive-semi-definite matrix
//! carrying position-velocity uncertainty in a given coordinate frame.
//! The type is defined and implemented in `principia`; this module
//! re-exports it so callers of `siderust::astro::dynamics` do not need
//! a direct `principia` dependency.

pub use principia::StateCovariance;
