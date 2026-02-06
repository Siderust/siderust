// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Math Core — Astronomy‑Agnostic Numerical Engine
//!
//! This module contains **all** reusable mathematical logic for scalar event
//! detection, root finding, extrema search, and interval assembly.  It is
//! intentionally astronomy‑agnostic: no coordinate frames, ephemeris libraries,
//! or target types appear here.  Every routine operates on scalar functions of a
//! single `f64` parameter (conceptually *time*) and uses [`qtty`] typed
//! quantities for tolerances and physical arguments where appropriate.
//!
//! ## Sub‑modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`root_finding`] | Brent & bisection solvers for f(t) = 0 |
//! | [`extrema`] | Golden‑section minimiser / maximiser; classify max vs min |
//! | [`intervals`] | Assemble "in‑range" intervals from roots of f(t)−h |
//! | [`bracketing`] | Seed / bracket generation policies (fixed step, adaptive) |

pub mod bracketing;
pub mod extrema;
pub mod intervals;
pub mod root_finding;
