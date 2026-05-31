// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # VSOP87 Planetary Theory
//!
//! ## Scientific scope
//!
//! VSOP87 (*Variations Séculaires des Orbites Planétaires* 1987) is a
//! semi-analytical planetary theory by Bretagnon and Francou (Bureau des
//! Longitudes, Paris).  It models the heliocentric and barycentric positions
//! of the eight major planets as trigonometric power series in time:
//!
//! ```text
//!   X = Σ_k  T^k · Σ_j  A_{kj} · cos(B_{kj} + C_{kj} · T)
//! ```
//!
//! where `T` is Julian millennia from J2000.0 (TDB),
//! coefficients `(A, B, C)` are tabulated in the generated data files,
//! and `X ∈ {X, Y, Z}` in astronomical units.
//!
//! Two variants are used here:
//!
//! - **VSOP87A** — heliocentric ecliptic rectangular coordinates, equinox J2000.0.
//! - **VSOP87E** — barycentric ecliptic rectangular coordinates, equinox J2000.0.
//!
//! Accuracy: better than 1″ over ±2000 years; better than 1 AU·mas over the
//! 1900–2100 interval for inner planets.
//!
//! ## Technical scope
//!
//! - [`VSOP87`] — trait exposing `vsop87a` and `vsop87e` on each planet type;
//!   both return a `Position<…, EclipticMeanJ2000, AstronomicalUnit>`.
//! - `vsop87a` / `vsop87e` — generated coefficient tables re-exported from
//!   `src/generated/vsop87{a,e}.rs`; compiled into `vsop_data` sub-modules.
//! - `vsop87_impl` — SIMD-accelerated series evaluation kernels (`position`,
//!   `velocity`, `position_velocity`), all operating on raw `f64` (AU and AU/day).
//!
//! ## References
//!
//! - Bretagnon, P., & Francou, G. (1988). "Planetary theories in rectangular
//!   and spherical variables: VSOP87 solutions".
//!   *Astronomy and Astrophysics* 202, 309–315.
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. Willmann-Bell.
//! - IMCCE VSOP87 coefficient archive:
//!   <https://www.imcce.fr/inpop/ephemerides/vsop87/>

mod vsop87_impl;

#[allow(clippy::all)]
mod vsop87a;
#[allow(clippy::all)]
mod vsop87e;

mod vsop87_trait;

pub use vsop87_trait::VSOP87;

use vsop87_impl::position;
use vsop87_impl::position_velocity;
use vsop87_impl::velocity;
