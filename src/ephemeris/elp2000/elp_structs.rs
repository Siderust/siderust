// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # ELP2000-82B Coefficient Structures
//!
//! ## Scientific scope
//!
//! The ELP2000-82B theory (Chapront-Touzé & Chapront 1988) partitions its
//! Poisson-series terms into three families depending on which gravitational
//! perturbations they model:
//!
//! - **Main Problem** — the Earth–Moon–Sun three-body problem without planetary
//!   perturbations.  Arguments are the four Delaunay angles D, M, Mʹ, F.
//! - **Earth Perturbations** — figure, tidal, and relativistic corrections
//!   involving the same four Delaunay angles plus an explicit power of T.
//! - **Planetary Perturbations** — gravitational influence of the eight planets
//!   (Mercury through Neptune), using eleven combined planetary/lunar arguments.
//!
//! Each family is represented by a distinct coefficient record type defined here.
//!
//! ## Technical scope
//!
//! - [`MainProblem`] — one sine term of the Main Problem series; amplitude `a`
//!   in 10⁻⁴ arc-seconds, integer Delaunay multipliers, and phase polynomial
//!   coefficients b₀…b₅.
//! - [`EarthPert`] — one Earth-perturbation term; amplitude `a` in 10⁻⁴
//!   arc-seconds, T-exponent `iz` (dimensionless), phase rate `o`
//!   (radians / Julian-millennia), and constant phase `p` (radians).
//! - [`PlanetPert`] — one planetary-perturbation term; amplitude `theta` in
//!   10⁻⁴ arc-seconds, eleven integer planetary multipliers, phase rate `o`
//!   (radians / Julian-millennia), and constant phase `p` (radians).
//!
//! All numeric fields use raw `f64`; the typed conversion to physical units
//! (arc-seconds → radians → km) is performed in [`super::elp_series`].
//!
//! ## References
//!
//! - Chapront-Touzé, M., & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical
//!   lunar ephemeris adequate for historical times".
//!   *Astronomy and Astrophysics* 190, 342–352.

#![allow(unreachable_pub)]

pub use crate::archive::elp::{EarthPert, MainProblem, PlanetPert};
