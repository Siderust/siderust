// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Earth-specific gravity providers layered on top of `principia`'s generic kernel.

pub mod egm_low_degree;

pub use egm_low_degree::{LowDegreeEarth, TwoBodyEarth};
pub use principia::gravity::{
    spherical_harmonic_acceleration, GravityConstants, GravityFieldProvider,
};
