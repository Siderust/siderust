// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Center-shift implementations for **position** types.
//!
//! ## Scientific scope
//!
//! A center shift (origin translation) moves the reference point from which
//! positions are measured. Key shifts in astrodynamics include:
//!
//! - **Geocentric → Topocentric**: subtract the observer's geocentric position
//!   vector (the *topocentric parallax* correction, ≈ 6 400 km maximum).
//! - **Heliocentric → Bodycentric**: shift from a heliocentric origin to the
//!   centre of a solar-system body.
//!
//! ## Technical scope
//!
//! Sub-modules implement the [`TransformCenter`](crate::coordinates::transform::centers::TransformCenter) trait for specific center pairs:
//!
//! - [`to_topocentric`]: Geocentric → Topocentric using site geodetic
//!   coordinates and Earth-orientation data.
//! - [`to_bodycentric`]: Heliocentric → bodycentric (for planetary positions).
//!
//! Center shifts only apply to **positions** (affine points); directions and
//! velocities are translation-invariant and are not handled here.
//!
//! ## References
//!
//! - Urban, S. E. & Seidelmann, P. K. (2013). *Explanatory Supplement to the
//!   Astronomical Almanac*, 3rd ed. §2.3 (parallax corrections).

pub mod to_bodycentric;
pub mod to_topocentric;

pub use to_topocentric::{to_topocentric_with, to_topocentric_with_ctx};
