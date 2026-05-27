// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Mission geometry
//!
//! ## Scientific scope
//!
//! Library-level mission-geometry kernels for observer-to-target line of
//! sight, body obstruction, occultation overlap, eclipse classification,
//! and local orbital timing. These are the typed primitives that
//! mission-analysis code and operational scheduling systems can build
//! on; this crate intentionally stops at the math layer. Mission-control
//! execution, scheduling, REST services, CLI/GUI front ends, and
//! operator-console workflows belong in **SatOps**.
//!
//! ## Technical scope
//!
//! All public surfaces accept typed inputs (typed centers/frames in
//! [`crate::coordinates`], typed quantities in [`qtty`]) and return
//! typed results:
//!
//! - [`AzElRange`] — observer-relative azimuth, elevation, slant range,
//!   range-rate, one-way light time, with explicit obstruction / mask
//!   diagnostics.
//! - [`azimuth_elevation_range`] — solve [`AzElRange`] from two state
//!   vectors expressed in the same Cartesian center / frame.
//! - [`line_of_sight_obstructed`] — strict spherical obstruction test for
//!   an occulting body sitting between observer and target.
//! - [`occultation_fraction`] — apparent disk overlap (0…1) between two
//!   bodies as seen from the observer.
//! - [`EclipseState`] / [`eclipse_state`] / [`solar_eclipsing`] — full /
//!   partial / no-eclipse classification.
//! - [`beta_angle`], [`local_solar_time`], [`ltan`] — orbit-relative
//!   geometry helpers.
//! - [`Fov`] / [`Instrument`] — field-of-view geometry and instrument
//!   metadata.
//! - [`TerrainMask`] — piece-wise linear terrain-mask elevation profile.
//! - [`LocalFrame`] — local east-north-up reference frame at an
//!   observer site.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §3.4, §3.7, §5.3, §11.3.
//! - Montenbruck, O., & Gill, E. (2000). *Satellite Orbits — Models,
//!   Methods, and Applications*, §3.4, §11.2.

// ─────────────────────────────────────────────────────────────────────────────
// Sub-modules
// ─────────────────────────────────────────────────────────────────────────────

pub mod az_el_range;
pub mod eclipse;
pub mod fov;
pub mod local_frame;
pub mod occultation;
pub mod orbit_relative;
pub mod terrain_mask;

// ─────────────────────────────────────────────────────────────────────────────
// Flat re-exports for ergonomic access via `mission::geometry::*`
// ─────────────────────────────────────────────────────────────────────────────

pub use az_el_range::{azimuth_elevation_range, AzElRange, LineOfSightStatus, MetersPerSecond};
pub use eclipse::{eclipse_state, solar_eclipsing, EclipseState};
pub use fov::{Fov, Instrument};
pub use local_frame::LocalFrame;
pub use occultation::{line_of_sight_obstructed, occultation_fraction};
pub use orbit_relative::{beta_angle, local_solar_time, ltan};
pub use terrain_mask::TerrainMask;
