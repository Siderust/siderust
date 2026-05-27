// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Mission-analysis building blocks
//!
//! ## Scientific scope
//!
//! Typed primitives for mission analysis: instrument field-of-view
//! geometry, terrain masking, observer-to-target line-of-sight solutions,
//! eclipse classification, orbit-relative geometry, and the runtime
//! context that aggregates these objects for operational use. These are
//! the math and metadata layers; mission-control execution, scheduling,
//! REST services, and operator-console workflows belong in **SatOps**.
//!
//! ## Technical scope
//!
//! - [`geometry`]: Line-of-sight azimuth/elevation/range, FoV
//!   containment, terrain masks, eclipse state, and orbit-relative
//!   geometry (beta angle, LST, LTAN).
//! - [`context`]: [`context::MissionContext`] — runtime mission-analysis
//!   context aggregating instrument and site catalogs.
//! - [`site`]: [`site::Location`] — observation-site / ground-station
//!   metadata.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §1.2, §3.4, §3.7.
//! - Wertz, J. R. (2011). *Mission Geometry; Orbit and Constellation
//!   Design and Management*. §10.

#![forbid(unsafe_code)]

pub mod context;
pub mod geometry;
pub mod site;
