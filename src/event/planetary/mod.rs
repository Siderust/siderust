// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Planetary physical-orientation geometry
//!
//! ## Scientific scope
//!
//! Light-time geometric sky-plane orientation of major planets: rotational-pole
//! and bright-limb position angles, plus Saturn ring openings. Observers may be
//! geocentric or topocentric; sky bearings may use true-of-date, mean J2000, or
//! ICRS equatorial frames. See [`orientation`] for the correction policy.
//!
//! ## Technical scope
//!
//! Generic entry point [`body_orientation`] over body, ephemeris, center, and
//! frame. Convenience wrappers preserve Mars/Saturn geocentric true-of-date
//! defaults. Saturn rings are layered on the shared geometry via
//! [`saturn_rings`] / [`saturn_ring_orientation`].
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Chapters 42 and 45.
//!   Willmann–Bell.
//! - Archinal, B. A. et al., IAU WGCCRE rotational-element reports.

pub mod orientation;

pub use orientation::{
    body_orientation, mars_orientation, saturn_ring_orientation, saturn_rings, BodyOrientation,
    MarsOrientation, PlanetaryObservationCenter, SaturnRingOrientation, SkyOrientationFrame,
};
