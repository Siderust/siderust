// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Planetary physical-orientation geometry
//!
//! ## Scientific scope
//!
//! Geocentric sky-plane orientation of Mars and Saturn for physical
//! observations: rotational-pole (and Saturn ring-axis) position angles,
//! Mars bright-limb bearing, and Saturn ring opening angles as seen from
//! Earth and from the Sun. Epochs are TT; planetary positions are retarded
//! for light-time. See [`orientation`] for the full correction policy.
//!
//! ## Technical scope
//!
//! Free functions [`mars_orientation`] and [`saturn_ring_orientation`] select
//! the ephemeris backend at compile time via [`PlanetaryEphemeris`].
//! Convenience wrappers live on [`Mars`](crate::bodies::solar_system::Mars) and
//! [`Saturn`](crate::bodies::solar_system::Saturn).
//!
//! ## References
//!
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Chapters 42 and 45.
//!   Willmann–Bell.
//! - Archinal, B. A. et al., IAU WGCCRE rotational-element reports.

pub mod orientation;

pub use orientation::{
    mars_orientation, saturn_ring_orientation, MarsOrientation, SaturnRingOrientation,
};
