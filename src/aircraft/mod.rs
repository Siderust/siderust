// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Aircraft — tracked atmospheric bodies
//!
//! ## Scientific scope
//!
//! Aircraft are physical objects whose state (geodetic position, inertial
//! velocity, magnetic heading, barometric altitude) evolves within Earth's
//! atmosphere. Their position must be composable with the rest of the
//! `siderust` coordinate model — in particular with the same typed
//! [`Geodetic`](crate::coordinates::centers::Geodetic) positions used by
//! observatories, and with the [`Trackable`](crate::targets::Trackable)
//! abstraction that lets observation-planning code treat aircraft and
//! astronomical bodies uniformly.
//!
//! This module does **not** model:
//! - flight planning, routing, or airspace queries (belongs in `constops`);
//! - aircraft performance or aerodynamics;
//! - ATC orchestration or replay (belongs in `constops`).
//!
//! ## Technical scope
//!
//! - [`Aircraft`] — static identity record: ICAO 24-bit address, callsign,
//!   and optional wake turbulence category.
//! - [`AircraftState`] — snapshot of position ([`Geodetic<ECEF>`]) plus
//!   geodetic altitude, horizontal velocity (ground-speed, track angle), and
//!   vertical rate, all typed via [`crate::qtty`].
//! - [`isa`] — ICAO International Standard Atmosphere model: pressure ↔
//!   geometric altitude, layer-by-layer temperature lapse, geopotential to
//!   geometric altitude conversion.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::aircraft::{Aircraft, AircraftState};
//! use siderust::coordinates::centers::Geodetic;
//! use siderust::coordinates::frames::ECEF;
//! use siderust::qtty::{Degrees, Meters};
//!
//! let ac = Aircraft::new(0x4CA2B5, "EIN104");
//!
//! let state = AircraftState::new(
//!     Geodetic::<ECEF>::new(
//!         Degrees::new(-6.270),
//!         Degrees::new(53.421),
//!         Meters::new(10_668.0),
//!     ),
//!     Degrees::new(275.0),
//! );
//!
//! assert_eq!(ac.callsign(), "EIN104");
//! let alt_m = state.barometric_altitude_m();
//! assert!(alt_m.value() > 0.0);
//! ```
//!
//! ## References
//!
//! - ICAO (2010). *Annex 2 — Rules of the Air*, 10th edition. International
//!   Civil Aviation Organization.
//! - ICAO (1993). *Manual of the ICAO Standard Atmosphere*, 3rd edition
//!   (Doc 7488). International Civil Aviation Organization.
//! - Blythe, D. (2002). "Decoding Mode S transponder signals." *The
//!   Aeronautical Journal* 106 (1062), pp. 423–432.

mod inner;
pub mod isa;

pub use inner::{Aircraft, AircraftState, AircraftTrack, MetersPerSecond, WakeCategory};
