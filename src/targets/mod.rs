// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astronomical target representation
//!
//! ## Scientific scope
//!
//! An *astronomical target* is the coupling of a position on the sky with the
//! epoch at which that position is asserted, together with an optional
//! proper-motion model that propagates the position to other epochs. This is
//! the standard formal model used by stellar catalogues (e.g. Hipparcos,
//! Gaia) and solar-system ephemerides: a position is only meaningful when
//! tagged with the time it refers to and, for non-stationary objects, with a
//! rule for how to move it forward or backward.
//!
//! The default coordinate frame in `siderust` for catalogue-style targets is
//! ICRS / Equatorial Mean J2000.0; proper-motion rates follow the IAU 2006
//! formalism `(μ_α* = μ_α · cos δ, μ_δ)` so that great-circle motion is
//! preserved. The model is valid for stars and small-body solar-system
//! objects whose proper motion is well determined; it does *not* model
//! parallax, light-time, or barycentric corrections — those live in
//! [`crate::astro`].
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - [`CoordinateWithPM<T>`] — a stamped coordinate snapshot pairing any
//!   typed coordinate `T` with a [`crate::time::JulianDate`] and an optional
//!   proper-motion model. Constructors are `const` where possible.
//! - [`Target<T>`] — backward-compatible alias for `CoordinateWithPM<T>`.
//! - [`Trackable`] — trait abstracting "anything that can produce a position
//!   at time *t*"; implemented for solar-system unit types, `Star`,
//!   `direction::ICRS`, and `CoordinateWithPM<T>` itself (identity).
//!
//! Time inputs are typed [`crate::time::JulianDate`] / [`crate::time::ModifiedJulianDate`];
//! coordinates retain their declared frame and unit at the type level.
//! Propagation kernels live in [`crate::astro::proper_motion`]; this module
//! only defines the data containers.
//!
//! ## Examples
//!
//! ```rust
//! use siderust::qtty::*;
//! use siderust::targets::CoordinateWithPM;
//! use siderust::time::ModifiedJulianDate;
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use siderust::coordinates::spherical::Direction;
//! use siderust::astro::proper_motion::ProperMotion;
//!
//! type MasPerYear = siderust::qtty::Per<siderust::qtty::MilliArcsecond, siderust::qtty::Year>;
//! type MasPerYearQ = siderust::qtty::Quantity<MasPerYear>;
//! let betelgeuse_pm = ProperMotion::from_mu_alpha_star::<MasPerYear>(
//!     MasPerYearQ::new(27.54),
//!     MasPerYearQ::new(10.86),
//! );
//! let betelgeuse = CoordinateWithPM::new(
//!     Direction::<EquatorialMeanJ2000>::new(88.792939*DEG, 7.407064*DEG),
//!     ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60200.0)).into(),
//!     betelgeuse_pm,
//! );
//!
//! let jupiter = CoordinateWithPM::new_static(
//!     Direction::<EquatorialMeanJ2000>::new(23.123*DEG, -5.321*DEG),
//!     ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60200.0)).into(),
//! );
//! ```
//!
//! ## References
//!
//! - International Astronomical Union (2006). Resolution B1.5 ("Definition
//!   of the Barycentric Celestial Reference System") and the IAU 2006
//!   proper-motion formalism. IAU Transactions XXVIB.
//! - Kovalevsky, J., & Seidelmann, P. K. (2004). *Fundamentals of
//!   Astrometry*. Cambridge University Press. ISBN 978-0-521-64216-7.
//! - ESA (1997). *The Hipparcos and Tycho Catalogues*. ESA SP-1200,
//!   §1.2 (proper-motion conventions).

mod target;
mod trackable;
mod transform;

pub use target::{CoordinateWithPM, Target};
pub use trackable::Trackable;
