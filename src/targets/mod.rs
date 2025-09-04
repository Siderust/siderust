//! Astronomical target representation
//!
//! This module defines [`Target`], a lightweight container that couples a position
//! with the time at which the target was seen at that position, together with an optional
//! proper‑motion model.  It is deliberately generic over the coordinate type so
//! it can be reused with Cartesian/Spherical equatorial/ecliptic/icrs.
//!
//! The design goals are:
//! * **Zero‑cost abstraction** – All helper constructors are `const`, allowing
//!   compile‑time evaluation when all arguments are known at build time.
//! * **Flexibility** – The generic parameter `T` lets client code choose any
//!   position type that implements the required operations.
//! * **Clarity** – Specific helpers (`new`, `new_static`, `new_raw`) make the
//!   author’s intent explicit: moving target, fixed target, or advanced manual
//!   construction respectively.
//!
//! ## Examples
//! ```rust
//! use siderust::units::*;
//! use siderust::targets::Target;
//! use siderust::astro::ModifiedJulianDate;
//! use siderust::coordinates::{spherical::Direction, frames::Equatorial, centers::Geocentric};
//! use siderust::astro::proper_motion::ProperMotion;
//!
//! // A star with known proper motion
//! let ra_in_mas_per_year = MilliArcseconds::new(-3.10) / DAY;
//! let dec_in_mas_per_year = MilliArcseconds::new(9.56) / DAY;
//! let betelgeuse_pm = ProperMotion::new(ra_in_mas_per_year, dec_in_mas_per_year);
//! let betelgeuse = Target::new(
//!     Direction::<Geocentric, Equatorial>::new(88.792939*DEG, 7.407064*DEG),
//!     ModifiedJulianDate::new(60200.0).to_julian_day(),
//!     betelgeuse_pm,
//! );
//!
//! // Jupiter’s geocentric position at a given epoch (no proper motion)
//! let jupiter = Target::new_static(
//!     Direction::<Geocentric, Equatorial>::new(23.123*DEG, -5.321*DEG),
//!     ModifiedJulianDate::new(60200.0).to_julian_day(),
//! );
//! ```
//!
//! The optional proper‑motion field allows the same `Target` API to represent
//! both sidereal objects (stars, galaxies) and solar‑system bodies or static
//! catalog entries.

mod target;
mod transform;

pub use target::Target;
