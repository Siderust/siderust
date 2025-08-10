//! Time‑stamped celestial targets.
//!
//! [`Target<T>`] couples a coordinate with the epoch at which it is valid and,
//! optionally, a linear proper‑motion model.  The generic coordinate allows the
//! same structure to describe stars, planets or static catalog entries while
//! remaining allocation‑free.
//!
//! Helper constructors (`new`, `new_static` and `new_raw`) are `const` so they
//! can participate in compile‑time calculations when all inputs are known.
//!
//! ## Example
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

mod transform;
mod target;

pub use target::Target;
