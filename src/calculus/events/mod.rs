//! # Events Module
//!
//! This module provides generic tools for detecting and analyzing astronomical events
//! such as extrema (maxima, minima), culminations, and other significant points in the
//! evolution of time-dependent celestial quantities.
//!
//! ## What Is an "Event" in Celestial Mechanics?
//!
//! In astronomy, an "event" typically refers to a moment in time when a function of
//! interest (e.g., altitude, distance, brightness) reaches a local maximum, minimum,
//! or crosses a threshold. Examples include:
//! - The culmination (highest point) of a celestial body in the sky.
//! - Closest approach (pericenter) or farthest distance (apocenter) in an orbit.
//! - Times of rise, set, or transit for planets and stars.
//! - Extremes in observable properties (e.g., brightness, angular separation).
//!
//! Detecting such events is fundamental for ephemeris generation, observation planning,
//! and scientific analysis.
//!
//! ## Where & Why Is It Used?
//!
//! * Planetarium and ephemeris software to list daily or nightly events.
//! * Automated telescope scheduling and observation planning.
//! * Scientific studies of orbital dynamics and light curves.
//! * Educational tools to visualize celestial phenomena in real time.
//!
//! ## Mathematical & Algorithmic Background
//!
//! The core problem is to find the time(s) `t` where a function `f(t)` has:
//! - A local maximum or minimum (extremum).
//! - A root (crosses a specific value).
//!
//! This module provides robust, efficient algorithms for:
//! - **Static extrema search:** Find extrema in a fixed interval for a known function.
//! - **Dynamic extrema search:** Track extrema as parameters or conditions change over time.
//!
//! ### Numerical Strategies
//!
//! - **Bracketing and root-finding:** Reliable detection of roots and extrema using
//!   bracketing, bisection, and interpolation.
//! - **Sampling and refinement:** Adaptive sampling to locate candidate events, followed
//!   by local refinement for high precision.
//! - **Generic interfaces:** Accepts user-supplied functions or closures, enabling
//!   application to a wide range of astronomical problems.
//!
//! ## Public API
//!
//! ```rust
//! use siderust::calculus::events::{find_static_extremas, find_dynamic_extremas, Culmination};
//! use siderust::bodies::catalog::SIRIUS;
//! use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
//! use qtty::*;
//! use siderust::astro::JulianDate;
//!
//! // Example: Find extrema of a function in a given interval
//! let extrema = find_static_extremas(
//!     &SIRIUS.target,                    // Target to be observed
//!     &ROQUE_DE_LOS_MUCHACHOS,           // Location of the observer
//!     JulianDate::J2000,                  // start time
//!     JulianDate::J2000 + Days::new(1.0), // end time
//! );
//!
//! for ext in extrema {
//!     match ext {
//!         Culmination::Upper { jd } => println!("Maximum culmination at t = {}", jd),
//!         Culmination::Lower { jd } => println!("Minimum culmination at t = {}", jd),
//!     }
//! }
//! ```
//!
//! ## Module Layout
//!
//! | Item                        | Visibility | Purpose                                 |
//! |-----------------------------|------------|-----------------------------------------|
//! | `Culmination`               | **pub**    | Enum for maxima/minima event reporting  |
//! | `find_static_extremas`      | **pub**    | Find extrema in a static interval       |
//! | `find_dynamic_extremas`     | **pub**    | Track extrema as conditions evolve      |
//!
//! ## Limitations & Future Work
//!
//! * Current algorithms are generic but may require tuning for highly oscillatory or
//!   discontinuous functions.
//! * No built-in support yet for event types beyond extrema (e.g., threshold crossings).
//! * PRs adding support for additional event types, SIMD acceleration, or higher-precision
//!   backends are welcome.
//!
//! This API is designed to be stable and extensible for future Siderust releases.

/// Represents a culmination event — the moment a celestial body crosses
/// the observer’s meridian.
///
/// - `Upper`: transit across the upper meridian (highest altitude).  
/// - `Lower`: transit across the lower meridian (lowest altitude).
///
/// The `jd` field stores the Julian Day of the event.
#[derive(Debug)]
pub enum Culmination {
    Upper { jd: crate::astro::JulianDate },
    Lower { jd: crate::astro::JulianDate },
}

pub mod altitude_periods;
mod find_dynamic_extremas;
mod find_static_extremas;

pub use altitude_periods::{
    find_altitude_periods, find_sun_above_altitude, find_sun_below_altitude,
    find_sun_in_altitude_range, sun_altitude_rad, twilight, AltitudeCondition, AltitudePeriod,
};
pub use find_dynamic_extremas::find_dynamic_extremas;
pub use find_static_extremas::find_static_extremas;
