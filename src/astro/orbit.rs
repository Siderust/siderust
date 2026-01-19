// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Keplerian Orbit Model
//!
//! This module defines the `Orbit` struct, which encapsulates the **six classical Keplerian orbital elements**
//! used to describe the motion of a celestial object around a central body, such as a planet around the Sun.
//!
//! These elements are:
//!
//! 1. **Semi-major axis (`a`)**  
//!    - Defines the size of the orbit.
//!    - It is half the longest diameter of the ellipse.
//!    - Expressed in astronomical units (AstronomicalUnits).
//!
//! 2. **Eccentricity (`e`)**  
//!    - Defines the shape of the orbit.
//!    - Values:
//!       - `e = 0`: circular  
//!       - `0 < e < 1`: elliptical  
//!       - `e = 1`: parabolic  
//!       - `e > 1`: hyperbolic
//!
//! 3. **Inclination (`i`)**  
//!    - The angle between the orbital plane and a reference plane (typically the ecliptic).
//!    - Describes the tilt of the orbit relative to the reference frame.
//!    - Expressed in degrees.
//!
//! 4. **Longitude of the ascending node (`Ω`)**  
//!    - Angle from a fixed reference direction (e.g., the vernal equinox) to the ascending node —
//!      the point where the orbit crosses the reference plane going north.
//!    - Expressed in degrees.
//!
//! 5. **Argument of perihelion (`ω`)**  
//!    - The angle from the ascending node to the perihelion (the point of closest approach).
//!    - Measured in the direction of motion.
//!    - Expressed in degrees.
//!
//! 6. **Mean anomaly at epoch (`M₀`)**  
//!    - Represents the position of the object along its orbit at a specific reference time (epoch).
//!    - It evolves linearly over time and is used to compute the true anomaly.
//!    - Expressed in degrees.
//!
//! ## Epoch
//!
//! The `epoch` is the reference point in time (given in Julian Day) at which the `mean_anomaly_at_epoch` applies.
//! From this point, the object's position can be propagated using Kepler’s equation.
//!
//! ## Coordinate Calculation
//!
//! The `Orbit::heliocentric_coordinates(jd)` method returns the **heliocentric ecliptic Cartesian coordinates**
//! of the orbiting body at a given Julian Day (`jd`), based on the orbital elements and epoch.
//! Internally, this method calls `calculate_orbit_position`, which solves Kepler's equation and performs
//! necessary coordinate transformations.
//!
//! ## Units
//!
//! This module assumes that:
//! - **Angles** are expressed in degrees (`Degrees`).
//! - **Distances** use astronomical units (`AstronomicalUnits`).
//! - **Time** is expressed as Julian Days (`JulianDate`).
//!
//! ## Usage Example
//!
//! This example computes Earth's position on a given Julian date.
//!
//! ```rust
//! use siderust::astro::orbit::Orbit;
//! use siderust::astro::JulianDate;
//! use qtty::*;
//!
//! let earth_orbit = Orbit::new(
//!     1.0*AU,                    // a
//!     0.0167,                    // e
//!     Degrees::new(0.00005),     // i
//!     Degrees::new(-11.26064),   // Ω
//!     Degrees::new(102.94719),   // ω
//!     Degrees::new(100.46435),   // M₀
//!     JulianDate::J2000,         // epoch (J2000)
//! );
//!
//! let coords = earth_orbit.kepler_position(JulianDate::new(2459200.5));
//! ```

use crate::astro::JulianDate;
use qtty::*;

/// Represents the Keplerian orbital elements of a celestial object.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    pub semi_major_axis: AstronomicalUnits, // Semi-major axis (AstronomicalUnits)
    pub eccentricity: f64,                  // Orbital eccentricity
    pub inclination: Degrees,               // Inclination (degrees)
    pub longitude_of_ascending_node: Degrees, // Longitude of ascending node (Ω)
    pub argument_of_perihelion: Degrees,    // Argument of perihelion (ω)
    pub mean_anomaly_at_epoch: Degrees,     // Mean anomaly at epoch (M₀)
    pub epoch: JulianDate,                  // Epoch (Julian Dat
}

impl Orbit {
    /// Creates a new set of orbital elements.
    pub const fn new(
        semi_major_axis: AstronomicalUnits,
        eccentricity: f64,
        inclination: Degrees,
        longitude_of_ascending_node: Degrees,
        argument_of_perihelion: Degrees,
        mean_anomaly_at_epoch: Degrees,
        epoch: JulianDate,
    ) -> Self {
        Self {
            semi_major_axis,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_perihelion,
            mean_anomaly_at_epoch,
            epoch,
        }
    }
}
