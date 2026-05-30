// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Pluto ephemeris computation using the Meeus/Williams abbreviated series.
//!
//! This module provides a single public function [`Pluto::get_heliocentric`] that
//! returns Pluto's heliocentric ecliptic position for a given Julian Day.
//!
//! The algorithm is taken from Jean Meeus, *Astronomical Algorithms*, 2nd ed.
//! (1998), chapter 36, which in turn reproduces the numerical series derived
//! by T. G. Williams (1991). Accuracy is on the order of **0.5 arc‑seconds in
//! longitude** and **0.2 arc‑seconds in latitude** over the interval 1885 – 2099,
//! while being dramatically faster than a full numerical integration or the
//! DE ephemerides.
//!
//! Coefficient tables live in `siderust_archive::pluto::pluto_data`.
//!
//! ```text
//! Step outline (see code for details):
//!  1. Convert the supplied Julian Day to Julian centuries T w.r.t. J2000.0.
//!  2. Compute mean longitudes λJ, λS, λP of Jupiter, Saturn and Pluto using
//!     linear expressions in T.
//!  3. For each of the 42 / 43 periodic terms:
//!       A. Form the argument  Ai = j·λJ + s·λS + p·λP  (degrees).
//!       B. Evaluate sin Ai and cos Ai.
//!       C. Accumulate contributions to ΣL, ΣB, ΣR via  A·sin Ai + B·cos Ai.
//!  4. Scale the sums and add constant/base terms to obtain L, B, R.
//!  5. Convert spherical ⟨L,B,R⟩ to Cartesian ⟨x,y,z⟩ in the ecliptic frame.
//! ```
//!
//! # References
//! * Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann‑Bell.
//! * Williams, T. G. (1991). "An optimized algorithm for Pluto", *Mem. Brit.
//!   Astron. Assoc.* **99** (2), 75–82.

use crate::coordinates::{cartesian, centers::Heliocentric, frames::EclipticMeanJ2000, spherical};
use crate::qtty::{AstronomicalUnit, Degrees, Radian, AU};
use crate::time::JulianDate;
use siderust_archive::pluto::pluto_data::{
    PLUTO_ARGUMENTS, PLUTO_LATITUDE_TERMS, PLUTO_LONGITUDE_TERMS, PLUTO_RADIUS_TERMS,
};

/// Marker struct for Pluto ephemeris computations via the Meeus/Williams series.
pub struct Pluto;

impl Pluto {
    /// Compute Pluto's heliocentric ecliptic rectangular position for the
    /// given Julian Day (TT time scale).
    pub fn get_heliocentric(
        jd: JulianDate,
    ) -> cartesian::Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        let t = jd.julian_centuries();
        let jupiter_lon = Degrees::new(34.35 + 3034.9057 * t);
        let saturn_lon = Degrees::new(50.08 + 1222.1138 * t);
        let pluto_lon = Degrees::new(238.96 + 144.9600 * t);

        // 3. Initialize sums for the periodic terms.
        let mut sum_longitude = Degrees::new(0.0);
        let mut sum_latitude = Degrees::new(0.0);
        let mut sum_radius = 0.0;

        // 4. Loop over all periodic terms.
        for i in 0..PLUTO_ARGUMENTS.len() {
            // Calculate the argument:
            let a = jupiter_lon * PLUTO_ARGUMENTS[i].j
                + saturn_lon * PLUTO_ARGUMENTS[i].s
                + pluto_lon * PLUTO_ARGUMENTS[i].p;

            // Convert 'a' from degrees to radians.
            let (sin_a, cos_a) = a.to::<Radian>().sin_cos();

            // Add periodic corrections for longitude, latitude, and radius.
            sum_longitude +=
                Degrees::new(PLUTO_LONGITUDE_TERMS[i].a * sin_a + PLUTO_LONGITUDE_TERMS[i].b * cos_a);
            sum_latitude += Degrees::new(PLUTO_LATITUDE_TERMS[i].a * sin_a + PLUTO_LATITUDE_TERMS[i].b * cos_a);
            sum_radius += PLUTO_RADIUS_TERMS[i].a * sin_a + PLUTO_RADIUS_TERMS[i].b * cos_a;
        }

        // 5. Calculate the final heliocentric spherical coordinates.
        // These base values and scale factors come from the chosen model (e.g., Meeus's data).
        let lon = sum_longitude * 0.000001 + Degrees::new(238.958116 + 144.96 * t);
        let lat = sum_latitude * 0.000001 - Degrees::new(3.908239);
        let rad = sum_radius * 0.0000001 + 40.7241346;

        // EclipticMeanJ2000: lon = azimuth, lat = polar
        spherical::Position::<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
            lon,
            lat,
            rad * AU,
        )
        .to_cartesian()
    }
}

#[cfg(test)]
mod tests {
    use crate::ephemeris::pluto::Pluto;
    use crate::qtty::AstronomicalUnits;

    #[test]
    fn pluto_heliocentric_position_j2000() {
        let pos = Pluto::get_heliocentric(crate::J2000);
        assert!(
            (pos.x() - AstronomicalUnits::new(-9.875333629852145))
                .abs()
                .value()
                < 1e-6
        );
        assert!(
            (pos.y() - AstronomicalUnits::new(-27.958786187190157))
                .abs()
                .value()
                < 1e-6
        );
        assert!(
            (pos.z() - AstronomicalUnits::new(5.850444258527083))
                .abs()
                .value()
                < 1e-6
        );
    }
}
