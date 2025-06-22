//! Pluto ephemeris computation using the Meeus/Williams abbreviated series.
//!
//! This module provides a single public function [`Pluto::get_heliocentric`] that
//! returns Pluto's heliocentric ecliptic position for a given Julian Day.
//!
//! The algorithm is taken from Jean Meeus, *Astronomical Algorithms*, 2nd ed.
//! (1998), chapter 36, which in turn reproduces the numerical series derived
//! by T. G. Williams (1991). Accuracy is on the order of **0.5 arc‑seconds in
//! longitude** and **0.2 arc‑seconds in latitude** over the interval 1885 – 2099,
//! while being dramatically faster than a full numerical integration or the
//! DE ephemerides.
//!
//! ```text
//! Step outline (see code for details):
//!  1. Convert the supplied Julian Day to Julian centuries T w.r.t. J2000.0.
//!  2. Compute mean longitudes λJ, λS, λP of Jupiter, Saturn and Pluto using
//!     linear expressions in T.
//!  3. For each of the 42 / 43 periodic terms:
//!       A. Form the argument  Ai = j·λJ + s·λS + p·λP  (degrees).
//!       B. Evaluate sin Ai and cos Ai.
//!       C. Accumulate contributions to ΣL, ΣB, ΣR via  A·sin Ai + B·cos Ai.
//!  4. Scale the sums and add constant/base terms to obtain L, B, R.
//!  5. Convert spherical ⟨L,B,R⟩ to Cartesian ⟨x,y,z⟩ in the ecliptic frame.
//! ```
//!
//! # References
//! * Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann‑Bell.
//! * Williams, T. G. (1991). "An optimized algorithm for Pluto", *Mem. Brit.
//!   Astron. Assoc.* **99** (2), 75–82.


use crate::units::{Degrees, JulianDay};
use crate::coordinates::{
    Vector, SphericalCoord,
    centers::Heliocentric,
    frames::Ecliptic
};

pub struct Pluto;

impl Pluto {

    pub fn get_heliocentric(jd: JulianDay) -> Vector<Heliocentric, Ecliptic> {
        let t = jd.julian_centuries().value();

        // 2. Calculate mean longitudes (in degrees) for Jupiter, Saturn, and Pluto.
        let jupiter_lon =  Degrees::new(34.35 + 3034.9057 * t);
        let saturn_lon  =  Degrees::new(50.08 + 1222.1138 * t);
        let pluto_lon   =  Degrees::new(238.96 + 144.9600 * t);

        // 3. Initialize sums for the periodic terms.
        let mut sum_longitude = Degrees::new(0.0);
        let mut sum_latitude  = Degrees::new(0.0);
        let mut sum_radius    = 0.0;
    
        // 4. Loop over all periodic terms.
        for i in 0..ARGUMENTS.len() {
            // Calculate the argument:
            let a = jupiter_lon * ARGUMENTS[i].j +
                    saturn_lon  * ARGUMENTS[i].s +
                    pluto_lon   * ARGUMENTS[i].p;
            
            // Convert 'a' from degrees to radians.
            let (sin_a, cos_a) = a.to_radians().sin_cos();

            // Add periodic corrections for longitude, latitude, and radius.
            sum_longitude += Degrees::new(LONGITUDE_TERMS[i].a * sin_a + LONGITUDE_TERMS[i].b * cos_a);
            sum_latitude  += Degrees::new(LATITUDE_TERMS[i].a  * sin_a + LATITUDE_TERMS[i].b  * cos_a);
            sum_radius    += RADIUS_TERMS[i].a    * sin_a + RADIUS_TERMS[i].b    * cos_a;
        }

        // 5. Calculate the final heliocentric spherical coordinates.
        // These base values and scale factors come from the chosen model (e.g., Meeus's data).
        SphericalCoord::<Heliocentric, Ecliptic>::new (
            sum_longitude * 0.000001 + Degrees::new(238.958116 + 144.96 * t),
            sum_latitude  * 0.000001 - Degrees::new(3.908239),
            sum_radius    * 0.0000001 + 40.7241346,
            JulianDay::new(t)
        ).to_cartesian()
    }
}

/// A single row of multipliers for the argument (the jupiter_lon, saturn_lon, pluto_lon coefficients).
#[derive(Debug, Clone, Copy)]
struct PlutoArgument {
    pub j: f64,
    pub s: f64,
    pub p: f64,
}

/// A single row of periodic coefficients (the A and B amplitudes).
#[derive(Debug, Clone, Copy)]
struct PlutoTerm {
    pub a: f64,
    pub b: f64,
}

const ARGUMENTS: &[PlutoArgument] = &[
    PlutoArgument { j: 0.0, s:  0.0, p:  1.0 },
    PlutoArgument { j: 0.0, s:  0.0, p:  2.0 },
    PlutoArgument { j: 0.0, s:  0.0, p:  3.0 },
    PlutoArgument { j: 0.0, s:  0.0, p:  4.0 },
    PlutoArgument { j: 0.0, s:  0.0, p:  5.0 },
    PlutoArgument { j: 0.0, s:  0.0, p:  6.0 },
    PlutoArgument { j: 0.0, s:  1.0, p: -1.0 },
    PlutoArgument { j: 0.0, s:  1.0, p:  0.0 },
    PlutoArgument { j: 0.0, s:  1.0, p:  1.0 },
    PlutoArgument { j: 0.0, s:  1.0, p:  2.0 },
    PlutoArgument { j: 0.0, s:  1.0, p:  3.0 },
    PlutoArgument { j: 0.0, s:  2.0, p: -2.0 },
    PlutoArgument { j: 0.0, s:  2.0, p: -1.0 },
    PlutoArgument { j: 0.0, s:  2.0, p:  0.0 },
    PlutoArgument { j: 1.0, s: -1.0, p:  0.0 },
    PlutoArgument { j: 1.0, s: -1.0, p:  1.0 },
    PlutoArgument { j: 1.0, s:  0.0, p: -3.0 },
    PlutoArgument { j: 1.0, s:  0.0, p: -2.0 },
    PlutoArgument { j: 1.0, s:  0.0, p: -1.0 },
    PlutoArgument { j: 1.0, s:  0.0, p:  0.0 },
    PlutoArgument { j: 1.0, s:  0.0, p:  1.0 },
    PlutoArgument { j: 1.0, s:  0.0, p:  2.0 },
    PlutoArgument { j: 1.0, s:  0.0, p:  3.0 },
    PlutoArgument { j: 1.0, s:  0.0, p:  4.0 },
    PlutoArgument { j: 1.0, s:  1.0, p: -3.0 },
    PlutoArgument { j: 1.0, s:  1.0, p: -2.0 },
    PlutoArgument { j: 1.0, s:  1.0, p: -1.0 },
    PlutoArgument { j: 1.0, s:  1.0, p:  0.0 },
    PlutoArgument { j: 1.0, s:  1.0, p:  1.0 },
    PlutoArgument { j: 1.0, s:  1.0, p:  3.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -6.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -5.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -4.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -3.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -2.0 },
    PlutoArgument { j: 2.0, s:  0.0, p: -1.0 },
    PlutoArgument { j: 2.0, s:  0.0, p:  0.0 },
    PlutoArgument { j: 2.0, s:  0.0, p:  1.0 },
    PlutoArgument { j: 2.0, s:  0.0, p:  2.0 },
    PlutoArgument { j: 2.0, s:  0.0, p:  3.0 },
    PlutoArgument { j: 3.0, s:  0.0, p: -2.0 },
    PlutoArgument { j: 3.0, s:  0.0, p: -1.0 },
    PlutoArgument { j: 3.0, s:  0.0, p:  0.0 }
];

const LONGITUDE_TERMS: &[PlutoTerm] = &[
    PlutoTerm { a: -19799805.0, b: 19850055.0},
    PlutoTerm { a:  897144.0,   b: -4954829.0},
    PlutoTerm { a:  611149.0,   b: 1211027.0},
    PlutoTerm { a: -341243.0,   b: -189585.0},
    PlutoTerm { a:  129287.0,   b: -34992.0},
    PlutoTerm { a: -38164.0,    b: 30893.0},
    PlutoTerm { a:  20442.0,    b: -9987.0},
    PlutoTerm { a: -4063.0,     b: -5071.0},
    PlutoTerm { a: -6016.0,     b: -3336.0},
    PlutoTerm { a: -3956.0,     b: 3039.0},
    PlutoTerm { a: -667.0,      b: 3572.0},
    PlutoTerm { a:  1276.0,     b: 501.0},
    PlutoTerm { a:  1152.0,     b: -917.0},
    PlutoTerm { a:  630.0,      b: -1277.0},
    PlutoTerm { a:  2571.0,     b: -459.0},
    PlutoTerm { a:  899.0,      b: -1449.0},
    PlutoTerm { a: -1016.0,     b: 1043.0},
    PlutoTerm { a: -2343.0,     b: -1012.0},
    PlutoTerm { a:  7042.0,     b: 788.0},
    PlutoTerm { a:  1199.0,     b: -338.0},
    PlutoTerm { a:  418.0,      b: -67.0},
    PlutoTerm { a:  120.0,      b: -274.0},
    PlutoTerm { a: -60.0,       b: -159.0},
    PlutoTerm { a: -82.0,       b: -29.0},
    PlutoTerm { a: -36.0,       b: -20.0},
    PlutoTerm { a: -40.0,       b: 7.0},
    PlutoTerm { a: -14.0,       b: 22.0},
    PlutoTerm { a:  4.0,        b: 13.0},
    PlutoTerm { a:  5.0,        b: 2.0},
    PlutoTerm { a: -1.0,        b: 0.0},
    PlutoTerm { a:  2.0,        b: 0.0},
    PlutoTerm { a: -4.0,        b: 5.0},
    PlutoTerm { a:  4.0,        b: -7.0},
    PlutoTerm { a:  14.0,       b: 24.0},
    PlutoTerm { a: -49.0,       b: -34.0},
    PlutoTerm { a:  163.0,      b: -48.0},
    PlutoTerm { a:  9.0,        b: 24.0},
    PlutoTerm { a: -4.0,        b: 1.0},
    PlutoTerm { a: -3.0,        b: 1.0},
    PlutoTerm { a:  1.0,        b: 3.0},
    PlutoTerm { a: -3.0,        b: -1.0},
    PlutoTerm { a:  5.0,        b: -3.0},
    PlutoTerm { a:  0.0,        b: 0.0}
];

const LATITUDE_TERMS: &[PlutoTerm] = &[
    PlutoTerm { a: -5452852.0, b: -14974862.0},
    PlutoTerm { a:  3527812.0, b:  1672790.0},
    PlutoTerm { a: -1050748.0, b:  327647.0},
    PlutoTerm { a:  178690.0,  b: -292153.0},
    PlutoTerm { a:  18650.0,   b:  100340.0},
    PlutoTerm { a: -30697.0,   b: -25823.0},
    PlutoTerm { a:  4878.0,    b:  11248.0},
    PlutoTerm { a:  226.0,     b: -64.0},
    PlutoTerm { a:  2030.0,    b: -836.0},
    PlutoTerm { a:  69.0,      b: -604.0},
    PlutoTerm { a: -247.0,     b: -567.0},
    PlutoTerm { a: -57.0,      b:  1.0},
    PlutoTerm { a: -122.0,     b:  175.0},
    PlutoTerm { a: -49.0,      b: -164.0},
    PlutoTerm { a: -197.0,     b:  199.0},
    PlutoTerm { a: -25.0,      b:  217.0},
    PlutoTerm { a:  589.0,     b: -248.0},
    PlutoTerm { a: -269.0,     b:  711.0},
    PlutoTerm { a:  185.0,     b:  193.0},
    PlutoTerm { a:  315.0,     b:  807.0},
    PlutoTerm { a: -130.0,     b: -43.0},
    PlutoTerm { a:  5.0,       b:  3.0},
    PlutoTerm { a:  2.0,       b:  17.0},
    PlutoTerm { a:  2.0,       b:  5.0},
    PlutoTerm { a:  2.0,       b:  3.0},
    PlutoTerm { a:  3.0,       b:  1.0},
    PlutoTerm { a:  2.0,       b: -1.0},
    PlutoTerm { a:  1.0,       b: -1.0},
    PlutoTerm { a:  0.0,       b: -1.0},
    PlutoTerm { a:  0.0,       b:  0.0},
    PlutoTerm { a:  0.0,       b: -2.0},
    PlutoTerm { a:  2.0,       b:  2.0},
    PlutoTerm { a: -7.0,       b:  0.0},
    PlutoTerm { a:  10.0,      b: -8.0},
    PlutoTerm { a: -3.0,       b:  20.0},
    PlutoTerm { a:  6.0,       b:  5.0},
    PlutoTerm { a:  14.0,      b:  17.0},
    PlutoTerm { a: -2.0,       b:  0.0},
    PlutoTerm { a:  0.0,       b:  0.0},
    PlutoTerm { a:  0.0,       b:  0.0},
    PlutoTerm { a:  0.0,       b:  1.0},
    PlutoTerm { a:  0.0,       b:  0.0},
    PlutoTerm { a:  1.0,       b:  0.0}
];

const RADIUS_TERMS: &[PlutoTerm] = &[
    PlutoTerm { a:  66865439.0, b:  68951812.0},
    PlutoTerm { a: -11827535.0, b: -332538.0},
    PlutoTerm { a:  1593179.0,  b: -1438890.0},
    PlutoTerm { a: -18444.0,    b:  483220.0},
    PlutoTerm { a: -65977.0,    b: -85431.0},
    PlutoTerm { a:  31174.0,    b: -6032.0},
    PlutoTerm { a: -5794.0,     b:  22161.0},
    PlutoTerm { a:  4601.0,     b:  4032.0},
    PlutoTerm { a: -1729.0,     b:  234.0},
    PlutoTerm { a: -415.0,      b:  702.0},
    PlutoTerm { a:  239.0,      b:  723.0},
    PlutoTerm { a:  67.0,       b: -67.0},
    PlutoTerm { a:  1034.0,     b: -451.0},
    PlutoTerm { a: -129.0,      b:  504.0},
    PlutoTerm { a:  480.0,      b: -231.0},
    PlutoTerm { a:  2.0,        b: -441.0},
    PlutoTerm { a: -3359.0,     b:  265.0},
    PlutoTerm { a:  7856.0,     b: -7832.0},
    PlutoTerm { a:  36.0,       b:  45763.0},
    PlutoTerm { a:  8663.0,     b:  8547.0},
    PlutoTerm { a: -809.0,      b: -769.0},
    PlutoTerm { a:  263.0,      b: -144.0},
    PlutoTerm { a: -126.0,      b:  32.0},
    PlutoTerm { a: -35.0,       b: -16.0},
    PlutoTerm { a: -19.0,       b: -4.0},
    PlutoTerm { a: -15.0,       b:  8.0},
    PlutoTerm { a: -4.0,        b:  12.0},
    PlutoTerm { a:  5.0,        b:  6.0},
    PlutoTerm { a:  3.0,        b:  1.0},
    PlutoTerm { a:  6.0,        b: -2.0},
    PlutoTerm { a:  2.0,        b:  2.0},
    PlutoTerm { a: -2.0,        b: -2.0},
    PlutoTerm { a:  14.0,       b:  13.0},
    PlutoTerm { a: -63.0,       b:  13.0},
    PlutoTerm { a:  136.0,      b: -236.0},
    PlutoTerm { a:  273.0,      b:  1065.0},
    PlutoTerm { a:  251.0,      b:  149.0},
    PlutoTerm { a: -25.0,       b: -9.0},
    PlutoTerm { a:  9.0,        b: -2.0},
    PlutoTerm { a: -8.0,        b:  7.0},
    PlutoTerm { a:  2.0,        b: -10.0},
    PlutoTerm { a:  19.0,       b:  35.0},
    PlutoTerm { a:  10.0,       b:  2.0}
];
