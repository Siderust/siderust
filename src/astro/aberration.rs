//! # Annual Aberration (Ron–Vondrák 2005)
//!
//! This module converts **mean (geometric) coordinates** of a
//! celestial object into **apparent coordinates** by accounting for
//! **annual aberration** – the apparent displacement of the source caused by
//! the velocity of the Earth with respect to the solar-system barycentre.
//!
//! ```text
//! maximum effect ≃ 20.5″  (on the ecliptic, at right angles to the apex)
//! typical precision   <   0.1 mas (1900–2100) using this implementation
//! ```
//!
//! ## What is aberration?
//! When an observer moves at velocity **v**, the direction of an incoming light
//! ray is aberrated by an angle \( \Delta\theta \approx v/c \) (for
//! non-relativistic speeds).  For the Earth, the dominant term is the
//! **annual** component (orbital speed ≃ 29.79 km s⁻¹).  High-precision star
//! catalogues list *mean* positions; to compare with a real observation we
//! must add aberration.
//!
//! ## Theory used
//! We implement the Ron–Vondrák trigonometric series published in
//! *IERS Conventions 2003* (TN 32, chap. 5, table 5.1) and kept unchanged in
//! the 2010/2020 editions.  The theory expresses the heliocentric velocity of
//! the Earth **v** as a sum of **36 terms**:
//!
//! ```text
//! v_x = Σ (S₁ₙ + S₂ₙ·T)·sin Aₙ + (C₁ₙ + C₂ₙ·T)·cos Aₙ  (likewise y, z)
//! ```
//! where *Aₙ* = Σ kᵢ·Φᵢ is a linear combination of **11 fundamental arguments**
//! (l₂ … F).  Coefficients are tabulated in **10⁻⁸ au d⁻¹**; the constant
//! _C_ below is the speed of light in the *same* units, so **v/c** is a
//! dimensionless ratio ~10⁻⁴.
//!
//! ## References
//! * IERS Technical Note 32 (2003), §5, table 5.1  
//! * IERS Technical Note 36 (2010) – unchanged  
//! * Kaplan, G.H. (2005) *USNO Circular 179*, eqs. (2.10)–(2.13)  
//! * *Astronomical Almanac* (2024), Sect. C, eqs. 22.1–22.4
//!
//! ## Implementation notes
//! - Fully vectorial formulation avoids singularities at ±90° declination
//!   (renormalisation error ~10⁻¹²).
//! - Fundamental angles are reduced modulo 2π with `rem_euclid(τ)` to maintain
//!   precision for epochs far from J2000.0.
//! - Large static coefficient tables are generated automatically from the
//!   IERS ASCII source.

use crate::units::JulianDay;
use crate::coordinates::{
    //spherical::Position, cartesian::Position,
    spherical, cartesian,
    centers::Geocentric, frames::Equatorial
};


//--------------------------------------------------------------------
// Physical & mathematical constants
//--------------------------------------------------------------------

/// Speed of light **c = 17314 463 350 × 10⁻⁸ au d⁻¹** (IAU 2012 resolution B2).
///
/// The factor 10⁻⁸ matches the published Ron–Vondrák coefficients.
pub const C_10E8: f64 = 17_314_463_350.0;

/// 2 π, useful for argument reduction.
pub const TAU: f64 = std::f64::consts::TAU;


const TERMS: usize = 36;


/// Integer multipliers of the fundamental arguments (l₂ … F).
#[derive(Copy, Clone, Debug)]
struct Arg {
    pub a_l2: i32,
    pub a_l3: i32,
    pub a_l4: i32,
    pub a_l5: i32,
    pub a_l6: i32,
    pub a_l7: i32,
    pub a_l8: i32,
    pub a_ll: i32,
    pub a_d:  i32,
    pub a_mm: i32,
    pub a_f:  i32,
}


/// Trigonometric coefficients (10⁻⁸ au d⁻¹).
#[derive(Copy, Clone, Debug)]
struct Xyz {
    pub sin1: i32,
    pub sin2: i32,
    pub cos1: i32,
    pub cos2: i32,
}

fn heliocentric_velocity_components(t: f64) -> cartesian::Position<Geocentric, Equatorial> {
    let mut vx = 0.0;
    let mut vy = 0.0;
    let mut vz = 0.0;

    let l2 = (3.176_146_7 + 1_021.328_554_6 * t).rem_euclid(TAU);
    let l3 = (1.753_470_3 +   628.307_584_9 * t).rem_euclid(TAU);
    let l4 = (6.203_480_9 +   334.061_243_1 * t).rem_euclid(TAU);
    let l5 = (0.599_546_4 +    52.969_096_5 * t).rem_euclid(TAU);
    let l6 = (0.874_016_8 +    21.329_909_095 * t).rem_euclid(TAU);
    let l7 = (5.481_293_9 +     7.478_159_9 * t).rem_euclid(TAU);
    let l8 = (5.311_886_3 +     3.813_303_6 * t).rem_euclid(TAU);
    let ll = (3.810_344_4 + 8_399.684_733_7 * t).rem_euclid(TAU);
    let  d = (5.198_466_7 + 7_771.377_148_6 * t).rem_euclid(TAU);
    let mm = (2.355_555_9 + 8_328.691_428_9 * t).rem_euclid(TAU);
    let  f = (1.627_905_2 + 8_433.466_160_1 * t).rem_euclid(TAU);

    for i in 0..TERMS {
        let arg =
              ARGUMENTS[i].a_l2 as f64 * l2
            + ARGUMENTS[i].a_l3 as f64 * l3
            + ARGUMENTS[i].a_l4 as f64 * l4
            + ARGUMENTS[i].a_l5 as f64 * l5
            + ARGUMENTS[i].a_l6 as f64 * l6
            + ARGUMENTS[i].a_l7 as f64 * l7
            + ARGUMENTS[i].a_l8 as f64 * l8
            + ARGUMENTS[i].a_ll as f64 * ll
            + ARGUMENTS[i].a_d  as f64 *  d
            + ARGUMENTS[i].a_mm as f64 * mm
            + ARGUMENTS[i].a_f  as f64 *  f;

        let (s, c) = arg.sin_cos();

        vx += (X_COEFFICIENTS[i].sin1 as f64 + X_COEFFICIENTS[i].sin2 as f64 * t) * s
            + (X_COEFFICIENTS[i].cos1 as f64 + X_COEFFICIENTS[i].cos2 as f64 * t) * c;
        vy += (Y_COEFFICIENTS[i].sin1 as f64 + Y_COEFFICIENTS[i].sin2 as f64 * t) * s
            + (Y_COEFFICIENTS[i].cos1 as f64 + Y_COEFFICIENTS[i].cos2 as f64 * t) * c;
        vz += (Z_COEFFICIENTS[i].sin1 as f64 + Z_COEFFICIENTS[i].sin2 as f64 * t) * s
            + (Z_COEFFICIENTS[i].cos1 as f64 + Z_COEFFICIENTS[i].cos2 as f64 * t) * c;
    }
    cartesian::Position::new(vx, vy, vz)
}

/// Add **annual aberration** to a mean equatorial Cartesian position.
///
/// * `mean` – Geocentric Cartesian coordinates referred to the true equator &
///   equinox of date (in astronomical units).
/// * `jd`   – Terrestrial Time (TT) in Julian Day.
///
/// # Returns
/// A new [`cartesian::Position`] whose x, y, z components include the effect of
/// the Earth's orbital velocity.
///
/// # Accuracy
/// Better than 0.1 mas for dates 1900-2100; dominated by the underlying
/// Ron–Vondrák theory.
#[must_use]
pub fn apply_aberration(
    mean: cartesian::Position<Geocentric, Equatorial>,
    jd:   JulianDay,
) -> cartesian::Position<Geocentric, Equatorial> {

    if mean.distance_from_origin() == 0.0 {
        // Don't look at your feet!
        return mean;
    }

    let t = jd.julian_centuries().value();

    //--------------------------------------------------------------------
    // 1. Heliocentric velocity components  (10⁻⁸ au d⁻¹)
    //--------------------------------------------------------------------
    let velocity = heliocentric_velocity_components(t);

    //--------------------------------------------------------------------
    // 2. Apply v/c to the unit vector of the star
    //--------------------------------------------------------------------
    // Normalize to obtain the unit vector
    let norm: cartesian::Position<Geocentric, Equatorial> = mean.normalize();
    let with_vc = norm + velocity / C_10E8;
    with_vc.normalize() * mean.distance_from_origin()
}

#[must_use]
pub fn remove_aberration(
    app: cartesian::Position<Geocentric, Equatorial>,
    jd:  JulianDay,
) -> cartesian::Position<Geocentric, Equatorial> {

    if app.distance_from_origin() == 0.0 {
        // Don't look at your feet!
        return app;
    }

    let t  = jd.julian_centuries().value();

    //--------------------------------------------------------------------
    // 1. Heliocentric velocity components  (10⁻⁸ au d⁻¹)
    //--------------------------------------------------------------------
    let velocity = heliocentric_velocity_components(t);

    //--------------------------------------------------------------------
    // 2. Apply -v/c to the unit vector of the star
    //--------------------------------------------------------------------
    let norm: cartesian::Position<Geocentric, Equatorial> = app.normalize();
    let with_vc = norm - velocity / C_10E8;
    with_vc.normalize() * app.distance_from_origin()
}


#[must_use]
pub fn apply_aberration_sph(
    mean: spherical::Position<Geocentric, Equatorial>,
    jd:   JulianDay,
) -> spherical::Position<Geocentric, Equatorial> {
    (&apply_aberration((&mean).into(), jd)).into()
}

const ARGUMENTS: [Arg; TERMS] = [
    Arg { a_l2: 0, a_l3: 1, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 2, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 1, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 1, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 3, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 1, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 1 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 1, a_d: 0, a_mm: 1, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 2, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 2, a_l4: 0, a_l5: -1, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 3, a_l4: -8, a_l5: 3, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 5, a_l4: -8, a_l5: 3, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 2, a_l3: -1, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 1, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 1, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 1, a_l4: 0, a_l5: -2, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 1, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 1, a_l4: 0, a_l5: 1, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 2, a_l3: -2, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 1, a_l4: 0, a_l5: -1, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 4, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 3, a_l4: 0, a_l5: -2, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 1, a_l3: -2, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 2, a_l3: -3, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 2, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 2, a_l3: 4, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 3, a_l4: -2, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 1, a_d: 2, a_mm: -1, a_f: 0 },
    Arg { a_l2: 8, a_l3: 12, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 8, a_l3: 14, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 2, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 3, a_l3: 4, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 2, a_l4: 0, a_l5: -2, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 3, a_l3: -3, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 2, a_l4: -2, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 0, a_d: 0, a_mm: 0, a_f: 0 },
    Arg { a_l2: 0, a_l3: 0, a_l4: 0, a_l5: 0, a_l6: 0, a_l7: 0, a_l8: 0, a_ll: 1, a_d: -2, a_mm: 0, a_f: 0 }
];

const X_COEFFICIENTS: [Xyz; TERMS] = [
    Xyz { sin1: -1719914, sin2: -2, cos1: -25, cos2: 0 },
    Xyz { sin1: 6434, sin2: 141, cos1: 28007, cos2: -107 },
    Xyz { sin1: 715, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 715, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 486, sin2: -5, cos1: -236, cos2: -4 },
    Xyz { sin1: 159, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 39, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 33, sin2: 0, cos1: -10, cos2: 0 },
    Xyz { sin1: 31, sin2: 0, cos1: 1, cos2: 0 },
    Xyz { sin1: 8, sin2: 0, cos1: -28, cos2: 0 },
    Xyz { sin1: 8, sin2: 0, cos1: -28, cos2: 0 },
    Xyz { sin1: 21, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -19, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 17, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 16, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 16, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 11, sin2: 0, cos1: -1, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -11, cos2: 0 },
    Xyz { sin1: -11, sin2: 0, cos1: -2, cos2: 0 },
    Xyz { sin1: -7, sin2: 0, cos1: -8, cos2: 0 },
    Xyz { sin1: -10, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -9, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -9, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -9, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -9, cos2: 0 },
    Xyz { sin1: 8, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 8, sin2: 0, cos1: 0, cos2: 0 }, 
    Xyz { sin1: -4, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: -4, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: -6, sin2: 0, cos1: -5, cos2: 0 },
    Xyz { sin1: -1, sin2: 0, cos1: -1, cos2: 0 },
    Xyz { sin1: 4, sin2: 0, cos1: -6, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: 5, sin2: 0, cos1: -5, cos2: 0 },
    Xyz { sin1: 5, sin2: 0, cos1: 0, cos2: 0 }
];

const Y_COEFFICIENTS: [Xyz; TERMS] = [
    Xyz { sin1: 25, sin2: -13, cos1: 1578089, cos2: 156 },
    Xyz { sin1: 25697, sin2: -95, cos1: -5904, cos2: -130 },
    Xyz { sin1: 6, sin2: 0, cos1: -657, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -656, cos2: 0 },
    Xyz { sin1: -216, sin2: -4, cos1: -446, cos2: 5 },
    Xyz { sin1: 2, sin2: 0, cos1: -147, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 26, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -36, cos2: 0 },
    Xyz { sin1: -9, sin2: 0, cos1: -30, cos2: 0 },
    Xyz { sin1: 1, sin2: 0, cos1: -28, cos2: 0 },
    Xyz { sin1: 25, sin2: 0, cos1: 8, cos2: 0 },
    Xyz { sin1: -25, sin2: 0, cos1: -8, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -19, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 17, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -16, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 15, cos2: 0 },
    Xyz { sin1: 1, sin2: 0, cos1: -15, cos2: 0 },
    Xyz { sin1: -1, sin2: 0, cos1: -10, cos2: 0 },
    Xyz { sin1: -10, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -2, sin2: 0, cos1: 9, cos2: 0 },
    Xyz { sin1: -8, sin2: 0, cos1: 6, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 9, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -9, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -8, cos2: 0 },
    Xyz { sin1: -8, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 8, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -8, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: -6, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: 6, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: -4, sin2: 0, cos1: 5, cos2: 0 },
    Xyz { sin1: -2, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: -5, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: -6, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -4, sin2: 0, cos1: -5, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -5, cos2: 0 }
];

const Z_COEFFICIENTS: [Xyz; TERMS] = [
    Xyz { sin1: 10, sin2: 32, cos1: 684185, cos2: -358 },
    Xyz { sin1: 11141, sin2: -48, cos1: -2559, cos2: -55 },
    Xyz { sin1: -15, sin2: 0, cos1: -282, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -285, cos2: 0 },
    Xyz { sin1: -94, sin2: 0, cos1: -193, cos2: 0 },
    Xyz { sin1: -6, sin2: 0, cos1: -61, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 59, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 16, cos2: 0 },
    Xyz { sin1: -5, sin2: 0, cos1: -13, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -12, cos2: 0 },
    Xyz { sin1: 11, sin2: 0, cos1: 3, cos2: 0 },
    Xyz { sin1: -11, sin2: 0, cos1: -3, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -8, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 8, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -7, cos2: 0 },
    Xyz { sin1: 1, sin2: 0, cos1: 7, cos2: 0 },
    Xyz { sin1: -3, sin2: 0, cos1: -6, cos2: 0 },
    Xyz { sin1: -1, sin2: 0, cos1: 5, cos2: 0 },
    Xyz { sin1: -4, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -1, sin2: 0, cos1: 4, cos2: 0 },
    Xyz { sin1: -3, sin2: 0, cos1: 3, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: 4, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: -3, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 3, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -3, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -3, cos2: 0 },
    Xyz { sin1: -3, sin2: 0, cos1: 2, cos2: 0 },
    Xyz { sin1: 3, sin2: 0, cos1: -2, cos2: 0 },
    Xyz { sin1: -2, sin2: 0, cos1: 2, cos2: 0 },
    Xyz { sin1: 1, sin2: 0, cos1: -4, cos2: 0 },
    Xyz { sin1: -2, sin2: 0, cos1: -2, cos2: 0 },
    Xyz { sin1: -3, sin2: 0, cos1: 0, cos2: 0 },
    Xyz { sin1: -2, sin2: 0, cos1: -2, cos2: 0 },
    Xyz { sin1: 0, sin2: 0, cos1: -2, cos2: 0 }
];



#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    use approx::assert_relative_eq;

    #[test]
    fn test_aberration_preserva_distance_and_epoch() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = spherical::Position::<Geocentric, Equatorial>::new(
            Degrees::new(10.0),
            Degrees::new(20.0),
            1.23
        );
        let out = apply_aberration_sph(mean, jd);

        assert_relative_eq!(out.distance, mean.distance, epsilon = 0.0);
    }

    #[test]
    fn test_aberration_introduces_shift() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = spherical::Position::<Geocentric, Equatorial>::new(
            Degrees::new(0.0),    // RA = 0°
            Degrees::new(0.0),    // Dec = 0°
            1.0
        );
        let out = apply_aberration_sph(mean, jd);

        let delta_ra = out.ra().diff_deg(mean.ra());
        let delta_dec = out.dec().diff_deg(mean.dec());
        assert!(delta_ra.as_f64() > 0.0 || delta_dec.as_f64() > 0.0,
            "Expected a change in RA or Dec");
        assert!(delta_ra.as_f64() < 0.01 && delta_dec.as_f64() < 0.01,
            "Shift is too large")
    }

    #[test]
    fn test_aberration_at_north_pole() {
        let jd = JulianDay::new(2451545.0);
        let mean = spherical::Position::<Geocentric, Equatorial>::new(
            Degrees::new(123.4),  // dummy RA
            Degrees::new(90.0),   // Dec = +90°
            1.0
        );
        let out = apply_aberration_sph(mean, jd);

        assert!(out.dec().as_f64() < 90.0, "Declination should decrease slightly at pole");
        assert!(!out.ra().as_f64().is_nan(), "RA must not be NaN at the pole");
    }

}
