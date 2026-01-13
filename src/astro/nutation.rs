//! # Earth Nutation Module
//!
//! **Nutation** describes short‑period wobbles of the Earth’s rotation axis that
//! ride on top of the much slower 26 000‑year precession.  The motion is driven
//! by the changing torque exerted by the Moon and the Sun on the Earth’s
//! equatorial bulge. Ignoring nutation would introduce errors of up to
//! ±17″ of arc in the position of a star, far larger than the field of view of
//! modern telescopes.
//!
//! The module supplies:
//! * **Δψ (longitude)**: the shift of the ecliptic along its own plane (°).
//! * **Δε (obliquity)**: the oscillation of the ecliptic’s tilt (°).
//! * **ε₀ (mean obliquity)**: the mean tilt of the ecliptic at the same epoch (°).
//!
//! With those three numbers you can rotate *mean* equatorial coordinates
//! (*RA*, *Dec*) of any object into *true* (apparent) coordinates valid for the
//! requested date.
//!
//! ## Numerical model
//! We implement the **IAU 1980** nutation theory (63 trigonometric terms). It is
//! still accurate to ≲ 0.1″ from year 1800 to 2050 and matches the recipe in
//! Chapter 22 of _Jean Meeus – Astronomical Algorithms_, 2nd ed. (1998).
//!
//! ```text
//! Δψ = Σ (A₁ + A₂·T) · sin(arg)   (in 0.0001″)
//! Δε = Σ (B₁ + B₂·T) · cos(arg)   (in 0.0001″)
//! arg = D·D + M·M + M′·M′ + F·F + Ω·Ω   (all in radians)
//! T   = (JDE − J2000) / 36525   (Julian centuries)
//! ```
//!
//! The fundamental arguments *D, M, M′, F, Ω* follow IERS 2003 expressions and
//! are evaluated in radians for numerical stability.
//!
//! ## API surface
//! * [`get_nutation`] → [`Nutation`].  Computes Δψ, Δε, ε₀ for a given
//!   [`JulianDate`].  All outputs are **degrees**.
//! * [`corrected_ra_with_nutation`]  
//!   Input: a mean [`Position`] (RA/Dec) and the same date.  Output: the
//!   apparent right ascension *αₜ* (°) after a 3‑1‑3 rotation using Δψ, Δε, ε₀.
//!
//! ## Quick example
//! ```rust
//! use chrono::prelude::*;
//! use siderust::astro::JulianDate;
//! use siderust::bodies::catalog::SIRIUS;
//! use siderust::astro::nutation::{get_nutation, corrected_ra_with_nutation};
//! use siderust::astro::precession::precess_from_j2000;
//!
//! let jd = JulianDate::from_utc(Utc::now());
//! let n = get_nutation(jd);
//! println!("Δψ = {:.4}°, Δε = {:.4}°", n.longitude, n.obliquity);
//!
//! let mean_of_date = precess_from_j2000(SIRIUS.target.get_position().clone(), jd);
//! let ra_app = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
//! println!("Apparent RA = {ra_app:.4}°");
//! ```
//!
//! ## Limitations
//! *Accuracy*: ≤ 0.1″ for 1800–2050; outside that span consider the full IAU 2000A
//! series (1365 terms) or IERS tabulated Δψ/Δε values.  

use crate::astro::dynamical_time::julian_ephemeris_day;
use crate::astro::JulianDate;
use crate::coordinates::spherical::direction::EquatorialMeanOfDate;
use affn::Rotation3;
use qtty::*;

/// Nutation components for a given epoch (all **degrees**).
#[derive(Debug)]
pub struct Nutation {
    /// Δψ: nutation in ecliptic longitude.
    pub longitude: Degrees,
    /// Δε: nutation in obliquity.
    pub obliquity: Degrees,
    /// ε₀: mean obliquity of the ecliptic.
    pub ecliptic: Degrees,
}

#[derive(Copy, Clone)]
struct NutationArguments {
    d: f64,
    m: f64,
    mm: f64,
    f: f64,
    o: f64,
}
#[derive(Copy, Clone)]
struct NutationCoefficients {
    longitude1: f64,
    longitude2: f64,
    obliquity1: f64,
    obliquity2: f64,
}

const TERMS: usize = 63;

/// Compute Δψ, Δε and ε₀ for the supplied Julian Day (JD).
#[inline]
pub fn get_nutation(jd: JulianDate) -> Nutation {
    let jde = julian_ephemeris_day(jd);
    let t = jde.julian_centuries().value();
    let t2 = t * t;
    let t3 = t2 * t;

    // Fundamental arguments (radians)
    let d = Degrees::new(297.850_36 + 445_267.111_480 * t - 0.001_914_2 * t2 + t3 / 189_474.0)
        .to::<Radian>();
    let m = Degrees::new(357.527_72 + 35_999.050_340 * t - 0.000_160_3 * t2 - t3 / 300_000.0)
        .to::<Radian>();
    let mp = Degrees::new(134.962_98 + 477_198.867_398 * t + 0.008_697_2 * t2 + t3 / 56_250.0)
        .to::<Radian>();
    let f = Degrees::new(93.271_91 + 483_202.017_538 * t - 0.003_682_5 * t2 + t3 / 327_270.0)
        .to::<Radian>();
    let om = Degrees::new(125.044_52 - 1_934.136_261 * t + 0.002_070_8 * t2 + t3 / 450_000.0)
        .to::<Radian>();

    // Evaluate trigonometric series (0.0001″ units)
    let mut dpsi = 0.0;
    let mut deps = 0.0;
    for i in 0..TERMS {
        let arg: Radians = d * ARGUMENTS[i].d
            + m * ARGUMENTS[i].m
            + mp * ARGUMENTS[i].mm
            + f * ARGUMENTS[i].f
            + om * ARGUMENTS[i].o;
        let a = COEFFICIENTS[i].longitude1 + COEFFICIENTS[i].longitude2 * t;
        let b = COEFFICIENTS[i].obliquity1 + COEFFICIENTS[i].obliquity2 * t;
        dpsi += a * arg.sin();
        deps += b * arg.cos();
    }

    // convert 0.0001″ → degrees
    dpsi /= 10000.0 * 3600.0;
    deps /= 10000.0 * 3600.0;

    // Mean obliquity ε₀ (″ → °)
    let eps0 = 23.0 + 26.0 / 60.0 + 21.448 / 3600.0 - 46.8150 / 3600.0 * t - 0.00059 / 3600.0 * t2
        + 0.001813 / 3600.0 * t3;

    Nutation {
        longitude: Degrees::new(dpsi),
        obliquity: Degrees::new(deps),
        ecliptic: Degrees::new(eps0),
    }
}

/// Rotate a mean position (RA, Dec) into **apparent** right ascension, applying nutation.
#[inline]
pub fn corrected_ra_with_nutation(target: &EquatorialMeanOfDate, jd: JulianDate) -> Degrees {
    let rot = nutation_rotation(jd);
    let cart = target.to_cartesian();
    let [x, y, _z] = rot.apply_array([cart.x(), cart.y(), cart.z()]);
    Degrees::new(y.atan2(x).to_degrees().rem_euclid(360.0))
}

/// Nutation rotation matrix from mean-of-date to true-of-date.
///
/// The `JulianDate` input is interpreted as TT.
pub fn nutation_rotation(jd: JulianDate) -> Rotation3 {
    let Nutation {
        longitude,
        obliquity,
        ecliptic,
    } = get_nutation(jd);

    let dpsi = longitude.to::<Radian>().value();
    let deps = obliquity.to::<Radian>().value();
    let eps0 = ecliptic.to::<Radian>().value();

    // R1(ε0+Δε) · R3(Δψ) · R1(−ε0)
    Rotation3::from_x_rotation(eps0 + deps)
        * Rotation3::from_z_rotation(dpsi)
        * Rotation3::from_x_rotation(-eps0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::spherical;
    use qtty::{Degrees, Radian};

    #[test]
    fn nutation_rotation_matches_ra_correction() {
        let jd = JulianDate::new(2_459_000.5);
        let mean = spherical::direction::EquatorialMeanOfDate::new(
            Degrees::new(123.0),
            Degrees::new(-10.0),
        );

        let rot = nutation_rotation(jd);
        let cart = mean.to_cartesian();
        let [x, y, z] = rot.apply_array([cart.x(), cart.y(), cart.z()]);

        let ra = y.atan2(x).to_degrees().rem_euclid(360.0);
        let ra = Degrees::new(ra);
        let ra_corr = corrected_ra_with_nutation(&mean, jd);

        let diff = ra.abs_separation(ra_corr).to::<Radian>().value();
        assert!(diff < 1e-10, "RA mismatch: {}", diff);

        let _ = z; // preserve unused warning if additional checks are added later
    }
}

const ARGUMENTS: [NutationArguments; TERMS] = [
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 1.0,
        mm: 0.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 1.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: -2.0,
        m: -1.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: -1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: -1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: -1.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 2.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: -2.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 2.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 2.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: -1.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 2.0,
        mm: 0.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: -1.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 2.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 1.0,
        mm: 0.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: -1.0,
        mm: 0.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 2.0,
        f: -2.0,
        o: 0.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: -1.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 1.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: -2.0,
        m: 1.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: -1.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 0.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 2.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: -2.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 2.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: -1.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: -1.0,
        mm: 0.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 2.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 0.0,
        mm: 2.0,
        f: 0.0,
        o: 1.0,
    },
    NutationArguments {
        d: -2.0,
        m: 1.0,
        mm: 0.0,
        f: 2.0,
        o: 1.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: -2.0,
        o: 0.0,
    },
    NutationArguments {
        d: -1.0,
        m: 0.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: -2.0,
        m: 1.0,
        mm: 0.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 1.0,
        m: 0.0,
        mm: 0.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 1.0,
        f: 2.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: -2.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: -1.0,
        m: -1.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: 1.0,
        mm: 1.0,
        f: 0.0,
        o: 0.0,
    },
    NutationArguments {
        d: 0.0,
        m: -1.0,
        mm: 1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 2.0,
        m: -1.0,
        mm: -1.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 0.0,
        m: 0.0,
        mm: 3.0,
        f: 2.0,
        o: 2.0,
    },
    NutationArguments {
        d: 2.0,
        m: -1.0,
        mm: 0.0,
        f: 2.0,
        o: 2.0,
    },
];

const COEFFICIENTS: [NutationCoefficients; TERMS] = [
    NutationCoefficients {
        longitude1: -171996.0,
        longitude2: -174.2,
        obliquity1: 92025.0,
        obliquity2: 8.9,
    },
    NutationCoefficients {
        longitude1: -13187.0,
        longitude2: -1.6,
        obliquity1: 5736.0,
        obliquity2: -3.1,
    },
    NutationCoefficients {
        longitude1: -2274.0,
        longitude2: -0.2,
        obliquity1: 977.0,
        obliquity2: -0.5,
    },
    NutationCoefficients {
        longitude1: 2062.0,
        longitude2: 0.2,
        obliquity1: -895.0,
        obliquity2: 0.5,
    },
    NutationCoefficients {
        longitude1: 1426.0,
        longitude2: -3.4,
        obliquity1: 54.0,
        obliquity2: -0.1,
    },
    NutationCoefficients {
        longitude1: 712.0,
        longitude2: 0.1,
        obliquity1: -7.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -517.0,
        longitude2: 1.2,
        obliquity1: 224.0,
        obliquity2: -0.6,
    },
    NutationCoefficients {
        longitude1: -386.0,
        longitude2: -0.4,
        obliquity1: 200.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -301.0,
        longitude2: 0.0,
        obliquity1: 129.0,
        obliquity2: -0.1,
    },
    NutationCoefficients {
        longitude1: 217.0,
        longitude2: -0.5,
        obliquity1: -95.0,
        obliquity2: 0.3,
    },
    NutationCoefficients {
        longitude1: -158.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 129.0,
        longitude2: 0.1,
        obliquity1: -70.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 123.0,
        longitude2: 0.0,
        obliquity1: -53.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 63.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 63.0,
        longitude2: 0.1,
        obliquity1: -33.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -59.0,
        longitude2: 0.0,
        obliquity1: 26.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -58.0,
        longitude2: -0.1,
        obliquity1: 32.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -51.0,
        longitude2: 0.0,
        obliquity1: 27.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 48.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 46.0,
        longitude2: 0.0,
        obliquity1: -24.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -38.0,
        longitude2: 0.0,
        obliquity1: 16.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -31.0,
        longitude2: 0.0,
        obliquity1: 13.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 29.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 29.0,
        longitude2: 0.0,
        obliquity1: -12.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 26.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -22.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 21.0,
        longitude2: 0.0,
        obliquity1: -10.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 17.0,
        longitude2: -0.1,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 16.0,
        longitude2: 0.0,
        obliquity1: -8.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -16.0,
        longitude2: 0.1,
        obliquity1: 7.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -15.0,
        longitude2: 0.0,
        obliquity1: 9.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -13.0,
        longitude2: 0.0,
        obliquity1: 7.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -12.0,
        longitude2: 0.0,
        obliquity1: 6.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 11.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -10.0,
        longitude2: 0.0,
        obliquity1: 5.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -8.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 7.0,
        longitude2: 0.0,
        obliquity1: -3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -7.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -7.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -7.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 6.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 6.0,
        longitude2: 0.0,
        obliquity1: -3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 6.0,
        longitude2: 0.0,
        obliquity1: -3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -6.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -6.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 5.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -5.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -5.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -5.0,
        longitude2: 0.0,
        obliquity1: 3.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -4.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: 3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
    NutationCoefficients {
        longitude1: -3.0,
        longitude2: 0.0,
        obliquity1: 0.0,
        obliquity2: 0.0,
    },
];
