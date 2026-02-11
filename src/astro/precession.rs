// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Equatorial Precession utilities
//!
//! This module implements **precession of the Earth’s mean equator and equinox**
//! following the *short series* formulation in chapter&nbsp;20 of Jean&nbsp;Meeus’
//! *Astronomical Algorithms* (2nd ed.).  The algorithm is sufficient for most
//! practical astronomical work: its errors remain below **0.1 arc‑second** for
//! epochs within ±5 Julian centuries of **J2000.0**. That is, the years
//! **1500 → 2500**.
//!
//! ## What is precession?
//! The rotation axis of the Earth is not fixed in inertial space: the torques
//! produced by the gravitational attraction of the Sun and the Moon on the
//! Earth’s equatorial bulge make the axis describe a **slow conical motion**, a
//! gyroscope under external torque.  The effect is known as **lunisolar
//! precession** and has a period of about **25 770 years** (*the Platonic year*).
//! A smaller contribution, **planetary precession**, is produced by the tidal
//! forces of the other planets.
//!
//! In equatorial coordinates (right ascension *α*, declination *δ*) this causes
//! the celestial poles and the equinox to *drift* at roughly **50″ · yr⁻¹**. If
//! the apparent or mean position of a star is referred to two different epochs
//! it must therefore be **precessed** to the desired date before it can be
//! compared with catalogued data or with another observation.
//!
//! ## Why do we need a numerical model?
//! Astronomical catalogues adopt a fixed **reference epoch** (e.g. B1950.0 or
//! J2000.0). Precise pointing, orbit determination or reduction of
//! observational data to standard coordinates all require a transformation
//! between the catalogue epoch and the observation epoch.  Performing that
//! transformation with an accuracy of a few milliarc‑seconds – small enough for
//! sub‑arc‑second telescopes and for most amateur‑level applications, is the
//! purpose of the present code.
//!
//! ## How is precession computed here?
//! Meeus expresses precession as three **Euler‑like rotation angles**
//! (ζ *zeta*, *z* and θ *theta*) around the axes *z → x → z*.  For a time span of
//! *t* Julian centuries (*T* denotes the starting epoch measured from J2000.0)
//! the angles are expanded as low‑order polynomials in *t* and *T* (**eqs.&nbsp;20.2
//! and 20.3**).  They are first given in arc‑seconds:
//!
//! ```text
//! ζ  = ( 2306.2181 + 1.39656 T − 0.000139 T² ) t
//!      + ( 0.30188 − 0.000344 T ) t² + 0.017998 t³
//! z  = ( 2306.2181 + 1.39656 T − 0.000139 T² ) t
//!      + ( 1.09468 + 0.000066 T ) t² + 0.018203 t³
//! θ  = ( 2004.3109 − 0.85330 T − 0.000217 T² ) t
//!      − ( 0.42665 + 0.000217 T ) t² − 0.041833 t³
//! ```
//!
//! Those angles are converted to radians and the **right‑hand rotation
//! (ζ, θ, −z)** is applied to the input equatorial coordinates using
//! Meeus eq.&nbsp;20.4.
//!
//! ## Precession vs Nutation
//! | Aspect            | **Precession**                                   | **Nutation**                                   |
//! |-------------------|--------------------------------------------------|------------------------------------------------|
//! | Physical cause    | Long‑term torque on the equatorial bulge by Sun and Moon, plus planetary tidal forces | Short‑term periodic torque variations mainly from the Moon’s changing distance, inclination and the Sun’s apparent motion |
//! | Character         | **Secular** (monotonic drift, ≈ 50″ · yr⁻¹)       | **Periodic** (dominant 18.6 yr term of ±9″, plus many smaller terms) |
//! | Time scale        | ~25 770 yr cycle                                 | Hours → decades (hundreds of terms)            |
//! | Typical magnitude | 0.014° per year                                  | up to 9 arc‑seconds peak‑to‑peak               |
//! | Modelling         | Low‑order polynomials (this module)              | Large trigonometric series (IAU 1980/2000B, ELP 2000) |
//! | When to apply?    | **Always** when comparing coordinates referred to different epochs | When sub‑arc‑second precision or true‑of‑date coordinates are required |
//!
//! In other words, **precession is the steady drift** of the reference frame;
//! **nutation is the superposed wobble**. The two effects must be added to
//! obtain the *true* or *apparent* equator and equinox of date.  This module
//! provides the precession step; nutation would be applied subsequently by a
//! dedicated routine.
//!
//! ## Public API
//! * [`precess_from_j2000`] – common convenience for the transformation
//!   *J2000.0 → date*.
//! * [`precess_equatorial`] – general *epoch → epoch* interface.
//!
//! ## Accuracy & limitations
//! The short series is adequate for most handbook‑level work (<0.1″ error up
//! to ±5 centuries).  For **modern astrometry below the milli‑arc‑second
//! level** the full IAU 2006 precession‑nutation model should be preferred.
//!
//! Equatorial precession utilities (Meeus, *Astronomical Algorithms*,
//! 2nd ed., ch. 20).  Accuracy of the **short model** is better than
//! 0.1″ for |T| ≤ 5 centuries.  All angles are handled via the crate’s
//! `Degrees` / `Radians` wrappers; right ascension is wrapped to
//! **0 ≤ α < 2π**.
//!
//! Public API
//! ----------
//! - [`precess_from_j2000`] – convenience layer for the common case
//!   “mean J2000.0 → given date”.
//! - [`precess_equatorial`] – precess between *any* two epochs.

use crate::coordinates::spherical::position;
use crate::time::JulianDate;
use affn::Rotation3;
use qtty::*;

/* -------------------------------------------------------------------------
 * Constants & small utilities
 * ---------------------------------------------------------------------- */

/// 85 ° in radians: threshold above which the `asin` formula for δ loses
/// precision.  Switching to the alternative formula earlier does no harm
/// and guarantees full numerical stability even closer to the pole.
const NEAR_POLE_LIMIT: Radians = Degrees::new(85.0).to_const::<Radian>();

/// Return (*t*, *t²*, *t³*).  Saves three multiplications when the caller
/// needs all powers.
#[inline]
fn t_powers(t: f64) -> (f64, f64, f64) {
    let t2 = t * t;
    (t, t2, t2 * t)
}

/* -------------------------------------------------------------------------
 * Precession coefficients (Meeus 20‑2/20‑3)
 * ---------------------------------------------------------------------- */

/// Compute the *short‑series* precession angles `ζ`, `z` and `θ` *(radians)*
/// for a span of `span` Julian centuries, starting at epoch
/// `epoch` (both measured from J2000.0).
///
/// Follows Meeus Eq. 20.2/20.3.  The trick of dividing **T** and *t* by
/// 3600 converts the published coefficients (arcseconds) into *degrees*,
/// avoiding one extra division.
#[inline]
#[allow(non_snake_case)]
fn short_series_coefficients(epoch: Centuries, span: Centuries) -> (Radians, Radians, Radians) {
    // Meeus uses full Julian centuries, *not* divided by 3600.
    let T = epoch.value();
    let (t, t2, t3) = t_powers(span.value());

    // 20.2 / 20.3 – all results in **arc-seconds**.
    let zeta_as = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (0.30188 - 0.000344 * T) * t2
        + 0.017_998 * t3;

    let z_as = (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (1.09468 + 0.000066 * T) * t2
        + 0.018_203 * t3;

    let theta_as = (2004.3109 - 0.85330 * T - 0.000217 * T * T) * t
        - (0.42665 + 0.000217 * T) * t2
        - 0.041_833 * t3;

    // arc-seconds → degrees → radians
    let as_to_rad = |a: f64| Degrees::new(a / 3600.0).to::<Radian>();
    (as_to_rad(zeta_as), as_to_rad(z_as), as_to_rad(theta_as))
}

#[inline]
fn precession_angles(from_jd: JulianDate, to_jd: JulianDate) -> (Radians, Radians, Radians) {
    let epoch = from_jd.julian_centuries();
    let span = to_jd.julian_centuries() - epoch;
    short_series_coefficients(epoch, span)
}

/* -------------------------------------------------------------------------
 * Core rotation (Meeus 20‑4)
 * ---------------------------------------------------------------------- */

/// Rotate a single equatorial coordinate by the precession angles.
/// Returns *(α, δ)* **in radians**.
#[inline]
fn rotate_equatorial(
    ra: Radians,
    dec: Radians,
    zeta: Radians,
    z: Radians,
    theta: Radians,
) -> (Radians, Radians) {
    // Pre‑compute repeated terms with `sin_cos` for a small perf gain.
    let (sin_ra_zeta, cos_ra_zeta) = (ra + zeta).sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();
    let (sin_th, cos_th) = theta.sin_cos();

    // Meeus Eq. 20.4
    let a = cos_dec * sin_ra_zeta;
    let b = cos_th * cos_dec * cos_ra_zeta - sin_th * sin_dec;
    let c = sin_th * cos_dec * cos_ra_zeta + cos_th * sin_dec;

    // New right ascension, wrapped to 0–2π.
    let new_ra = (Radians::new(a.atan2(b)) + z).normalize();

    // New declination: pole‑safe formula when |δ| > 85 °.
    let new_dec = if dec.abs() > NEAR_POLE_LIMIT {
        let mut d = (a * a + b * b).sqrt().acos();
        if dec < Radians::new(0.0) {
            d = -d;
        }
        d
    } else {
        c.asin()
    };

    (new_ra, Radians::new(new_dec))
}

#[inline]
fn rotation_y(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
}

#[inline]
fn rotation_z(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
}

/* -------------------------------------------------------------------------
 * Public API
 * ---------------------------------------------------------------------- */

/// **J2000.0 → to_jd** shortcut.
///
/// Transforms mean right ascension & declination supplied at the standard
/// reference epoch (JD 2 451 545.0) to the mean equator & equinox of `to_jd`.
/// The `JulianDate` input is interpreted as TT.
#[inline]
pub fn precess_from_j2000<U: LengthUnit>(
    mean_position: position::EquatorialMeanJ2000<U>,
    to_jd: JulianDate,
) -> position::EquatorialMeanOfDate<U> {
    let mean_of_date = position::EquatorialMeanOfDate::<U>::new(
        mean_position.ra(),
        mean_position.dec(),
        mean_position.distance,
    );
    precess_equatorial(mean_of_date, JulianDate::J2000, to_jd)
}

/// **Epoch → epoch** precession.
///
/// * `position` – (α, δ, r) referred to the mean equator/equinox of `from_jd`.
/// * `from_jd`, `to_jd` – Julian Days of source and target epochs.
///
/// Returns the coordinates referred to `to_jd`.
/// The `JulianDate` inputs are interpreted as TT.
pub fn precess_equatorial<U: LengthUnit>(
    position: position::EquatorialMeanOfDate<U>,
    from_jd: JulianDate,
    to_jd: JulianDate,
) -> position::EquatorialMeanOfDate<U> {
    let ra0 = position.ra().to::<Radian>();
    let dec0 = position.dec().to::<Radian>();

    let (zeta, z, theta) = precession_angles(from_jd, to_jd);

    let (ra, dec) = rotate_equatorial(ra0, dec0, zeta, z, theta);

    position::EquatorialMeanOfDate::<U>::new(
        ra.to::<Degree>(),
        dec.to::<Degree>(),
        position.distance,
    )
}

/// Precession rotation matrix from `from_jd` mean equator/equinox to `to_jd`.
///
/// The input and output frames are **mean equator/equinox** (nutation removed).
/// The `JulianDate` inputs are interpreted as TT.
pub fn precession_rotation(from_jd: JulianDate, to_jd: JulianDate) -> Rotation3 {
    let (zeta, z, theta) = precession_angles(from_jd, to_jd);

    // Meeus Eq. 20.4 corresponds to Rz(z) * Ry(-theta) * Rz(zeta).
    rotation_z(z) * rotation_y(-theta) * rotation_z(zeta)
}

/// Convenience precession rotation from J2000.0 mean equator/equinox to `to_jd`.
///
/// The `JulianDate` input is interpreted as TT.
#[inline]
pub fn precession_rotation_from_j2000(to_jd: JulianDate) -> Rotation3 {
    precession_rotation(JulianDate::J2000, to_jd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{cartesian, spherical};
    use qtty::{AstronomicalUnit, Degrees};

    #[test]
    fn sirius_2025() {
        // Sirius: α0 ≈ 101.287°, δ0 ≈ −16.716° (J2000.0, Hipparcos)
        use crate::bodies::catalog::SIRIUS;

        // Target epoch: 2025‑05‑12 (JD 2 469 807.5)
        let prec = precess_from_j2000(
            SIRIUS.target.get_position().clone(),
            JulianDate::new(2_460_807.5),
        );

        // Expected (Meeus short model): α ≈ 101.84557°, δ ≈ −16.77182°
        assert!(
            (prec.ra() - Degrees::new(101.57047)).abs() < Degrees::new(3e-5),
            "current α ≈ {}",
            prec.ra()
        );
        assert!(
            (prec.dec() - Degrees::new(-16.74409)).abs() < Degrees::new(1e-5),
            "current δ ≈ {}",
            prec.dec()
        );
    }

    #[test]
    fn rotate_equatorial_handles_near_pole_branch() {
        use qtty::Radians;

        let (ra, dec) = super::rotate_equatorial(
            Radians::new(0.1),
            Radians::new(-1.5), // |δ| > 85° to trigger special path
            Radians::new(0.01),
            Radians::new(0.02),
            Radians::new(0.03),
        );

        assert!(ra.is_finite());
        assert!(dec < 0.0);
    }

    #[test]
    fn precession_rotation_matches_spherical() {
        let jd = JulianDate::new(2_460_000.5);
        let pos = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            Degrees::new(120.0),
            Degrees::new(-30.0),
            1.0,
        );

        let prec = precess_from_j2000(pos.clone(), jd);
        let dir = pos.direction().to_cartesian();
        let rot = precession_rotation_from_j2000(jd);
        let [x, y, z] = rot.apply_array([dir.x(), dir.y(), dir.z()]);
        let dir_rot = cartesian::direction::EquatorialMeanOfDate::normalize(x, y, z);
        let sph_rot = spherical::Direction::from_cartesian(&dir_rot);

        let ra_diff = sph_rot.ra().abs_separation(prec.ra());
        let dec_diff = (sph_rot.dec() - prec.dec()).abs();

        assert!(ra_diff < Degrees::new(1e-10), "RA mismatch: {}", ra_diff);
        assert!(dec_diff < Degrees::new(1e-10), "Dec mismatch: {}", dec_diff);
    }
}
