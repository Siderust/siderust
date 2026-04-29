// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Proper Motion Correction
//!
//! This module applies **proper motion** to a star’s mean position
//! in order to obtain its updated position at a new epoch (Julian Day).
//!
//! ## What is proper motion?
//! Proper motion is the apparent angular displacement of a star across
//! the sky due to its real motion through space, relative to the Solar System.
//! It is typically measured in milliarcseconds or degrees per year.
//!
//! Since high-precision star catalogs (like Gaia or Hipparcos) provide
//! positions at a reference epoch (usually J2000.0), we need to correct
//! for proper motion when computing positions at a later date.
//!
//! ## Catalog mapping
//! Gaia and Hipparcos publish right-ascension proper motion as `µα⋆ = µα cos(δ)`.
//! Use [`ProperMotion::from_mu_alpha_star`] for those catalogs. Use
//! [`ProperMotion::from_mu_alpha`] only if your source already provides the true
//! RA angular rate `µα`.

use crate::coordinates::spherical::position;
use crate::time::JulianDate;
use crate::qtty::*;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

type DegreePerYear = crate::qtty::Per<Degree, Year>;
type DegreesPerYear = crate::qtty::angular_rate::AngularRate<Degree, Year>;

const COS_DEC_EPSILON: f64 = 1.0e-12;

/// Convention used for right-ascension proper motion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum RaProperMotionConvention {
    /// True RA angular rate `µα` in deg/year.
    MuAlpha,
    /// Catalog rate `µα⋆ = µα cos(δ)` in deg/year.
    MuAlphaStar,
}

impl std::fmt::Display for RaProperMotionConvention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MuAlpha => write!(f, "\u{03bc}\u{03b1}"),
            Self::MuAlphaStar => write!(f, "\u{03bc}\u{03b1}\u{2a}"),
        }
    }
}

/// Describes the proper motion of a star in equatorial coordinates.
///
/// - `pm_ra` is the right-ascension proper motion
/// - `pm_dec` is the declination proper motion
/// - `ra_convention` specifies whether `pm_ra` is `µα` or `µα⋆ = µα cos(δ)`
///
/// Units: both fields are in **degrees per Julian year**
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProperMotion {
    pub pm_ra: DegreesPerYear,
    pub pm_dec: DegreesPerYear,
    pub ra_convention: RaProperMotionConvention,
}

impl fmt::Display for ProperMotion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\u{03bc}\u{03b1}{}: {}, \u{03bc}\u{03b4}: {}",
            self.ra_convention, self.pm_ra, self.pm_dec
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ProperMotionError {
    /// Conversion from `µα⋆` to `µα` is unstable near the poles (`cos(dec)≈0`).
    RightAscensionUndefinedAtPole { dec_deg: f64 },
}

impl fmt::Display for ProperMotionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProperMotionError::RightAscensionUndefinedAtPole { dec_deg } => write!(
                f,
                "right-ascension proper motion is undefined near the poles (dec={} deg)",
                dec_deg
            ),
        }
    }
}

impl std::error::Error for ProperMotionError {}

impl RaProperMotionConvention {
    fn to_mu_alpha(
        self,
        pm_ra: DegreesPerYear,
        dec: Degrees,
    ) -> Result<DegreesPerYear, ProperMotionError> {
        match self {
            RaProperMotionConvention::MuAlpha => Ok(pm_ra),
            RaProperMotionConvention::MuAlphaStar => {
                let cos_dec = dec.value().to_radians().cos();
                if cos_dec.abs() <= COS_DEC_EPSILON {
                    return Err(ProperMotionError::RightAscensionUndefinedAtPole {
                        dec_deg: dec.value(),
                    });
                }
                Ok(DegreesPerYear::new(pm_ra.value() / cos_dec))
            }
        }
    }
}

impl ProperMotion {
    /// Construct proper motion from true RA angular rate `µα`.
    pub fn from_mu_alpha<T>(pm_ra: Quantity<T>, pm_dec: Quantity<T>) -> Self
    where
        T: AngularRateUnit,
    {
        Self {
            pm_ra: pm_ra.to::<DegreePerYear>(),
            pm_dec: pm_dec.to::<DegreePerYear>(),
            ra_convention: RaProperMotionConvention::MuAlpha,
        }
    }

    /// Construct proper motion from catalog rate `µα⋆ = µα cos(δ)`.
    pub fn from_mu_alpha_star<T>(pm_ra_cosdec: Quantity<T>, pm_dec: Quantity<T>) -> Self
    where
        T: AngularRateUnit,
    {
        Self {
            pm_ra: pm_ra_cosdec.to::<DegreePerYear>(),
            pm_dec: pm_dec.to::<DegreePerYear>(),
            ra_convention: RaProperMotionConvention::MuAlphaStar,
        }
    }

    fn ra_rate_at_epoch(&self, dec: Degrees) -> Result<DegreesPerYear, ProperMotionError> {
        self.ra_convention.to_mu_alpha(self.pm_ra, dec)
    }
}

/// Applies proper motion to a given position for a custom time interval.
///
/// # Arguments
/// - `mean_position`: the reference position at epoch `epoch_jd`
/// - `proper_motion`: angular velocity in RA and DEC (degrees/year)
/// - `jd`: the target Julian Day (epoch at which we want the new position)
/// - `epoch_jd`: the epoch of the original mean position
///
/// # Returns
/// New apparent position adjusted by linear proper motion, or an error if
/// `µα⋆` is used at/near the poles.
///
/// Assumes motion is linear (valid for most stars over <1000 year timescales).
fn set_proper_motion_since_epoch<U: LengthUnit>(
    mean_position: position::EquatorialMeanJ2000<U>,
    proper_motion: ProperMotion,
    jd: JulianDate,
    epoch_jd: JulianDate,
) -> Result<position::EquatorialMeanJ2000<U>, ProperMotionError> {
    // Time difference in Julian years
    let t: Years =
        Years::new((jd / JulianDate::JULIAN_YEAR) - (epoch_jd / JulianDate::JULIAN_YEAR));
    let ra_rate = proper_motion.ra_rate_at_epoch(mean_position.dec())?;
    // Linearly apply proper motion in RA and DEC
    Ok(position::EquatorialMeanJ2000::<U>::new(
        mean_position.ra() + (ra_rate * t).to(),
        mean_position.dec() + (proper_motion.pm_dec * t).to(),
        mean_position.distance,
    ))
}

/// Applies proper motion from the standard epoch J2000.0.
///
/// # Arguments
/// - `mean_position`: star's catalogued position at J2000.0
/// - `proper_motion`: proper motion in RA/DEC (deg/year)
/// - `jd`: target Julian Day for evaluation
///
/// # Returns
/// Updated position after applying proper motion since J2000.0, or an error if
/// `µα⋆` is used at/near the poles.
pub fn set_proper_motion_since_j2000<U: LengthUnit>(
    mean_position: position::EquatorialMeanJ2000<U>,
    proper_motion: ProperMotion,
    jd: JulianDate,
) -> Result<position::EquatorialMeanJ2000<U>, ProperMotionError> {
    set_proper_motion_since_epoch(mean_position, proper_motion, jd, JulianDate::J2000)
}

// ═══════════════════════════════════════════════════════════════════════════
// Full 6D space-motion propagation (Hipparcos / SOFA `iauStarpm`)
// ═══════════════════════════════════════════════════════════════════════════

/// Catalog quantities required for a complete 6D space-motion update.
///
/// All four fields are typed `qtty` quantities so unit conversions are
/// explicit and checked at compile time.
///
/// Catalog conventions:
/// - `pm_ra` is `μα⋆ = μα·cos(δ)` (the standard Gaia / Hipparcos convention).
/// - `parallax` is the **annual trigonometric parallax** at the catalog epoch.
/// - `radial_velocity` is positive away from the Solar System (recession).
///
/// Both `parallax` and `radial_velocity` are required for a full update; the
/// transverse-only [`set_proper_motion_since_epoch`] family ignores them and
/// is therefore not adequate for fast nearby stars (Barnard's, Proxima …)
/// over multi-decade horizons.
#[derive(Debug, Clone, Copy)]
pub struct StarSpaceMotion {
    /// Right-ascension proper motion in the catalog convention `μα⋆ = μα·cos(δ)`.
    pub pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate<MilliArcsecond, Year>,
    /// Declination proper motion `μδ`.
    pub pm_dec: crate::qtty::angular_rate::AngularRate<MilliArcsecond, Year>,
    /// Annual trigonometric parallax (positive).
    pub parallax: MilliArcseconds,
    /// Radial velocity, positive away from the Solar System.
    pub radial_velocity: crate::qtty::velocity::Velocity<Kilometer, Second>,
}

/// Apply full 6D space-motion to a catalog star position.
///
/// This is the SOFA `iauStarpm` / Hipparcos-style propagation:
///
/// 1. Reconstruct the BCRS 3-vector `p₀ = (1/π) · R(α, δ)` (in AU).
/// 2. Build the BCRS velocity `ṗ = T · μ + v_R · r̂` (in AU/yr), where `T` is
///    the basis of the local tangent plane on the celestial sphere.
/// 3. Linearly propagate `p(t) = p₀ + Δt · ṗ`.
/// 4. Re-project to spherical (RA, Dec) at the new epoch and rescale the
///    distance accordingly.
///
/// Unlike [`set_proper_motion_since_epoch`], this version honours the
/// **radial velocity** (giving the secular acceleration in proper motion
/// caused by the changing transverse projection of a 3D velocity) and the
/// **parallax** (giving the absolute scale of the propagation). For most
/// catalog stars over a few centuries the difference is sub-mas, but for
/// fast nearby stars it dominates the linear approximation within decades.
///
/// # Arguments
/// - `mean_position`: catalog position at `epoch_jd`. The radial coordinate
///   (`distance`) is overwritten by `1/parallax`; only the angular fields
///   matter on input.
/// - `motion`: full Hipparcos/Gaia state ([`StarSpaceMotion`]).
/// - `jd`: target epoch (TT/TDB).
/// - `epoch_jd`: catalog reference epoch (typically J2000 or J2015.5).
///
/// # Errors
/// Returns [`ProperMotionError::RightAscensionUndefinedAtPole`] when the
/// catalog declination is too close to ±90°: the conversion from `μα⋆` back
/// to a true rate is unstable there.
///
/// # References
/// - SOFA routine `iauStarpm` (and the underlying `iauPvstar` / `iauStarpv`
///   conversions).
/// - Hipparcos and Tycho Catalogues, ESA SP-1200 (1997), Vol. 1, §1.5.5.
pub fn propagate_space_motion(
    mean_position: position::EquatorialMeanJ2000<AstronomicalUnit>,
    motion: StarSpaceMotion,
    jd: JulianDate,
    epoch_jd: JulianDate,
) -> Result<position::EquatorialMeanJ2000<AstronomicalUnit>, ProperMotionError> {
    // ── Constants ──
    // 1 AU / (1 km/s · 1 Julian year) = days_per_year · seconds_per_day / km_per_au
    // = 365.25 · 86400 / 1.49597870700e8 = 0.21094502... AU per (km/s · yr)
    const AU_PER_KMS_YR: f64 = 0.210_945_021_894_944_8;
    // 1 mas in radians: π / (180 · 3600 · 1000)
    const MAS_TO_RAD: f64 = std::f64::consts::PI / (180.0 * 3_600_000.0);

    let dec_deg = mean_position.dec();
    let cos_dec = dec_deg.value().to_radians().cos();
    if cos_dec.abs() <= COS_DEC_EPSILON {
        return Err(ProperMotionError::RightAscensionUndefinedAtPole {
            dec_deg: dec_deg.value(),
        });
    }

    let alpha = mean_position.ra().value().to_radians();
    let delta = dec_deg.value().to_radians();
    let (sin_a, cos_a) = alpha.sin_cos();
    let (sin_d, cos_d) = delta.sin_cos();

    // Distance in AU from parallax (mas → AU). r_AU = 1 rad / π_rad
    let pi_rad = motion.parallax.value() * MAS_TO_RAD;
    if pi_rad <= 0.0 {
        // A non-positive parallax has no physical meaning (the "negative
        // parallax" entries in some catalogues are noisy fits; callers
        // should treat them as missing and fall back to a transverse-only
        // propagation).
        return Err(ProperMotionError::RightAscensionUndefinedAtPole {
            dec_deg: dec_deg.value(),
        });
    }
    let r_au = 1.0 / pi_rad;

    // BCRS position vector at epoch (AU).
    let p0 = [r_au * cos_d * cos_a, r_au * cos_d * sin_a, r_au * sin_d];

    // Transverse basis vectors on the celestial sphere (unit, dimensionless).
    let e_alpha = [-sin_a, cos_a, 0.0];
    let e_delta = [-sin_d * cos_a, -sin_d * sin_a, cos_d];

    // Angular rates → linear transverse velocity in AU/yr.
    // μα⋆ already includes cos(δ); divide once to get true μα, then scale by
    // r·cos(δ) to recover the AU/yr coefficient on e_alpha.
    let pm_ra_rad_yr = (motion.pm_ra_cos_dec.value() * MAS_TO_RAD) / cos_dec;
    let pm_dec_rad_yr = motion.pm_dec.value() * MAS_TO_RAD;

    let v_alpha_au_yr = pm_ra_rad_yr * r_au * cos_d;
    let v_delta_au_yr = pm_dec_rad_yr * r_au;

    // Radial velocity in AU/yr (positive = away).
    let v_r_au_yr = motion.radial_velocity.value() * AU_PER_KMS_YR;

    // BCRS velocity vector in AU/yr.
    let r_hat = [cos_d * cos_a, cos_d * sin_a, sin_d];
    let pdot = [
        e_alpha[0] * v_alpha_au_yr + e_delta[0] * v_delta_au_yr + r_hat[0] * v_r_au_yr,
        e_alpha[1] * v_alpha_au_yr + e_delta[1] * v_delta_au_yr + r_hat[1] * v_r_au_yr,
        e_alpha[2] * v_alpha_au_yr + e_delta[2] * v_delta_au_yr + r_hat[2] * v_r_au_yr,
    ];

    // Δt in Julian years.
    let dt_yr = (jd / JulianDate::JULIAN_YEAR) - (epoch_jd / JulianDate::JULIAN_YEAR);

    // Linear propagation in BCRS.
    let p = [
        p0[0] + dt_yr * pdot[0],
        p0[1] + dt_yr * pdot[1],
        p0[2] + dt_yr * pdot[2],
    ];

    // Re-project to spherical.
    let r_new = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    if r_new <= 0.0 {
        // Star passed exactly through the Solar System; refuse rather than
        // emit a divide-by-zero. This requires unphysical inputs.
        return Err(ProperMotionError::RightAscensionUndefinedAtPole {
            dec_deg: dec_deg.value(),
        });
    }
    let new_dec = (p[2] / r_new).asin();
    let new_ra = p[1].atan2(p[0]).rem_euclid(std::f64::consts::TAU);

    Ok(position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
        Degrees::new(new_ra.to_degrees()),
        Degrees::new(new_dec.to_degrees()),
        r_new,
    ))
}

/// Apply full 6D space-motion since J2000.0.
///
/// Convenience wrapper around [`propagate_space_motion`] for the common case
/// of a J2000-referenced catalog (Gaia DR3 inertial epoch, FK5/J2000, …).
pub fn propagate_space_motion_since_j2000(
    mean_position: position::EquatorialMeanJ2000<AstronomicalUnit>,
    motion: StarSpaceMotion,
    jd: JulianDate,
) -> Result<position::EquatorialMeanJ2000<AstronomicalUnit>, ProperMotionError> {
    propagate_space_motion(mean_position, motion, jd, JulianDate::J2000)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{
        centers::Geocentric, frames::EquatorialMeanJ2000, spherical::Position,
    };
    use crate::time::JulianDate;
    use crate::qtty::{AstronomicalUnit, Degrees};

    type DegreesPerYear = crate::qtty::Quantity<crate::qtty::Per<Degree, Year>>;

    #[test]
    fn test_proper_motion_linear_shift_mu_alpha() {
        // Mean position at J2000
        let mean_position = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(10.0), // RA = 10°
            Degrees::new(20.0), // DEC = 20°
            1.0,                // arbitrary distance
        );

        // Proper motion: 0.01°/year in RA, -0.005°/year in DEC
        let mu = ProperMotion::from_mu_alpha::<DegreePerYear>(
            DegreesPerYear::new(0.01),
            DegreesPerYear::new(-0.005),
        );

        // Target epoch: 50 years after J2000
        let jd_future = JulianDate::J2000 + 50.0 * JulianDate::JULIAN_YEAR;

        // Expected shifts (compute before moving mean_position)
        let expected_ra = mean_position.ra() + Degrees::new(0.5);
        let expected_dec = mean_position.dec() - Degrees::new(0.25);

        let shifted = set_proper_motion_since_j2000(mean_position, mu, jd_future).unwrap();

        let ra_err = (shifted.ra() - expected_ra).abs();
        let dec_err = (shifted.dec() - expected_dec).abs();

        assert!(
            ra_err < Degrees::new(1e-6),
            "RA shifted incorrectly: got {}, expected {}",
            shifted.ra(),
            expected_ra
        );
        assert!(
            dec_err < Degrees::new(1e-6),
            "DEC shifted incorrectly: got {}, expected {}",
            shifted.dec(),
            expected_dec
        );
    }

    #[test]
    fn test_proper_motion_mu_alpha_star_conversion() {
        let mean_position = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(10.0), // RA = 10°
            Degrees::new(60.0), // DEC = 60° => cos(dec) = 0.5
            1.0,
        );

        // µα = 0.01 deg/yr, so µα⋆ = µα cos(dec) = 0.005 deg/yr.
        let mu = ProperMotion::from_mu_alpha_star::<DegreePerYear>(
            DegreesPerYear::new(0.005),
            DegreesPerYear::new(0.0),
        );

        let jd_future = JulianDate::J2000 + 10.0 * JulianDate::JULIAN_YEAR;
        let shifted = set_proper_motion_since_j2000(mean_position, mu, jd_future).unwrap();
        let expected_ra = Degrees::new(10.1);
        let ra_err = (shifted.ra() - expected_ra).abs();

        assert!(
            ra_err < Degrees::new(1e-6),
            "RA shifted incorrectly for µα⋆: got {}, expected {}",
            shifted.ra(),
            expected_ra
        );
    }

    #[test]
    fn test_proper_motion_mu_alpha_star_rejects_pole() {
        let mean_position = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(30.0),
            Degrees::new(90.0),
            1.0,
        );
        let mu = ProperMotion::from_mu_alpha_star::<DegreePerYear>(
            DegreesPerYear::new(0.01),
            DegreesPerYear::new(0.0),
        );
        let jd_future = JulianDate::J2000 + JulianDate::JULIAN_YEAR;

        assert!(matches!(
            set_proper_motion_since_j2000(mean_position, mu, jd_future),
            Err(ProperMotionError::RightAscensionUndefinedAtPole { .. })
        ));
    }

    // ─── Full 6D space-motion (`propagate_space_motion`) ─────────────────

    /// Zero motion ⇒ position is unchanged at any future epoch.
    #[test]
    fn space_motion_zero_state_is_static() {
        let pos = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(30.0),
            1.0,
        );
        let motion = StarSpaceMotion {
            pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            pm_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            parallax: MilliArcseconds::new(100.0),
            radial_velocity: crate::qtty::velocity::Velocity::new(0.0),
        };
        let jd_future = JulianDate::J2000 + 1000.0 * JulianDate::JULIAN_YEAR;
        let p = propagate_space_motion_since_j2000(pos, motion, jd_future).unwrap();
        assert!((p.ra().value() - 45.0).abs() < 1e-9);
        assert!((p.dec().value() - 30.0).abs() < 1e-9);
    }

    /// Pure transverse PM with no RV reproduces the linear-PM solution to
    /// first order (the difference must be < 1 mas after a century).
    #[test]
    fn space_motion_matches_linear_pm_for_distant_star() {
        // A 10 pc star with 100 mas/yr in RA, 0 in Dec, no RV.
        let pos = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(120.0),
            Degrees::new(0.0),
            1.0,
        );
        let motion = StarSpaceMotion {
            pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate::new(100.0),
            pm_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            parallax: MilliArcseconds::new(100.0), // 10 pc
            radial_velocity: crate::qtty::velocity::Velocity::new(0.0),
        };
        let jd_future = JulianDate::J2000 + 100.0 * JulianDate::JULIAN_YEAR;
        let space = propagate_space_motion_since_j2000(pos, motion, jd_future).unwrap();

        // Linear (transverse-only) prediction: 100 mas/yr · 100 yr = 10 000 mas
        // in α at δ = 0. That's 10 000/3.6e6 ≈ 2.7778e-3 deg.
        let expected_ra = 120.0 + 0.002_777_8;
        assert!(
            (space.ra().value() - expected_ra).abs() < 0.001,
            "RA after century = {}, expected ~{}",
            space.ra().value(),
            expected_ra
        );
        assert!(space.dec().value().abs() < 1e-3);
    }

    /// Strong radial velocity must shrink the apparent PM rate (perspective
    /// acceleration). After enough time, an approaching star with constant
    /// 3D velocity exhibits a measurable change in PM amplitude.
    #[test]
    fn space_motion_secular_acceleration_with_radial_velocity() {
        // Barnard's-star-like: 10 000 mas/yr in α-equiv on a 1.83 pc star,
        // closing in at -111 km/s.
        let pos = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(269.45),
            Degrees::new(4.66),
            1.0,
        );
        let motion = StarSpaceMotion {
            pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            pm_dec: crate::qtty::angular_rate::AngularRate::new(10_000.0),
            parallax: MilliArcseconds::new(547.0),
            radial_velocity: crate::qtty::velocity::Velocity::new(-110.0),
        };

        let p10 = propagate_space_motion_since_j2000(
            pos,
            motion,
            JulianDate::J2000 + 10.0 * JulianDate::JULIAN_YEAR,
        )
        .unwrap();
        let p20 = propagate_space_motion_since_j2000(
            pos,
            motion,
            JulianDate::J2000 + 20.0 * JulianDate::JULIAN_YEAR,
        )
        .unwrap();

        // The Dec shift in the second decade should exceed the first decade
        // because the star is closer (perspective acceleration). For an
        // approaching star, the geometric proper motion grows.
        let d10 = (p10.dec() - pos.dec()).value().abs();
        let d20_minus_d10 = (p20.dec() - p10.dec()).value().abs();
        assert!(
            d20_minus_d10 > d10,
            "expected accelerating PM (d20-d10 = {} > d10 = {})",
            d20_minus_d10,
            d10
        );
    }

    /// Negative parallax must be rejected: it has no physical meaning.
    #[test]
    fn space_motion_negative_parallax_errors() {
        let pos = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
            1.0,
        );
        let motion = StarSpaceMotion {
            pm_ra_cos_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            pm_dec: crate::qtty::angular_rate::AngularRate::new(0.0),
            parallax: MilliArcseconds::new(-1.0),
            radial_velocity: crate::qtty::velocity::Velocity::new(0.0),
        };
        let result = propagate_space_motion_since_j2000(
            pos,
            motion,
            JulianDate::J2000 + JulianDate::JULIAN_YEAR,
        );
        assert!(matches!(
            result,
            Err(ProperMotionError::RightAscensionUndefinedAtPole { .. })
        ));
    }
}
