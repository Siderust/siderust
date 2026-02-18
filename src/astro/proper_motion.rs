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
use qtty::*;
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

type DegreePerYear = qtty::Per<Degree, Year>;
type DegreesPerYear = qtty::frequency::Frequency<Degree, Year>;

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
        T: FrequencyUnit,
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
        T: FrequencyUnit,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{
        centers::Geocentric, frames::EquatorialMeanJ2000, spherical::Position,
    };
    use crate::time::JulianDate;
    use qtty::{AstronomicalUnit, Degrees};

    type DegreesPerYear = qtty::Quantity<qtty::Per<Degree, Year>>;

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
            ra_err < 1e-6,
            "RA shifted incorrectly: got {}, expected {}",
            shifted.ra(),
            expected_ra
        );
        assert!(
            dec_err < 1e-6,
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
            ra_err < 1e-6,
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
}
