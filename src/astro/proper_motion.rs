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

use crate::units::*;
use crate::coordinates::spherical::position;

/// Describes the proper motion of a star in equatorial coordinates.
///
/// - `ra_μ` is the proper motion in **right ascension**, typically including `cos(δ)`
/// - `dec_μ` is the proper motion in **declination**
///
/// Units: both fields are in **degrees per Julian year**
#[derive(Debug, Clone)]
pub struct ProperMotion {
    pub ra_μ: DmsPerYear,
    pub dec_μ: DmsPerYear,
}

impl ProperMotion {
    pub fn from_mas_per_year(ra: f64, dec: f64) -> ProperMotion {
        ProperMotion{
            ra_μ: DmsPerYear(DMS::from_milliarcseconds(ra), YEAR),
            dec_μ: DmsPerYear(DMS::from_milliarcseconds(dec), YEAR)
        }
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
/// New apparent position adjusted by linear proper motion.
///
/// Assumes motion is linear (valid for most stars over <1000 year timescales).
fn set_proper_motion_since_epoch<U: LengthUnit>(
    mean_position: position::Equatorial<U>,
    proper_motion: ProperMotion,
    jd: JulianDay,
    epoch_jd: JulianDay
) -> position::Equatorial<U> {
    // Time difference in Julian years
    let t: Years = Years::new((jd - epoch_jd) / JulianDay::JULIAN_YEAR);
    // Linearly apply proper motion in RA and DEC
    position::Equatorial::new::<Quantity<U>>(
        mean_position.ra() + (proper_motion.ra_μ * t).to_degrees().normalize(),
        (mean_position.dec() + (proper_motion.dec_μ * t).to_degrees()).normalize(),
        mean_position.distance.unwrap(),
    )
}

/// Applies proper motion from the standard epoch J2000.0.
///
/// # Arguments
/// - `mean_position`: star's catalogued position at J2000.0
/// - `proper_motion`: proper motion in RA/DEC (deg/year)
/// - `jd`: target Julian Day for evaluation
///
/// # Returns
/// Updated position after applying proper motion since J2000.0
pub fn set_proper_motion_since_j2000<U: LengthUnit>(
    mean_position: position::Equatorial<U>,
    proper_motion: ProperMotion,
    jd: JulianDay
) -> position::Equatorial<U> {
    set_proper_motion_since_epoch(mean_position, proper_motion, jd, JulianDay::J2000)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{Degrees, DmsPerYear, JulianDay, AstronomicalUnit};
    use crate::coordinates::{
        spherical::Position,
        centers::Geocentric,
        frames::Equatorial
    };

    #[test]
    fn test_proper_motion_linear_shift() {
        // Mean position at J2000
        let mean_position = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(
            Degrees::new(10.0),   // RA = 10°
            Degrees::new(20.0),  // DEC = 20°
            1.0      // arbitrary distance
        );

        // Proper motion: 0.01°/year in RA, -0.005°/year in DEC
        let mu = ProperMotion {
            ra_μ: DmsPerYear::from_decimal(0.01),
            dec_μ: DmsPerYear::from_decimal(-0.005),
        };

        // Target epoch: 50 years after J2000
        let jd_future = JulianDay::J2000 + 50.0 * JulianDay::JULIAN_YEAR;

        let shifted = set_proper_motion_since_j2000(mean_position, mu, jd_future);

        // Expected shifts
        let expected_ra = mean_position.ra()   + Degrees::new(0.5);
        let expected_dec = mean_position.dec() - Degrees::new(0.25);

        let ra_err = (shifted.ra() - expected_ra).abs();
        let dec_err = (shifted.dec() - expected_dec).abs();

        assert!(ra_err.as_f64() < 1e-6, "RA shifted incorrectly: got {}, expected {}", shifted.ra(), expected_ra);
        assert!(dec_err.as_f64() < 1e-6, "DEC shifted incorrectly: got {}, expected {}", shifted.dec(), expected_dec);
    }
}
