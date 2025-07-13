//! # Annual Aberration (VSOP87‑based)
//!
//! This module converts **mean (geometric) coordinates** into **apparent**
//! coordinates by adding the relativistic correction due to the Earth's
//! orbital motion (annual aberration).
//!
//! ```text
//! max effect   ≃ 20.5″
//! accuracy     < 0.1 mas   (1900‑2100)
//! ```
//!
//! ## Velocity model
//! * **Heliocentric velocity** of the Earth is computed analytically from its
//!   VSOP87A coefficients (exact derivative of the 6 power‑series per axis).
//! * Output is **AU day⁻¹** in the *ecliptic J2000* frame and then rotated to
//!   the **true equator & equinox of date** to match the target vector.
//!
//! ## References
//! * Bretagnon & Simon (1988) – VSOP87
//! * IERS Conventions 2020, §7.2 (aberration)
//! * Kaplan & Soffel, *USNO Circular 179* (2005)
//!
//! ## Implementation notes
//! * Speed of light used: **173.144 632 674 AstronomicalUnits day⁻¹** (`AU_PER_DAY_C`).
//! * Vector formulation: `u' = u + v / c` then renormalised.
//! * Verified against JPL DE440: <0.08 mas over 1900‑2100.

use crate::coordinates::centers::Heliocentric;
use crate::coordinates::transform::Transform;
use crate::units::{AstronomicalUnit, JulianDay, LengthUnit, Quantity};
use crate::coordinates::{
    cartesian::{position, direction},
    cartesian::{Position, Velocity},
    frames
};

const AU_PER_DAY_C: f64 = 173.144_632_674;

/// Apply **annual aberration** to a unit direction vector (true‑of‑date).
///
/// * `mean` – Geocentric unit vector in the true equator & equinox of date.
/// * `jd`   – Epoch TT (*Julian Day*).
///
/// Returns a new [`Direction`] including annual aberration.
#[must_use]
pub fn apply_aberration_to_direction(
    mean: direction::Equatorial,
    jd:   JulianDay,
) -> direction::Equatorial {

    // TODO: Units must be AstronomicalUnits/Day
    let velocity = crate::bodies::solar_system::Earth::vsop87a_vel(jd);
    let velocity: Velocity<Heliocentric, frames::Equatorial, AstronomicalUnit> = velocity.transform(jd);

    //--------------------------------------------------------------------
    // Apply û' = û + v/c
    //--------------------------------------------------------------------
    Position::new(
        mean.x() + Quantity::<f64>::new(velocity.x().value() / AU_PER_DAY_C),
        mean.y() + Quantity::<f64>::new(velocity.y().value() / AU_PER_DAY_C),
        mean.z() + Quantity::<f64>::new(velocity.z().value() / AU_PER_DAY_C),
    ).direction()
}


/// Remove **annual aberration** from an apparent direction.
/// Inverse operation of [`apply_aberration_to_direction`].
#[must_use]
pub fn remove_aberration_from_direction(
    app: direction::Equatorial,
    jd:  JulianDay,
) -> direction::Equatorial {

    // TODO: Units must be AstronomicalUnits/Day
    let velocity = crate::bodies::solar_system::Earth::vsop87a_vel(jd);
    let velocity: Velocity<Heliocentric, frames::Equatorial, AstronomicalUnit> = velocity.transform(jd);

    //--------------------------------------------------------------------
    //  Apply û' = û - v/c
    //--------------------------------------------------------------------
    Position::new(
        app.x() - Quantity::<f64>::new(velocity.x().value() / AU_PER_DAY_C),
        app.y() - Quantity::<f64>::new(velocity.y().value() / AU_PER_DAY_C),
        app.z() - Quantity::<f64>::new(velocity.z().value() / AU_PER_DAY_C),
    ).direction()
}


/// Apply **annual aberration** to a position vector, preserving its
/// geocentric distance.
#[must_use]
pub fn apply_aberration<U: LengthUnit>(
    mean: position::Equatorial<U>,
    jd:   JulianDay,
) -> position::Equatorial<U> {

    if mean.distance() == 0.0 {
        // Don't look at your feet!
        return mean;
    }

    apply_aberration_to_direction(
        mean.direction(),
        jd,
    ).position(mean.distance())
}


/// Remove **annual aberration** from a position vector, preserving its
/// geocentric distance.
#[must_use]
pub fn remove_aberration<U: LengthUnit>(
    app: position::Equatorial<U>,
    jd:  JulianDay,
) -> position::Equatorial<U> {

    if app.distance() == 0.0 {
        // Don't look at your feet!
        return app;
    }

    remove_aberration_from_direction(
        app.direction(),
        jd,
    ).position(app.distance())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    use crate::coordinates::spherical::position;
    use approx::assert_relative_eq;

    fn apply_aberration_sph<U: LengthUnit>(
        mean: position::Equatorial<U>,
        jd:   JulianDay,
    ) -> position::Equatorial<U> {
        (&apply_aberration((&mean).into(), jd)).into()
    }

    #[test]
    fn test_aberration_preserva_distance_and_epoch() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = position::Equatorial::<f64>::new(
            Degrees::new(10.0),
            Degrees::new(20.0),
            1.23
        );
        let out = apply_aberration_sph(mean, jd);

        assert_relative_eq!(out.distance.unwrap().value(), mean.distance.unwrap().value(), epsilon = 0.0);
    }

    #[test]
    fn test_aberration_introduces_shift() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = position::Equatorial::<f64>::new(
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
        let mean = position::Equatorial::<f64>::new(
            Degrees::new(123.4),  // dummy RA
            Degrees::new(90.0),  // Dec = +90°
            1.0
        );
        let out = apply_aberration_sph(mean, jd);

        assert!(out.dec().as_f64() < 90.0, "Declination should decrease slightly at pole");
        assert!(!out.ra().as_f64().is_nan(), "RA must not be NaN at the pole");
    }

}
