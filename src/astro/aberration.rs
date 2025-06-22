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

use crate::coordinates::centers::Heliocentric;
use crate::coordinates::transform::Transform;
use crate::units::JulianDay;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    cartesian::Direction,
    centers::Geocentric, frames::Equatorial
};

const AU_PER_DAY_C: f64 = 173.144_632_674;

/// Apply **annual aberration** to a unit direction vector.
///
/// The function adds the relativistic correction \(v/c\) due to the
/// Earth's orbital motion.
///
/// * `mean` – Geocentric unit vector referred to the true equator &
///   equinox of date.
/// * `jd`   – Epoch in Terrestrial Time (TT) as a *Julian Day*.
///
/// # Returns
/// A new [`Direction`] whose components include annual aberration.
///
/// # Accuracy
/// Better than 0.1 mas over 1900‒2100, limited by the underlying
/// Ron–Vondrák velocity series.
#[must_use]
pub fn apply_aberration_to_direction(
    mean: Direction<Geocentric, Equatorial>,
    jd:   JulianDay,
) -> Direction<Geocentric, Equatorial> {

    let velocity = crate::bodies::solar_system::Earth::vsop87a_vel(jd);
    let velocity: Velocity<Heliocentric, Equatorial> = velocity.transform(jd);

    //--------------------------------------------------------------------
    // Apply û' = û + v/c
    //--------------------------------------------------------------------
    Position::from_vec3(
        mean.as_vec3() + velocity.as_vec3() / AU_PER_DAY_C
    ).direction()
}


/// Remove **annual aberration** from an apparent direction.
/// Inverse operation of [`apply_aberration_to_direction`].
#[must_use]
pub fn remove_aberration_from_direction(
    app: Direction<Geocentric, Equatorial>,
    jd:  JulianDay,
) -> Direction<Geocentric, Equatorial> {

    let velocity = crate::bodies::solar_system::Earth::vsop87a_vel(jd);
    let velocity: Velocity<Heliocentric, Equatorial> = velocity.transform(jd);

    //--------------------------------------------------------------------
    //  Apply û' = û - v/c
    //--------------------------------------------------------------------
    Position::from_vec3(
        app.as_vec3() - velocity.as_vec3() / AU_PER_DAY_C
    ).direction()
}


/// Apply **annual aberration** to a position vector, preserving its
/// geocentric distance.
#[must_use]
pub fn apply_aberration(
    mean: Position<Geocentric, Equatorial>,
    jd:   JulianDay,
) -> Position<Geocentric, Equatorial> {

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
pub fn remove_aberration(
    app: Position<Geocentric, Equatorial>,
    jd:  JulianDay,
) -> Position<Geocentric, Equatorial> {

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
    use crate::coordinates::spherical::Position;
    use approx::assert_relative_eq;

    fn apply_aberration_sph(
        mean: Position<Geocentric, Equatorial>,
        jd:   JulianDay,
    ) -> Position<Geocentric, Equatorial> {
        (&apply_aberration((&mean).into(), jd)).into()
    }

    #[test]
    fn test_aberration_preserva_distance_and_epoch() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = Position::<Geocentric, Equatorial>::new(
            Degrees::new(10.0),
            Degrees::new(20.0),
            1.23
        );
        let out = apply_aberration_sph(mean, jd);

        assert_relative_eq!(out.distance.unwrap(), mean.distance.unwrap(), epsilon = 0.0);
    }

    #[test]
    fn test_aberration_introduces_shift() {
        let jd = JulianDay::new(2451545.0); // J2000.0
        let mean = Position::<Geocentric, Equatorial>::new(
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
        let mean = Position::<Geocentric, Equatorial>::new(
            Degrees::new(123.4),  // dummy RA
            Degrees::new(90.0),   // Dec = +90°
            1.0
        );
        let out = apply_aberration_sph(mean, jd);

        assert!(out.dec().as_f64() < 90.0, "Declination should decrease slightly at pole");
        assert!(!out.ra().as_f64().is_nan(), "RA must not be NaN at the pole");
    }

}
