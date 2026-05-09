// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Trackable, Zero-Cost Object Abstraction
//!
//! The [`Trackable`] trait represents "anything that can produce coordinates
//! at an arbitrary time *t*".  It separates the *object* (a star, a planet,
//! a fixed direction …) from the *coordinate sample* ([`CoordinateWithPM`]).
//!
//! ## Design
//!
//! * **Associated type** `Coords`, each implementor chooses its natural output
//!   type, so the trait is zero-cost through monomorphization.
//! * **No `&mut self`**, tracking is a pure, stateless query.
//! * **Static dispatch**, callers write `fn foo<T: Trackable>(t: &T)` and pay
//!   no dynamic-dispatch overhead.
//!
//! ## Implementors
//!
//! | Type | `Coords` | Backing engine |
//! |------|----------|----------------|
//! | `direction::ICRS` | `direction::ICRS` | identity (fixed direction) |
//! | `CoordinateWithPM<T>` | `T` | identity (returns stored position) |
//! | `Star` | `direction::ICRS` | extracts RA/Dec from stored coordinate |
//! | `Sun`, `Mercury` … `Neptune` | `CoordinateWithPM<Position<Bary, Ecl, AU>>` | VSOP87 |
//! | `Moon` | `Position<Geocentric, Ecl, Km>` | ELP2000 |

use crate::coordinates::spherical::direction;
use crate::time::JulianDate;

use super::CoordinateWithPM;

// =============================================================================
// Trait Definition
// =============================================================================

/// Anything that can produce coordinates at an arbitrary epoch.
///
/// This is the central abstraction that separates *objects* (which evolve in
/// time) from *coordinate samples* (a single snapshot).
///
/// # Example
///
/// ```rust
/// use siderust::targets::Trackable;
/// use siderust::bodies::solar_system::Mars;
/// use siderust::time::JulianDate;
///
/// let pos = Mars.track(JulianDate::J2000);
/// println!("Mars at J2000: {:?}", pos);
/// ```
pub trait Trackable {
    /// The coordinate type produced by this object.
    type Coords;

    /// Produce the object's coordinates at Julian Date `jd`.
    fn track(&self, jd: JulianDate) -> Self::Coords;
}

// =============================================================================
// Implementation: direction::ICRS  (fixed direction, identity)
// =============================================================================

impl Trackable for direction::ICRS {
    type Coords = direction::ICRS;

    #[inline]
    fn track(&self, _jd: JulianDate) -> Self::Coords {
        *self
    }
}

// =============================================================================
// Implementation: CoordinateWithPM<T>  (identity, returns stored position)
// =============================================================================

impl<T: Clone> Trackable for CoordinateWithPM<T> {
    type Coords = T;

    #[inline]
    fn track(&self, _jd: JulianDate) -> Self::Coords {
        self.position.clone()
    }
}

// =============================================================================
// Tests
// ============================================================================= 

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::catalog;
    use crate::bodies::solar_system;
    use crate::qtty::*;
    use crate::time::JulianDate;

    #[test]
    fn icrs_direction_is_time_invariant() {
        let dir = direction::ICRS::new(
            crate::qtty::Degrees::new(101.287),
            crate::qtty::Degrees::new(-16.716),
        );
        let at_j2000 = dir.track(JulianDate::J2000);
        let at_j2050 = dir.track(JulianDate::J2000 + crate::qtty::Days::new(365.25 * 50.0));
        assert_eq!(at_j2000.ra(), at_j2050.ra());
        assert_eq!(at_j2000.dec(), at_j2050.dec());
    }

    #[test]
    fn star_produces_icrs_direction() {
        let sirius = &catalog::SIRIUS;
        let dir = sirius.track(JulianDate::J2000);
        // Sirius: RA ≈ 101.287°, Dec ≈ −16.716°
        let ra_diff = dir.ra() - Degrees::new(101.287);
        assert!(ra_diff <= Degrees::new(0.01) && ra_diff >= -Degrees::new(0.01));
        let dec_diff = dir.dec() + Degrees::new(16.716);
        assert!(dec_diff <= Degrees::new(0.01) && dec_diff >= -Degrees::new(0.01));
    }

    #[test]
    fn sun_produces_barycentric_position() {
        let pos = solar_system::Sun.track(JulianDate::J2000);
        // Sun's barycentric distance from SSB is small (< 0.02 AU)
        let dist = pos.position.distance();
        assert!(
            dist < AstronomicalUnits::new(0.02),
            "Sun should be near SSB, got {}",
            dist
        );
    }

    #[test]
    fn earth_changes_with_time() {
        let p1 = solar_system::Earth.track(JulianDate::J2000);
        let p2 = solar_system::Earth.track(JulianDate::J2000 + crate::qtty::Days::new(182.625));
        // After half a year, Earth should be on the opposite side (~2 AU apart)
        let sep = p1.position.distance_to(&p2.position);
        assert!(
            sep > AstronomicalUnits::new(1.0),
            "Half-year separation should be > 1 AU, got {}",
            sep
        );
    }

    #[test]
    fn moon_produces_geocentric_position() {
        let pos = solar_system::Moon.track(JulianDate::J2000);
        // Moon is ~384 400 km from Earth center
        let dist = pos.distance();
        assert!(
            dist >= Kilometers::new(350_000.0) && dist <= Kilometers::new(410_000.0),
            "Moon distance should be ~384 400 km, got {}",
            dist
        );
    }

    #[test]
    fn planets_produce_nonzero_positions() {
        let jd = JulianDate::J2000;
        let mars = solar_system::Mars.track(jd);
        let dist = mars.position.distance();
        assert!(
            dist > AstronomicalUnits::new(1.0),
            "Mars should be > 1 AU from SSB, got {}",
            dist
        );
    }

    #[test]
    fn coordinate_with_pm_returns_stored_position() {
        use crate::coordinates::spherical::position;
        let pos = position::EquatorialMeanJ2000::<crate::qtty::LightYear>::new(
            crate::qtty::Degrees::new(88.0),
            crate::qtty::Degrees::new(7.0),
            crate::qtty::LightYears::new(548.0),
        );
        let sample = CoordinateWithPM::new_static(pos, JulianDate::J2000);
        let result = sample.track(JulianDate::J2000 + crate::qtty::Days::new(365.25));
        assert_eq!(result.ra(), crate::qtty::Degrees::new(88.0));
        assert_eq!(result.dec(), crate::qtty::Degrees::new(7.0));
    }

    #[test]
    fn generic_trackable_function() {
        fn track_anything<T: Trackable>(obj: &T, jd: JulianDate) -> T::Coords {
            obj.track(jd)
        }

        let dir = direction::ICRS::new(
            crate::qtty::Degrees::new(45.0),
            crate::qtty::Degrees::new(30.0),
        );
        let result = track_anything(&dir, JulianDate::J2000);
        assert_eq!(result.ra(), crate::qtty::Degrees::new(45.0));
        assert_eq!(result.dec(), crate::qtty::Degrees::new(30.0));
    }
}
