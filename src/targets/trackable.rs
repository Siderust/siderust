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
//!
//! The implementations for `Star` and solar-system bodies live in the `bodies`
//! module to keep `targets` free of `bodies` dependencies.

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
// Implementation: Star  (extracts RA/Dec from stored coordinate)
// =============================================================================
// NOTE: impl Trackable for Star<'_> lives in bodies/stars.rs to avoid a
// bodies ↔ targets circular dependency.

// =============================================================================
// Implementations: Solar-system unit types  (VSOP87 / ELP2000)
// =============================================================================
// NOTE: impl_trackable_vsop87! and impl Trackable for Moon live in
// bodies/solar_system.rs for the same reason.

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::JulianDate;

    #[test]
    fn icrs_direction_is_time_invariant() {
        let dir = direction::ICRS::new(
            crate::qtty::Degrees::new(101.287),
            crate::qtty::Degrees::new(-16.716),
        );
        let at_j2000 = dir.track(JulianDate::J2000);
        let at_j2050 = dir.track(crate::time::jd(
            JulianDate::J2000.raw() + crate::qtty::Days::new(365.25 * 50.0),
        ));
        assert_eq!(at_j2000.ra(), at_j2050.ra());
        assert_eq!(at_j2000.dec(), at_j2050.dec());
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
        let result = sample.track(crate::time::jd(
            JulianDate::J2000.raw() + crate::qtty::Days::new(365.25),
        ));
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
