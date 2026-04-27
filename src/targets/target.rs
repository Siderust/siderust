// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::proper_motion::ProperMotion;
use crate::time::JulianDate;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A coordinate sample tagged with an epoch and optional proper motion.
///
/// `CoordinateWithPM<T>` is intentionally lightweight: it stores one coordinate
/// value of type `T`, the epoch at which that value is valid, and optionally a
/// proper-motion model that can be used by higher-level tracking code.
///
/// This is a sample container, not a full physical object model. Use
/// [`crate::targets::Trackable`] for "things that can produce coordinates at a
/// requested epoch".
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CoordinateWithPM<T> {
    /// Position of the object at epoch `time`.
    pub position: T,
    /// Epoch of the position expressed as a [Julian Day].
    pub time: JulianDate,
    /// Proper‑motion model (e.g. mas yr⁻¹); `None` for objects assumed static.
    pub proper_motion: Option<ProperMotion>,
}

impl<T> CoordinateWithPM<T> {
    /// Creates a coordinate sample with an explicit proper-motion model.
    ///
    /// Use this for catalog entries or observations that are expected to drift
    /// over time.
    pub const fn new(position: T, time: JulianDate, proper_motion: ProperMotion) -> Self {
        CoordinateWithPM {
            position,
            time,
            proper_motion: Some(proper_motion),
        }
    }

    /// Creates a coordinate sample with no proper motion.
    ///
    /// This is appropriate for targets that are treated as fixed over the
    /// caller's time horizon.
    pub const fn new_static(position: T, time: JulianDate) -> Self {
        CoordinateWithPM {
            position,
            time,
            proper_motion: None,
        }
    }

    /// Creates a sample from raw fields.
    ///
    /// This is the most direct constructor and is mainly useful when the caller
    /// already has an `Option<ProperMotion>`.
    pub const fn new_raw(
        position: T,
        time: JulianDate,
        proper_motion: Option<ProperMotion>,
    ) -> Self {
        CoordinateWithPM {
            position,
            time,
            proper_motion,
        }
    }

    /// Returns the stored coordinate sample.
    pub const fn get_position(&self) -> &T {
        &self.position
    }

    /// Returns the stored proper-motion model, if any.
    pub const fn get_proper_motion(&self) -> Option<&ProperMotion> {
        self.proper_motion.as_ref()
    }

    /// Returns the epoch associated with the stored coordinate.
    pub const fn get_time(&self) -> &JulianDate {
        &self.time
    }

    /// Overwrite the position *and* epoch in one operation.
    ///
    /// This helper is convenient when a new observation arrives.  It leaves the
    /// proper‑motion model untouched so that later extrapolations still work.
    pub fn update(&mut self, position: T, time: JulianDate) {
        self.time = time;
        self.position = position;
    }
}

/// Backward-compatibility type alias, existing code that says `Target<T>`
/// continues to compile without changes.  New code should prefer
/// [`CoordinateWithPM`].
pub type Target<T> = CoordinateWithPM<T>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::proper_motion::{ProperMotion, RaProperMotionConvention};
    use crate::bodies::catalog::ALDEBARAN;
    use crate::coordinates::spherical::position::GCRS;
    use crate::time::JulianDate;
    use crate::qtty::*;

    type MilliArcsecondPerDay = crate::qtty::Per<crate::qtty::MilliArcsecond, crate::qtty::Day>;
    type MilliArcsecondsPerDay = crate::qtty::Quantity<MilliArcsecondPerDay>;
    type DegreesPerYear =
        crate::qtty::angular_rate::AngularRate<crate::qtty::unit::Degree, crate::qtty::unit::Year>;

    #[test]
    fn test_target_new() {
        let target = Target::new_static(*ALDEBARAN.coordinate.get_position(), JulianDate::J2000);

        assert_eq!(
            target.position.ra(),
            ALDEBARAN.coordinate.get_position().ra()
        );
        assert_eq!(
            target.position.dec(),
            ALDEBARAN.coordinate.get_position().dec()
        );
        assert_eq!(target.time, JulianDate::J2000);
    }

    #[test]
    fn test_target_new_with_proper_motion() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(45.0), crate::qtty::Degrees::new(30.0), 100.0);
        let proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(10.0),
            MilliArcsecondsPerDay::new(5.0),
        );
        let target = Target::new(position, JulianDate::J2000, proper_motion);

        assert_eq!(target.position.ra().value(), 45.0);
        assert_eq!(target.position.dec().value(), 30.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_some());
    }

    #[test]
    fn test_target_new_static() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(60.0), crate::qtty::Degrees::new(45.0), 200.0);
        let target = Target::new_static(position, JulianDate::J2000);

        assert_eq!(target.position.ra().value(), 60.0);
        assert_eq!(target.position.dec().value(), 45.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_new_raw() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(90.0), crate::qtty::Degrees::new(60.0), 300.0);
        let proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(15.0),
            MilliArcsecondsPerDay::new(8.0),
        );

        // Test with Some(proper_motion)
        let target = Target::new_raw(position, JulianDate::J2000, Some(proper_motion));
        assert_eq!(target.position.ra().value(), 90.0);
        assert_eq!(target.position.dec().value(), 60.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_some());

        // Test with None proper_motion
        let target = Target::new_raw(position, JulianDate::J2000, None);
        assert_eq!(target.position.ra().value(), 90.0);
        assert_eq!(target.position.dec().value(), 60.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_get_position() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(120.0), crate::qtty::Degrees::new(75.0), 400.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let retrieved_position = target.get_position();
        assert_eq!(retrieved_position.ra(), Degrees::new(120.0));
        assert_eq!(retrieved_position.dec(), Degrees::new(75.0));
        assert_eq!(retrieved_position.distance, AstronomicalUnits::new(400.0));
    }

    #[test]
    fn test_target_get_proper_motion() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(150.0), crate::qtty::Degrees::new(80.0), 500.0);
        let proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(20.0),
            MilliArcsecondsPerDay::new(12.0),
        );

        // Test with proper motion
        let target = Target::new(position, JulianDate::J2000, proper_motion);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_some());
        if let Some(pm) = retrieved_pm {
            assert_eq!(pm.pm_ra, DegreesPerYear::new(0.0020291249999999997));
            assert_eq!(pm.pm_dec, DegreesPerYear::new(0.001217475));
            assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
        }

        // Test without proper motion
        let target = Target::new_static(position, JulianDate::J2000);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_none());
    }

    #[test]
    fn test_target_get_time() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(180.0), crate::qtty::Degrees::new(85.0), 600.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let retrieved_time = target.get_time();
        assert_eq!(*retrieved_time, JulianDate::J2000);
    }

    #[test]
    fn test_target_update() {
        let initial_position =
            GCRS::<Au>::new(crate::qtty::Degrees::new(200.0), crate::qtty::Degrees::new(90.0), 700.0);
        let proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(25.0),
            MilliArcsecondsPerDay::new(15.0),
        );
        let mut target = Target::new(initial_position, JulianDate::J2000, proper_motion);

        // Update position and time
        let new_position =
            GCRS::<Au>::new(crate::qtty::Degrees::new(220.0), crate::qtty::Degrees::new(85.0), 800.0);
        let new_time = JulianDate::J2000 + crate::qtty::Days::new(365.25);

        target.update(new_position, new_time);

        // Check that position and time were updated
        assert_eq!(target.position.ra(), Degrees::new(220.0));
        assert_eq!(target.position.dec(), Degrees::new(85.0));
        assert_eq!(target.position.distance, AstronomicalUnits::new(800.0));
        assert_eq!(target.time, new_time);

        // Check that proper motion was preserved
        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert_eq!(pm.pm_ra, DegreesPerYear::new(0.0025364062499999996));
            assert_eq!(pm.pm_dec, DegreesPerYear::new(0.0015218437499999998));
            assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
        }
    }

    #[test]
    fn test_target_debug() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(240.0), crate::qtty::Degrees::new(80.0), 900.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let debug_str = format!("{:?}", target);
        assert!(debug_str.contains("CoordinateWithPM"));
    }

    #[test]
    fn test_target_clone() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(260.0), crate::qtty::Degrees::new(75.0), 1000.0);
        let proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(30.0),
            MilliArcsecondsPerDay::new(18.0),
        );
        let target1 = Target::new(position, JulianDate::J2000, proper_motion);

        let target2 = target1.clone();

        // Check that all fields were cloned correctly
        assert_eq!(target1.position.ra(), target2.position.ra());
        assert_eq!(
            target1.position.dec().value(),
            target2.position.dec().value()
        );
        assert_eq!(target1.position.distance, target2.position.distance);
        assert_eq!(target1.time, target2.time);
        assert_eq!(
            target1.proper_motion.is_some(),
            target2.proper_motion.is_some()
        );
    }

    #[test]
    fn test_target_edge_cases() {
        // Test with zero coordinates
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(0.0), crate::qtty::Degrees::new(0.0), 0.0);
        let target = Target::new_static(position, JulianDate::J2000);
        assert_eq!(target.position.ra(), Degrees::new(0.0));
        assert_eq!(target.position.dec(), Degrees::new(0.0));
        assert_eq!(target.position.distance, AstronomicalUnits::new(0.0));

        // Test with very large coordinates
        let position =
            GCRS::<Au>::new(crate::qtty::Degrees::new(359.999), crate::qtty::Degrees::new(89.999), 1e6);
        let target = Target::new_static(position, JulianDate::J2000);
        assert!((target.position.ra().value() - 359.999).abs() < 1e-6);
        assert!((target.position.dec().value() - 89.999).abs() < 1e-6);
        assert_eq!(target.position.distance, AstronomicalUnits::new(1e6));
    }

    #[test]
    fn test_target_zero_proper_motion() {
        let position = GCRS::<Au>::new(crate::qtty::Degrees::new(280.0), crate::qtty::Degrees::new(70.0), 1100.0);
        let zero_proper_motion = ProperMotion::from_mu_alpha_star::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(0.0),
            MilliArcsecondsPerDay::new(0.0),
        );
        let target = Target::new(position, JulianDate::J2000, zero_proper_motion);

        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert_eq!(pm.pm_ra, DegreesPerYear::new(0.0));
            assert_eq!(pm.pm_dec, DegreesPerYear::new(0.0));
            assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
        }
    }
}
