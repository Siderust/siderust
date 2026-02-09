// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::astro::proper_motion::ProperMotion;
use crate::time::JulianDate;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Target<T> {
    /// Position of the object at epoch `time`.
    pub position: T,
    /// Epoch of the position expressed as a [Julian Day].
    pub time: JulianDate,
    /// Proper‑motion model (e.g. mas yr⁻¹); `None` for objects assumed static.
    pub proper_motion: Option<ProperMotion>,
}

impl<T> Target<T> {
    pub const fn new(position: T, time: JulianDate, proper_motion: ProperMotion) -> Self {
        Target {
            position,
            time,
            proper_motion: Some(proper_motion),
        }
    }

    pub const fn new_static(position: T, time: JulianDate) -> Self {
        Target {
            position,
            time,
            proper_motion: None,
        }
    }

    pub const fn new_raw(
        position: T,
        time: JulianDate,
        proper_motion: Option<ProperMotion>,
    ) -> Self {
        Target {
            position,
            time,
            proper_motion,
        }
    }

    pub const fn get_position(&self) -> &T {
        &self.position
    }

    pub const fn get_proper_motion(&self) -> Option<&ProperMotion> {
        self.proper_motion.as_ref()
    }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::proper_motion::ProperMotion;
    use crate::bodies::catalog::ALDEBARAN;
    use crate::coordinates::spherical::position::GCRS;
    use crate::time::JulianDate;
    use qtty::*;

    type MilliArcsecondPerDay = qtty::Per<qtty::MilliArcsecond, qtty::Day>;
    type MilliArcsecondsPerDay = qtty::Quantity<MilliArcsecondPerDay>;

    #[test]
    fn test_target_new() {
        let target = Target::new_static(ALDEBARAN.target.get_position().clone(), JulianDate::J2000);

        assert_eq!(target.position.ra(), ALDEBARAN.target.get_position().ra());
        assert_eq!(target.position.dec(), ALDEBARAN.target.get_position().dec());
        assert_eq!(target.time, JulianDate::J2000);
    }

    #[test]
    fn test_target_new_with_proper_motion() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(45.0), qtty::Degrees::new(30.0), 100.0);
        let proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(10.0),
            MilliArcsecondsPerDay::new(5.0),
        );
        let target = Target::new(position, JulianDate::J2000, proper_motion);

        assert_eq!(target.position.ra(), 45.0);
        assert_eq!(target.position.dec(), 30.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_some());
    }

    #[test]
    fn test_target_new_static() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(60.0), qtty::Degrees::new(45.0), 200.0);
        let target = Target::new_static(position, JulianDate::J2000);

        assert_eq!(target.position.ra(), 60.0);
        assert_eq!(target.position.dec(), 45.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_new_raw() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(90.0), qtty::Degrees::new(60.0), 300.0);
        let proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(15.0),
            MilliArcsecondsPerDay::new(8.0),
        );

        // Test with Some(proper_motion)
        let target = Target::new_raw(position.clone(), JulianDate::J2000, Some(proper_motion));
        assert_eq!(target.position.ra(), 90.0);
        assert_eq!(target.position.dec(), 60.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_some());

        // Test with None proper_motion
        let target = Target::new_raw(position, JulianDate::J2000, None);
        assert_eq!(target.position.ra(), 90.0);
        assert_eq!(target.position.dec(), 60.0);
        assert_eq!(target.time, JulianDate::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_get_position() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(120.0), qtty::Degrees::new(75.0), 400.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let retrieved_position = target.get_position();
        assert_eq!(retrieved_position.ra(), 120.0);
        assert_eq!(retrieved_position.dec(), 75.0);
        assert_eq!(retrieved_position.distance(), 400.0);
    }

    #[test]
    fn test_target_get_proper_motion() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(150.0), qtty::Degrees::new(80.0), 500.0);
        let proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(20.0),
            MilliArcsecondsPerDay::new(12.0),
        );

        // Test with proper motion
        let target = Target::new(position.clone(), JulianDate::J2000, proper_motion);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_some());
        if let Some(pm) = retrieved_pm {
            assert_eq!(pm.ra_μ, 0.0020291249999999997);
            assert_eq!(pm.dec_μ, 0.001217475);
        }

        // Test without proper motion
        let target = Target::new_static(position, JulianDate::J2000);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_none());
    }

    #[test]
    fn test_target_get_time() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(180.0), qtty::Degrees::new(85.0), 600.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let retrieved_time = target.get_time();
        assert_eq!(*retrieved_time, JulianDate::J2000);
    }

    #[test]
    fn test_target_update() {
        let initial_position =
            GCRS::<Au>::new(qtty::Degrees::new(200.0), qtty::Degrees::new(90.0), 700.0);
        let proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(25.0),
            MilliArcsecondsPerDay::new(15.0),
        );
        let mut target = Target::new(initial_position, JulianDate::J2000, proper_motion);

        // Update position and time
        let new_position =
            GCRS::<Au>::new(qtty::Degrees::new(220.0), qtty::Degrees::new(85.0), 800.0);
        let new_time = JulianDate::J2000 + qtty::Days::new(365.25);

        target.update(new_position, new_time);

        // Check that position and time were updated
        assert_eq!(target.position.ra(), 220.0);
        assert_eq!(target.position.dec(), 85.0);
        assert_eq!(target.position.distance(), 800.0);
        assert_eq!(target.time, new_time);

        // Check that proper motion was preserved
        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert_eq!(pm.ra_μ, 0.0025364062499999996);
            assert_eq!(pm.dec_μ, 0.0015218437499999998);
        }
    }

    #[test]
    fn test_target_debug() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(240.0), qtty::Degrees::new(80.0), 900.0);
        let target = Target::new_static(position, JulianDate::J2000);

        let debug_str = format!("{:?}", target);
        assert!(debug_str.contains("Target"));
    }

    #[test]
    fn test_target_clone() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(260.0), qtty::Degrees::new(75.0), 1000.0);
        let proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
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
        assert_eq!(target1.position.distance(), target2.position.distance());
        assert_eq!(target1.time, target2.time);
        assert_eq!(
            target1.proper_motion.is_some(),
            target2.proper_motion.is_some()
        );
    }

    #[test]
    fn test_target_edge_cases() {
        // Test with zero coordinates
        let position = GCRS::<Au>::new(qtty::Degrees::new(0.0), qtty::Degrees::new(0.0), 0.0);
        let target = Target::new_static(position, JulianDate::J2000);
        assert_eq!(target.position.ra(), 0.0);
        assert_eq!(target.position.dec(), 0.0);
        assert_eq!(target.position.distance(), 0.0);

        // Test with very large coordinates
        let position =
            GCRS::<Au>::new(qtty::Degrees::new(359.999), qtty::Degrees::new(89.999), 1e6);
        let target = Target::new_static(position, JulianDate::J2000);
        assert!((target.position.ra().value() - 359.999).abs() < 1e-6);
        assert!((target.position.dec().value() - 89.999).abs() < 1e-6);
        assert_eq!(target.position.distance(), 1e6);
    }

    #[test]
    fn test_target_zero_proper_motion() {
        let position = GCRS::<Au>::new(qtty::Degrees::new(280.0), qtty::Degrees::new(70.0), 1100.0);
        let zero_proper_motion = ProperMotion::new::<MilliArcsecondPerDay>(
            MilliArcsecondsPerDay::new(0.0),
            MilliArcsecondsPerDay::new(0.0),
        );
        let target = Target::new(position, JulianDate::J2000, zero_proper_motion);

        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert_eq!(pm.ra_μ, 0.0);
            assert_eq!(pm.dec_μ, 0.0);
        }
    }
}
