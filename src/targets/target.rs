use crate::units::JulianDay;
use crate::astro::proper_motion::ProperMotion;

#[derive(Debug, Clone)]
pub struct Target<T> {
    /// Position of the object at epoch `time`.
    pub position: T,
    /// Epoch of the position expressed as a [Julian Day].
    pub time: JulianDay,
    /// Proper‑motion model (e.g. mas yr⁻¹); `None` for objects assumed static.
    pub proper_motion: Option<ProperMotion>,
}

impl<T> Target<T> {
    pub const fn new(position: T, time: JulianDay, proper_motion: ProperMotion) -> Self {
        Target { position, time, proper_motion: Some(proper_motion) }
    }

    pub const fn new_static(position: T, time: JulianDay) -> Self {
        Target { position, time, proper_motion: None }
    }

    pub const fn new_raw(position: T, time: JulianDay, proper_motion: Option<ProperMotion>) -> Self {
        Target { position, time, proper_motion }
    }

    pub const fn get_position(&self) -> &T {
        &self.position
    }

    pub const fn get_proper_motion(&self) -> Option<&ProperMotion> {
        self.proper_motion.as_ref()
    }

    pub const fn get_time(&self) -> &JulianDay {
        &self.time
    }


    /// Overwrite the position *and* epoch in one operation.
    ///
    /// This helper is convenient when a new observation arrives.  It leaves the
    /// proper‑motion model untouched so that later extrapolations still work.
    pub  fn update(&mut self, position: T, time: JulianDay) {
        self.time = time;
        self.position = position;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::JulianDay;
    use crate::bodies::catalog::ALDEBARAN;
    use crate::astro::proper_motion::ProperMotion;
    use crate::coordinates::spherical::ICRSGeocentricSphericalCoord;

    #[test]
    fn test_target_new() {
        let target = Target::new_static(ALDEBARAN.target.get_position().clone(), JulianDay::J2000);

        assert_eq!(target.position.ra(), ALDEBARAN.target.get_position().ra());
        assert_eq!(target.position.dec(), ALDEBARAN.target.get_position().dec());
        assert_eq!(target.time, JulianDay::J2000);
    }

    #[test]
    fn test_target_new_with_proper_motion() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(45.0),
            crate::units::Degrees::new(30.0),
            100.0
        );
        let proper_motion = ProperMotion::from_mas_per_year(10.0, 5.0);
        let target = Target::new(position.clone(), JulianDay::J2000, proper_motion);

        assert_eq!(target.position.ra().as_f64(), 45.0);
        assert_eq!(target.position.dec().as_f64(), 30.0);
        assert_eq!(target.time, JulianDay::J2000);
        assert!(target.proper_motion.is_some());
    }

    #[test]
    fn test_target_new_static() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(60.0),
            crate::units::Degrees::new(45.0),
            200.0
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);

        assert_eq!(target.position.ra().as_f64(), 60.0);
        assert_eq!(target.position.dec().as_f64(), 45.0);
        assert_eq!(target.time, JulianDay::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_new_raw() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(90.0),
            crate::units::Degrees::new(60.0),
            300.0
        );
        let proper_motion = ProperMotion::from_mas_per_year(15.0, 8.0);

        // Test with Some(proper_motion)
        let target = Target::new_raw(position.clone(), JulianDay::J2000, Some(proper_motion));
        assert_eq!(target.position.ra().as_f64(), 90.0);
        assert_eq!(target.position.dec().as_f64(), 60.0);
        assert_eq!(target.time, JulianDay::J2000);
        assert!(target.proper_motion.is_some());

        // Test with None proper_motion
        let target = Target::new_raw(position.clone(), JulianDay::J2000, None);
        assert_eq!(target.position.ra().as_f64(), 90.0);
        assert_eq!(target.position.dec().as_f64(), 60.0);
        assert_eq!(target.time, JulianDay::J2000);
        assert!(target.proper_motion.is_none());
    }

    #[test]
    fn test_target_get_position() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(120.0),
            crate::units::Degrees::new(75.0),
            400.0
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);

        let retrieved_position = target.get_position();
        assert_eq!(retrieved_position.ra().as_f64(), 120.0);
        assert_eq!(retrieved_position.dec().as_f64(), 75.0);
        assert_eq!(retrieved_position.distance, Some(400.0));
    }

    #[test]
    fn test_target_get_proper_motion() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(150.0),
            crate::units::Degrees::new(80.0),
            500.0
        );
        let proper_motion = ProperMotion::from_mas_per_year(20.0, 12.0);

        // Test with proper motion
        let target = Target::new(position.clone(), JulianDay::J2000, proper_motion);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_some());
        if let Some(pm) = retrieved_pm {
            assert!((pm.ra_μ.value() - 5.555555555555556e-6).abs() < 1e-15);
            assert!((pm.dec_μ.value() - 3.333333333333333e-6).abs() < 1e-15);
        }

        // Test without proper motion
        let target = Target::new_static(position.clone(), JulianDay::J2000);
        let retrieved_pm = target.get_proper_motion();
        assert!(retrieved_pm.is_none());
    }

    #[test]
    fn test_target_get_time() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(180.0),
            crate::units::Degrees::new(85.0),
            600.0
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);

        let retrieved_time = target.get_time();
        assert_eq!(*retrieved_time, JulianDay::J2000);
    }

    #[test]
    fn test_target_update() {
        let initial_position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(200.0),
            crate::units::Degrees::new(90.0),
            700.0
        );
        let proper_motion = ProperMotion::from_mas_per_year(25.0, 15.0);
        let mut target = Target::new(initial_position.clone(), JulianDay::J2000, proper_motion);

        // Update position and time
        let new_position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(220.0),
            crate::units::Degrees::new(85.0),
            800.0
        );
        let new_time = JulianDay::J2000 + crate::units::Days::new(365.25);

        target.update(new_position.clone(), new_time);

        // Check that position and time were updated
        assert_eq!(target.position.ra().as_f64(), 220.0);
        assert_eq!(target.position.dec().as_f64(), 85.0);
        assert_eq!(target.position.distance, Some(800.0));
        assert_eq!(target.time, new_time);

        // Check that proper motion was preserved
        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert!((pm.ra_μ.value() - 6.944444444444445e-6).abs() < 1e-15);
            assert!((pm.dec_μ.value() - 4.166666666666667e-6).abs() < 1e-15);
        }
    }

    #[test]
    fn test_target_debug() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(240.0),
            crate::units::Degrees::new(80.0),
            900.0
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);

        let debug_str = format!("{:?}", target);
        assert!(debug_str.contains("Target"));
    }

    #[test]
    fn test_target_clone() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(260.0),
            crate::units::Degrees::new(75.0),
            1000.0
        );
        let proper_motion = ProperMotion::from_mas_per_year(30.0, 18.0);
        let target1 = Target::new(position.clone(), JulianDay::J2000, proper_motion);

        let target2 = target1.clone();

        // Check that all fields were cloned correctly
        assert_eq!(target1.position.ra().as_f64(), target2.position.ra().as_f64());
        assert_eq!(target1.position.dec().as_f64(), target2.position.dec().as_f64());
        assert_eq!(target1.position.distance, target2.position.distance);
        assert_eq!(target1.time, target2.time);
        assert_eq!(target1.proper_motion.is_some(), target2.proper_motion.is_some());
    }

    #[test]
    fn test_target_edge_cases() {
        // Test with zero coordinates
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(0.0),
            crate::units::Degrees::new(0.0),
            0.0
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);
        assert_eq!(target.position.ra().as_f64(), 0.0);
        assert_eq!(target.position.dec().as_f64(), 0.0);
        assert_eq!(target.position.distance, Some(0.0));

        // Test with very large coordinates
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(359.999),
            crate::units::Degrees::new(89.999),
            1e6
        );
        let target = Target::new_static(position.clone(), JulianDay::J2000);
        assert!((target.position.ra().as_f64() - 359.999).abs() < 1e-6);
        assert!((target.position.dec().as_f64() - 89.999).abs() < 1e-6);
        assert_eq!(target.position.distance, Some(1e6));
    }

    #[test]
    fn test_target_zero_proper_motion() {
        let position = ICRSGeocentricSphericalCoord::new(
            crate::units::Degrees::new(280.0),
            crate::units::Degrees::new(70.0),
            1100.0
        );
        let zero_proper_motion = ProperMotion::from_mas_per_year(0.0, 0.0);
        let target = Target::new(position.clone(), JulianDay::J2000, zero_proper_motion);

        assert!(target.proper_motion.is_some());
        if let Some(pm) = target.get_proper_motion() {
            assert_eq!(pm.ra_μ.value(), 0.0);
            assert_eq!(pm.dec_μ.value(), 0.0);
        }
    }
}
