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

    #[test]
    fn test_target_new() {
        let target = Target::new_static(ALDEBARAN.target.get_position().clone(), JulianDay::J2000);

        assert_eq!(target.position.ra(), ALDEBARAN.target.get_position().ra());
        assert_eq!(target.position.dec(), ALDEBARAN.target.get_position().dec());
        assert_eq!(target.time, JulianDay::J2000);
    }
}
