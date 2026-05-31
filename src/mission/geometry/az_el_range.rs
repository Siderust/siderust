// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Azimuth-elevation-range observation result and solver.

use core::f64::consts::TAU;

use affn::cartesian::{
    line_of_sight_with_distance, Direction as CartesianDirection, Position, Velocity,
};
use affn::frames::Horizontal;
use affn::spherical::Direction;
use affn::{ReferenceCenter, ReferenceFrame};
use qtty::length::{Meter, Meters};
use qtty::time::{Second, Seconds};
use qtty::{Per, Quantity};

use super::local_frame::LocalFrame;

/// Scalar metres per second (for range-rate and similar scalar fields).
pub type MetersPerSecond = Quantity<Per<Meter, Second>>;

/// Velocity unit: metres per second.
type MeterPerSecond = Per<Meter, Second>;

/// Speed of light in vacuum, m·s⁻¹. CODATA 2018 exact value.
const C_M_PER_S: f64 = 299_792_458.0;

/// Status of the geometric line of sight between observer and target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineOfSightStatus {
    /// Direct geometric line of sight, no body in between.
    Clear,
    /// A third body (typically the central body) lies between observer
    /// and target.
    BodyObstructed,
    /// The target is below the observer's terrain mask.
    BelowTerrainMask,
}

/// Full mission-geometry observation result for a single epoch.
#[derive(Debug, Clone, Copy)]
pub struct AzElRange {
    /// Local horizontal pointing direction.
    pub pointing: Direction<Horizontal>,
    /// Slant range from observer to target.
    pub range: Meters,
    /// Line-of-sight range rate (positive = receding).
    pub range_rate: MetersPerSecond,
    /// One-way light time `range / c`.
    pub light_time: Seconds,
    /// Obstruction / terrain-mask diagnostics.
    pub status: LineOfSightStatus,
}

impl PartialEq for AzElRange {
    fn eq(&self, other: &Self) -> bool {
        self.pointing.alt() == other.pointing.alt()
            && self.pointing.az() == other.pointing.az()
            && self.range == other.range
            && self.range_rate == other.range_rate
            && self.light_time == other.light_time
            && self.status == other.status
    }
}

/// Solve [`AzElRange`] from observer / target positions and velocities.
///
/// All four vectors must be expressed in the same Cartesian center and
/// frame. The `frame` supplies the local east/north/up basis used to
/// project the line of sight into azimuth/elevation.
pub fn azimuth_elevation_range<C, F>(
    observer_pos: &Position<C, F, Meter>,
    observer_vel: &Velocity<F, MeterPerSecond>,
    target_pos: &Position<C, F, Meter>,
    target_vel: &Velocity<F, MeterPerSecond>,
    frame: LocalFrame,
) -> AzElRange
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let (dir, rng): (CartesianDirection<F>, Meters) =
        line_of_sight_with_distance(observer_pos, target_pos);
    let (east, north, up) = frame.basis::<F>();
    let e_comp = dir.dot(&east);
    let n_comp = dir.dot(&north);
    let u_comp = dir.dot(&up);

    let mut az = e_comp.atan2(n_comp);
    if az < 0.0 {
        az += TAU;
    }
    let el = u_comp.asin();

    // Range rate = (v_target - v_observer) · LOS direction.
    let rel_v = target_vel - observer_vel;
    let rng_rate =
        rel_v.x().value() * dir.x() + rel_v.y().value() * dir.y() + rel_v.z().value() * dir.z();

    AzElRange {
        pointing: Direction::<Horizontal>::new(Quantity::new(el), Quantity::new(az)),
        range: rng,
        range_rate: Quantity::new(rng_rate),
        light_time: Quantity::new(rng.value() / C_M_PER_S),
        status: LineOfSightStatus::Clear,
    }
}

#[cfg(test)]
mod tests {
    use super::LocalFrame;
    use super::*;
    use affn::cartesian::{Position, Velocity};
    use affn::DeriveReferenceCenter;
    use affn::DeriveReferenceFrame;
    use approx::assert_abs_diff_eq;
    use core::f64::consts::PI;
    use qtty::angular::Radians;
    use qtty::length::Meter;
    use qtty::Quantity;

    #[derive(Debug, Clone, Copy, DeriveReferenceCenter)]
    struct C;

    #[derive(Debug, Clone, Copy, DeriveReferenceFrame)]
    struct F;

    type P = Position<C, F, Meter>;
    type V = Velocity<F, MeterPerSecond>;

    fn p(x: f64, y: f64, z: f64) -> P {
        P::new(x, y, z)
    }
    fn v(x: f64, y: f64, z: f64) -> V {
        V::new(x, y, z)
    }
    fn rad(x: f64) -> Radians {
        Quantity::new(x)
    }

    #[test]
    fn azel_observer_below_target() {
        // Observer at North pole (lat=π/2, lon=0). Target directly above.
        let r = 6_378_137.0;
        let obs = p(0.0, 0.0, r);
        let tgt = p(0.0, 0.0, r + 1000.0);
        let zero = v(0.0, 0.0, 0.0);
        let result = azimuth_elevation_range(
            &obs,
            &zero,
            &tgt,
            &zero,
            LocalFrame::new(rad(PI / 2.0), rad(0.0)),
        );
        assert_abs_diff_eq!(result.pointing.alt().value(), PI / 2.0, epsilon = 1e-9);
        assert_abs_diff_eq!(result.pointing.az().value(), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(result.range.value(), 1000.0, epsilon = 1e-6);
    }

    #[test]
    fn azel_observer_below_horizon() {
        // Observer at equator, lon=0. Target on opposite side of Earth.
        let r = 6_378_137.0;
        let obs = p(r, 0.0, 0.0);
        let tgt = p(-r, 0.0, 0.0);
        let zero = v(0.0, 0.0, 0.0);
        let result = azimuth_elevation_range(
            &obs,
            &zero,
            &tgt,
            &zero,
            LocalFrame::new(rad(0.0), rad(0.0)),
        );
        // Target is in -up direction → elevation = -π/2.
        assert!(result.pointing.alt().value() < 0.0);
    }
}
