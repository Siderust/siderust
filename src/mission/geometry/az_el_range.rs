// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Azimuth-elevation-range observation result and solver.

use core::f64::consts::TAU;

use qtty::angular::Radians;
use qtty::length::{Meter, Meters};
use qtty::time::{Second, Seconds};
use qtty::{Per, Quantity};

use super::local_frame::LocalFrame;

/// Velocity expressed as metres per second.
pub type MetersPerSecond = Quantity<Per<Meter, Second>>;

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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AzElRange {
    /// Azimuth measured east of true north in `[0, 2π)`.
    pub azimuth: Radians,
    /// Elevation above the local horizon in `[-π/2, π/2]`.
    pub elevation: Radians,
    /// Slant range from observer to target.
    pub range: Meters,
    /// Line-of-sight range rate (positive = receding).
    pub range_rate: MetersPerSecond,
    /// One-way light time `range / c`.
    pub light_time: Seconds,
    /// Obstruction / terrain-mask diagnostics.
    pub status: LineOfSightStatus,
}

/// Solve [`AzElRange`] from observer / target positions and velocities.
///
/// `observer_pos`, `target_pos` and the velocities are interpreted as
/// vectors in a single shared Cartesian frame, in metres / metres per
/// second. The `frame` supplies the local east/north/up basis used to
/// project the line of sight into azimuth/elevation.
pub fn azimuth_elevation_range(
    observer_pos: [Meters; 3],
    observer_vel: [MetersPerSecond; 3],
    target_pos: [Meters; 3],
    target_vel: [MetersPerSecond; 3],
    frame: LocalFrame,
) -> AzElRange {
    let op = [
        observer_pos[0].value(),
        observer_pos[1].value(),
        observer_pos[2].value(),
    ];
    let tp = [
        target_pos[0].value(),
        target_pos[1].value(),
        target_pos[2].value(),
    ];
    let ov = [
        observer_vel[0].value(),
        observer_vel[1].value(),
        observer_vel[2].value(),
    ];
    let tv = [
        target_vel[0].value(),
        target_vel[1].value(),
        target_vel[2].value(),
    ];

    let los = super::sub(tp, op);
    let rng = super::norm(los);
    let dir = super::unit(los);

    let (east, north, up) = frame.basis();
    let e_comp = super::dot(dir, east);
    let n_comp = super::dot(dir, north);
    let u_comp = super::dot(dir, up);

    let mut az = e_comp.atan2(n_comp);
    if az < 0.0 {
        az += TAU;
    }
    let el = u_comp.asin();

    // Range rate = (v_target - v_observer) · LOS direction.
    let rel_v = super::sub(tv, ov);
    let rng_rate = super::dot(rel_v, dir);

    AzElRange {
        azimuth: Quantity::new(az),
        elevation: Quantity::new(el),
        range: Quantity::new(rng),
        range_rate: Quantity::new(rng_rate),
        light_time: Quantity::new(rng / C_M_PER_S),
        status: LineOfSightStatus::Clear,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use core::f64::consts::PI;

    fn m(x: f64) -> Meters {
        Quantity::new(x)
    }
    fn mps(x: f64) -> MetersPerSecond {
        Quantity::new(x)
    }
    fn rad(x: f64) -> Radians {
        Quantity::new(x)
    }

    #[test]
    fn azel_observer_below_target() {
        // Observer at North pole (lat=π/2, lon=0). Target directly above.
        let r = 6_378_137.0;
        let obs = [m(0.0), m(0.0), m(r)];
        let tgt = [m(0.0), m(0.0), m(r + 1000.0)];
        let zero = [mps(0.0), mps(0.0), mps(0.0)];
        let result = azimuth_elevation_range(
            obs,
            zero,
            tgt,
            zero,
            LocalFrame::new(rad(PI / 2.0), rad(0.0)),
        );
        assert_abs_diff_eq!(result.elevation.value(), PI / 2.0, epsilon = 1e-9);
        assert_abs_diff_eq!(result.range.value(), 1000.0, epsilon = 1e-6);
    }

    #[test]
    fn azel_observer_below_horizon() {
        // Observer at equator, lon=0. Target on opposite side of Earth.
        let r = 6_378_137.0;
        let obs = [m(r), m(0.0), m(0.0)];
        let tgt = [m(-r), m(0.0), m(0.0)];
        let zero = [mps(0.0), mps(0.0), mps(0.0)];
        let result =
            azimuth_elevation_range(obs, zero, tgt, zero, LocalFrame::new(rad(0.0), rad(0.0)));
        // Target is in -up direction → elevation = -π/2.
        assert!(result.elevation.value() < 0.0);
    }
}
