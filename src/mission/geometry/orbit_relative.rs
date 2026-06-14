// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orbit-relative geometry helpers: beta angle, local solar time, LTAN.

use core::f64::consts::TAU;

use affn::cartesian::{Direction, Position, Velocity};
use affn::{ReferenceCenter, ReferenceFrame};
use qtty::angular::Radians;
use qtty::length::Meter;
use qtty::time::{Second, Seconds};
use qtty::{Per, Quantity};

/// Velocity unit: metres per second.
type MeterPerSecond = Per<Meter, Second>;

/// Beta angle: angle between the orbit plane and the satellite-to-Sun
/// vector.
///
/// Returns a signed angle in `[-π/2, π/2]`. Positive values mean the Sun
/// lies on the angular-momentum side of the orbit plane.
pub fn beta_angle<C, F>(
    sat_pos: &Position<C, F, Meter>,
    sat_vel: &Velocity<F, MeterPerSecond>,
    sun_dir: &Direction<F>,
) -> Radians
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let r_hat = sat_pos.direction_unchecked();
    let v_hat = Direction::<F>::new(
        sat_vel.x().value(),
        sat_vel.y().value(),
        sat_vel.z().value(),
    );
    let h = r_hat.cross(&v_hat).expect("non-degenerate orbit");
    // Beta = π/2 − angle(h, s)  ⇒  sin(beta) = h · s
    let sb = h.dot(sun_dir).clamp(-1.0, 1.0);
    Quantity::new(sb.asin())
}

/// Local apparent solar time at a sub-satellite longitude, in hours `[0, 24)`.
///
/// `sub_sat_lon` and `sun_lon` are the apparent longitudes of the
/// sub-satellite point and the sub-solar point in the same body-fixed
/// frame. Returns the local solar time expressed as [`Seconds`].
pub fn local_solar_time(sub_sat_lon: Radians, sun_lon: Radians) -> Seconds {
    let diff = (sub_sat_lon - sun_lon).rem_euclid(TAU).value();
    // Solar noon at the sub-solar point. Convert from radians to seconds.
    let hours: f64 = 12.0 + diff * 24.0 / TAU;
    let h = hours.rem_euclid(24.0);
    Quantity::new(h * 3600.0)
}

/// Local time of the ascending node, in seconds since local midnight
/// at the equator.
///
/// `raan` is the right ascension of the ascending node and
/// `sun_ra` is the right ascension of the Sun in the same inertial
/// frame; both are typed radians. Returns the LTAN expressed as
/// [`Seconds`].
pub fn ltan(raan: Radians, sun_ra: Radians) -> Seconds {
    let diff = (raan - sun_ra).rem_euclid(TAU).value();
    let hours: f64 = (12.0 + diff * 24.0 / TAU).rem_euclid(24.0);
    Quantity::new(hours * 3600.0)
}

#[cfg(test)]
mod tests {
    use super::MeterPerSecond;
    use super::*;
    use affn::cartesian::{Direction, Position, Velocity};
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

    fn pos(x: f64, y: f64, z: f64) -> P {
        P::new(x, y, z)
    }
    fn vel(x: f64, y: f64, z: f64) -> V {
        V::new(x, y, z)
    }
    fn rad(x: f64) -> Radians {
        Quantity::new(x)
    }

    #[test]
    fn beta_angle_in_plane_is_zero() {
        // Satellite orbit in XY plane, Sun in XY plane → β = 0.
        let r = pos(7e6, 0.0, 0.0);
        let v = vel(0.0, 7500.0, 0.0);
        let s = Direction::<F>::new(1.0, 0.0, 0.0);
        let beta = beta_angle(&r, &v, &s);
        assert_abs_diff_eq!(beta.value(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn beta_angle_perpendicular_is_pi_over_2() {
        let r = pos(7e6, 0.0, 0.0);
        let v = vel(0.0, 7500.0, 0.0);
        let s = Direction::<F>::new(0.0, 0.0, 1.0);
        let beta = beta_angle(&r, &v, &s);
        assert_abs_diff_eq!(beta.value(), PI / 2.0, epsilon = 1e-12);
    }

    #[test]
    fn local_solar_time_noon_at_subsolar() {
        let t = local_solar_time(rad(1.234), rad(1.234));
        assert_abs_diff_eq!(t.value(), 12.0 * 3600.0, epsilon = 1e-6);
    }

    #[test]
    fn local_solar_time_wraps_around_midnight() {
        // 12 hours offset → local midnight (0 or 24).
        let t = local_solar_time(rad(0.0), rad(PI));
        let h = t.value() / 3600.0;
        // Either 0 or 24 (modular).
        assert!(h < 1e-9 || (h - 24.0).abs() < 1e-9);
    }

    #[test]
    fn ltan_offset_progresses_with_raan() {
        let t1 = ltan(rad(0.0), rad(0.0));
        let t2 = ltan(rad(PI / 2.0), rad(0.0));
        // Quarter of a day later = 6 hours.
        assert_abs_diff_eq!(t2.value() - t1.value(), 6.0 * 3600.0, epsilon = 1e-6);
    }
}
