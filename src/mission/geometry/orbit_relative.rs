// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orbit-relative geometry helpers: beta angle, local solar time, LTAN.

use core::f64::consts::TAU;

use qtty::angular::Radians;
use qtty::length::Meters;
use qtty::time::Seconds;
use qtty::Quantity;

use super::az_el_range::MetersPerSecond;

/// Beta angle: angle between the orbit plane and the satellite-to-Sun
/// vector.
///
/// Returns a signed angle in `[-π/2, π/2]`. Positive values mean the Sun
/// lies on the angular-momentum side of the orbit plane.
pub fn beta_angle(
    sat_pos: [Meters; 3],
    sat_vel: [MetersPerSecond; 3],
    sun_dir: [f64; 3],
) -> Radians {
    let r = [sat_pos[0].value(), sat_pos[1].value(), sat_pos[2].value()];
    let v = [sat_vel[0].value(), sat_vel[1].value(), sat_vel[2].value()];
    let h = super::unit(super::cross(r, v));
    let s = super::unit(sun_dir);
    // Beta = π/2 − angle(h, s)  ⇒  sin(beta) = h · s
    let sb = super::dot(h, s).clamp(-1.0, 1.0);
    Quantity::new(sb.asin())
}

/// Local apparent solar time at a sub-satellite longitude, in hours `[0, 24)`.
///
/// `sub_sat_lon` and `sun_lon` are the apparent longitudes of the
/// sub-satellite point and the sub-solar point in the same body-fixed
/// frame. Returns the local solar time expressed as [`Seconds`].
pub fn local_solar_time(sub_sat_lon: Radians, sun_lon: Radians) -> Seconds {
    let diff = (sub_sat_lon.value() - sun_lon.value()).rem_euclid(TAU);
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
    let diff = (raan.value() - sun_ra.value()).rem_euclid(TAU);
    let hours: f64 = (12.0 + diff * 24.0 / TAU).rem_euclid(24.0);
    Quantity::new(hours * 3600.0)
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
    fn beta_angle_in_plane_is_zero() {
        // Satellite orbit in XY plane, Sun in XY plane → β = 0.
        let r = [m(7e6), m(0.0), m(0.0)];
        let v = [mps(0.0), mps(7500.0), mps(0.0)];
        let s = [1.0, 0.0, 0.0];
        let beta = beta_angle(r, v, s);
        assert_abs_diff_eq!(beta.value(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn beta_angle_perpendicular_is_pi_over_2() {
        let r = [m(7e6), m(0.0), m(0.0)];
        let v = [mps(0.0), mps(7500.0), mps(0.0)];
        let s = [0.0, 0.0, 1.0];
        let beta = beta_angle(r, v, s);
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
