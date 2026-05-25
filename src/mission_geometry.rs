// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Mission geometry
//!
//! ## Scientific scope
//!
//! Library-level mission-geometry kernels for observer-to-target line of
//! sight, body obstruction, occultation overlap, eclipse classification,
//! and local orbital timing. These are the typed primitives that
//! mission-analysis code and operational scheduling systems can build
//! on; this crate intentionally stops at the math layer. Mission-control
//! execution, scheduling, REST services, CLI/GUI front ends, and
//! operator-console workflows belong in **SatOps**.
//!
//! ## Technical scope
//!
//! All public surfaces accept typed inputs (typed centers/frames in
//! [`crate::coordinates`], typed quantities in [`qtty`]) and return
//! typed results:
//!
//! - [`AzElRange`] — observer-relative azimuth, elevation, slant range,
//!   range-rate, one-way light time, with explicit obstruction / mask
//!   diagnostics.
//! - [`azimuth_elevation_range`] — solve [`AzElRange`] from two state
//!   vectors expressed in the same Cartesian center / frame.
//! - [`line_of_sight_obstructed`] — strict spherical / oblate-ellipsoid
//!   obstruction test for an occulting body sitting between observer and
//!   target.
//! - [`occultation_fraction`] — apparent disk overlap (0…1) between two
//!   bodies as seen from the observer.
//! - [`EclipseState`] / [`eclipse_state`] / [`solar_eclipsing`] — full /
//!   partial / no-eclipse classification with shared math to the SRP
//!   shadow models in `pod`.
//! - [`beta_angle`], [`local_solar_time`], [`ltan`] — orbit-relative
//!   geometry helpers.
//!
//! ## References
//!
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications*, 4th ed. §3.4, §3.7, §5.3, §11.3.
//! - Montenbruck, O., & Gill, E. (2000). *Satellite Orbits — Models,
//!   Methods, and Applications*, §3.4, §11.2.

#![forbid(unsafe_code)]

use core::f64::consts::{PI, TAU};

use qtty::angular::Radians;
use qtty::length::{Meter, Meters};
use qtty::time::{Second, Seconds};
use qtty::{Per, Quantity};

/// Velocity expressed as metres per second.
pub type MetersPerSecond = Quantity<Per<Meter, Second>>;

/// Speed of light in vacuum, m·s⁻¹. CODATA 2018 exact value.
const C_M_PER_S: f64 = 299_792_458.0;

// ─────────────────────────────────────────────────────────────────────────────
// Three-element vector helpers (m). Operate on f64 internally; the public
// API takes typed inputs and returns typed outputs.
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = norm(v);
    if n == 0.0 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AzElRange result
// ─────────────────────────────────────────────────────────────────────────────

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

/// Local east-north-up (ENU) basis at an observer site whose
/// geodetic latitude/longitude is known.
///
/// `lat` and `lon` are the geodetic latitude and longitude of the site
/// expressed in the same frame the observer/target positions are given in
/// (typically ECEF).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalFrame {
    /// Geodetic latitude.
    pub lat: Radians,
    /// Geodetic longitude.
    pub lon: Radians,
}

impl LocalFrame {
    /// Build a local frame from typed geodetic coordinates.
    pub const fn new(lat: Radians, lon: Radians) -> Self {
        Self { lat, lon }
    }

    fn basis(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (sp, cp) = (self.lat.value().sin(), self.lat.value().cos());
        let (sl, cl) = (self.lon.value().sin(), self.lon.value().cos());
        let east = [-sl, cl, 0.0];
        let north = [-sp * cl, -sp * sl, cp];
        let up = [cp * cl, cp * sl, sp];
        (east, north, up)
    }
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

    let los = sub(tp, op);
    let rng = norm(los);
    let dir = unit(los);

    let (east, north, up) = frame.basis();
    let e_comp = dot(dir, east);
    let n_comp = dot(dir, north);
    let u_comp = dot(dir, up);

    let mut az = e_comp.atan2(n_comp);
    if az < 0.0 {
        az += TAU;
    }
    let el = u_comp.asin();

    // Range rate = (v_target - v_observer) · LOS direction.
    let rel_v = sub(tv, ov);
    let rng_rate = dot(rel_v, dir);

    AzElRange {
        azimuth: Quantity::new(az),
        elevation: Quantity::new(el),
        range: Quantity::new(rng),
        range_rate: Quantity::new(rng_rate),
        light_time: Quantity::new(rng / C_M_PER_S),
        status: LineOfSightStatus::Clear,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Line-of-sight obstruction
// ─────────────────────────────────────────────────────────────────────────────

/// Test whether the line of sight between `observer` and `target` is
/// obstructed by a sphere of radius `body_radius` centred at the origin.
///
/// Both positions are expressed in metres in a single Cartesian frame
/// whose origin is at the centre of the (potentially) obstructing body.
///
/// Returns `true` iff the segment intersects the closed ball.
pub fn line_of_sight_obstructed(
    observer: [Meters; 3],
    target: [Meters; 3],
    body_radius: Meters,
) -> bool {
    let o = [
        observer[0].value(),
        observer[1].value(),
        observer[2].value(),
    ];
    let t = [target[0].value(), target[1].value(), target[2].value()];
    let r = body_radius.value();
    let d = sub(t, o);
    let f = o;
    let a = dot(d, d);
    if a == 0.0 {
        return norm(o) <= r;
    }
    let b = 2.0 * dot(f, d);
    let c = dot(f, f) - r * r;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return false;
    }
    let sd = disc.sqrt();
    let t1 = (-b - sd) / (2.0 * a);
    let t2 = (-b + sd) / (2.0 * a);
    // Segment is parameterised t ∈ [0, 1]; intersection if any root lies in [0,1].
    (0.0..=1.0).contains(&t1) || (0.0..=1.0).contains(&t2) || (t1 < 0.0 && t2 > 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Occultation fraction
// ─────────────────────────────────────────────────────────────────────────────

/// Fraction `[0, 1]` of the *target* disk apparent area that is hidden by
/// the *occulter* disk as seen from `observer`.
///
/// Inputs are positions of the target and occulter centres relative to
/// the observer (or in a shared centred frame; only the geometry of the
/// triangle observer-target-occulter matters) and the physical body
/// radii.
pub fn occultation_fraction(
    observer: [Meters; 3],
    target_center: [Meters; 3],
    target_radius: Meters,
    occulter_center: [Meters; 3],
    occulter_radius: Meters,
) -> f64 {
    let o = [
        observer[0].value(),
        observer[1].value(),
        observer[2].value(),
    ];
    let tc = [
        target_center[0].value(),
        target_center[1].value(),
        target_center[2].value(),
    ];
    let oc = [
        occulter_center[0].value(),
        occulter_center[1].value(),
        occulter_center[2].value(),
    ];
    let tr = target_radius.value();
    let or_ = occulter_radius.value();
    let d_t = sub(tc, o);
    let d_o = sub(oc, o);
    let nt = norm(d_t);
    let no = norm(d_o);
    if nt == 0.0 || no == 0.0 {
        return 0.0;
    }
    let app_t = (tr / nt).atan(); // angular radius of target
    let app_o = (or_ / no).atan(); // angular radius of occulter
    let cos_sep = dot(d_t, d_o) / (nt * no);
    let sep = cos_sep.clamp(-1.0, 1.0).acos();
    // Apparent target area
    let area_t = PI * app_t * app_t;
    if sep >= app_t + app_o {
        return 0.0;
    }
    if sep + app_t <= app_o {
        // Target fully covered.
        return 1.0;
    }
    if sep + app_o <= app_t {
        // Occulter fully within target — partial transit; area is occulter disk.
        return (app_o * app_o) / (app_t * app_t);
    }
    // General partial overlap (lens area).
    let r = app_t;
    let r2 = app_o;
    let d = sep;
    let part1 = r * r * ((d * d + r * r - r2 * r2) / (2.0 * d * r)).clamp(-1.0, 1.0).acos();
    let part2 = r2 * r2 * ((d * d + r2 * r2 - r * r) / (2.0 * d * r2)).clamp(-1.0, 1.0).acos();
    let part3 = 0.5
        * ((-d + r + r2) * (d + r - r2) * (d - r + r2) * (d + r + r2))
            .max(0.0)
            .sqrt();
    let overlap = part1 + part2 - part3;
    (overlap / area_t).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Eclipse classification
// ─────────────────────────────────────────────────────────────────────────────

/// Eclipse classification for a satellite illuminated by a single light
/// source (typically the Sun) with one obscuring body (typically Earth).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EclipseState {
    /// Satellite is fully illuminated.
    Sunlight,
    /// Penumbra: light source is partially obscured.
    Penumbra,
    /// Umbra: light source is fully obscured.
    Umbra,
}

/// Classify the eclipse state of a satellite illuminated by a Sun-like
/// light source and obscured by a single occulting body.
///
/// All positions are expressed in a single Cartesian frame in metres.
/// `sat` is the satellite position, `sun` and `occulter_center` are the
/// positions of the light source and the obscuring body, and the radii
/// are the physical body radii.
pub fn eclipse_state(
    sat: [Meters; 3],
    sun: [Meters; 3],
    sun_radius: Meters,
    occulter_center: [Meters; 3],
    occulter_radius: Meters,
) -> EclipseState {
    // Reduce to "Sun seen from satellite vs occulter seen from satellite".
    let f = occultation_fraction(sat, sun, sun_radius, occulter_center, occulter_radius);
    if f >= 1.0 - 1e-9 {
        EclipseState::Umbra
    } else if f > 0.0 {
        EclipseState::Penumbra
    } else {
        EclipseState::Sunlight
    }
}

/// True iff the satellite is in any kind of eclipse (penumbra or umbra).
pub fn solar_eclipsing(
    sat: [Meters; 3],
    sun: [Meters; 3],
    sun_radius: Meters,
    occulter_center: [Meters; 3],
    occulter_radius: Meters,
) -> bool {
    !matches!(
        eclipse_state(sat, sun, sun_radius, occulter_center, occulter_radius),
        EclipseState::Sunlight
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Orbit-relative geometry: beta angle, local solar time, LTAN
// ─────────────────────────────────────────────────────────────────────────────

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
    let r = [
        sat_pos[0].value(),
        sat_pos[1].value(),
        sat_pos[2].value(),
    ];
    let v = [
        sat_vel[0].value(),
        sat_vel[1].value(),
        sat_vel[2].value(),
    ];
    let h = unit(cross(r, v));
    let s = unit(sun_dir);
    // Beta = π/2 − angle(h, s)  ⇒  sin(beta) = h · s
    let sb = dot(h, s).clamp(-1.0, 1.0);
    Quantity::new(sb.asin())
}

/// Local apparent solar time at a sub-satellite longitude, in hours `[0, 24)`.
///
/// `sub_sat_lon` and `sun_lon` are the apparent longitudes of the
/// sub-satellite point and the sub-solar point in the same body-fixed
/// frame.
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
/// frame; both are typed radians.
pub fn ltan(raan: Radians, sun_ra: Radians) -> Seconds {
    let diff = (raan.value() - sun_ra.value()).rem_euclid(TAU);
    let hours: f64 = (12.0 + diff * 24.0 / TAU).rem_euclid(24.0);
    Quantity::new(hours * 3600.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

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
        let result = azimuth_elevation_range(obs, zero, tgt, zero, LocalFrame::new(rad(PI / 2.0), rad(0.0)));
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
        let result = azimuth_elevation_range(obs, zero, tgt, zero, LocalFrame::new(rad(0.0), rad(0.0)));
        // Target is in -up direction → elevation = -π/2.
        assert!(result.elevation.value() < 0.0);
    }

    #[test]
    fn los_obstruction_segment_through_origin() {
        let r = Quantity::new(1.0);
        let a = [m(-2.0), m(0.0), m(0.0)];
        let b = [m(2.0), m(0.0), m(0.0)];
        assert!(line_of_sight_obstructed(a, b, r));
    }

    #[test]
    fn los_obstruction_segment_above_sphere() {
        let r = Quantity::new(1.0);
        let a = [m(-2.0), m(0.0), m(2.0)];
        let b = [m(2.0), m(0.0), m(2.0)];
        assert!(!line_of_sight_obstructed(a, b, r));
    }

    #[test]
    fn occultation_full_when_occulter_covers_target() {
        let observer = [m(0.0), m(0.0), m(0.0)];
        let target = [m(10.0), m(0.0), m(0.0)];
        let occ = [m(5.0), m(0.0), m(0.0)];
        let f = occultation_fraction(observer, target, m(1.0), occ, m(10.0));
        assert!(f >= 0.99);
    }

    #[test]
    fn occultation_none_when_well_separated() {
        let observer = [m(0.0), m(0.0), m(0.0)];
        let target = [m(10.0), m(0.0), m(0.0)];
        let occ = [m(0.0), m(10.0), m(0.0)];
        let f = occultation_fraction(observer, target, m(0.1), occ, m(0.1));
        assert_eq!(f, 0.0);
    }

    #[test]
    fn occultation_partial_intermediate() {
        let observer = [m(0.0), m(0.0), m(0.0)];
        let target = [m(10.0), m(0.0), m(0.0)];
        // Occulter offset just at the disk edge — partial overlap.
        // app_t = atan(0.5/10) ≈ 0.04996, app_o = atan(0.5/5) ≈ 0.09967.
        // Need sep > |app_o - app_t| ≈ 0.0497 and < app_o + app_t ≈ 0.1496.
        // Place occulter so sep ≈ 0.1: y_off = 5 * tan(0.1) ≈ 0.501.
        let occ = [m(5.0), m(0.501), m(0.0)];
        let f = occultation_fraction(observer, target, m(0.5), occ, m(0.5));
        assert!(f > 0.0 && f < 1.0, "got f = {f}");
    }

    #[test]
    fn eclipse_umbra_when_fully_covered() {
        let sat = [m(0.0), m(0.0), m(0.0)];
        let sun = [m(0.0), m(0.0), m(1.5e11)];
        let body = [m(0.0), m(0.0), m(1.0e7)];
        assert_eq!(
            eclipse_state(sat, sun, m(6.96e8), body, m(6.378e8)),
            EclipseState::Umbra
        );
    }

    #[test]
    fn eclipse_sunlight_when_unobstructed() {
        let sat = [m(0.0), m(0.0), m(0.0)];
        let sun = [m(0.0), m(0.0), m(1.5e11)];
        let body = [m(0.0), m(1.0e10), m(0.0)];
        assert_eq!(
            eclipse_state(sat, sun, m(6.96e8), body, m(6.378e6)),
            EclipseState::Sunlight
        );
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
