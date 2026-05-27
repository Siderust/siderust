// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Line-of-sight obstruction and occultation-fraction geometry.

use core::f64::consts::PI;

use qtty::length::Meters;

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
    let d = super::sub(t, o);
    let f = o;
    let a = super::dot(d, d);
    if a == 0.0 {
        return super::norm(o) <= r;
    }
    let b = 2.0 * super::dot(f, d);
    let c = super::dot(f, f) - r * r;
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
    let d_t = super::sub(tc, o);
    let d_o = super::sub(oc, o);
    let nt = super::norm(d_t);
    let no = super::norm(d_o);
    if nt == 0.0 || no == 0.0 {
        return 0.0;
    }
    let app_t = (tr / nt).atan(); // angular radius of target
    let app_o = (or_ / no).atan(); // angular radius of occulter
    let cos_sep = super::dot(d_t, d_o) / (nt * no);
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
    let part1 = r
        * r
        * ((d * d + r * r - r2 * r2) / (2.0 * d * r))
            .clamp(-1.0, 1.0)
            .acos();
    let part2 = r2
        * r2
        * ((d * d + r2 * r2 - r * r) / (2.0 * d * r2))
            .clamp(-1.0, 1.0)
            .acos();
    let part3 = 0.5
        * ((-d + r + r2) * (d + r - r2) * (d - r + r2) * (d + r + r2))
            .max(0.0)
            .sqrt();
    let overlap = part1 + part2 - part3;
    (overlap / area_t).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Quantity;

    fn m(x: f64) -> Meters {
        Quantity::new(x)
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
}
