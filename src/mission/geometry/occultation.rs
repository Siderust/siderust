// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Line-of-sight obstruction and occultation-fraction geometry.

use core::f64::consts::PI;

use affn::cartesian::{Displacement, Position};
use affn::{ReferenceCenter, ReferenceFrame};
use qtty::length::{Meter, Meters};

/// Test whether the line of sight between `observer` and `target` is
/// obstructed by a sphere of radius `body_radius` centred at the origin of
/// the shared Cartesian center/frame.
///
/// Returns `true` iff the segment intersects the closed ball.
pub fn line_of_sight_obstructed<C, F>(
    observer: &Position<C, F, Meter>,
    target: &Position<C, F, Meter>,
    body_radius: Meters,
) -> bool
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let d: Displacement<F, Meter> = target - observer;
    let f = Displacement::<F, Meter>::new(observer.x(), observer.y(), observer.z());
    let r = body_radius.value();
    let a = d.magnitude_squared().value();
    if a == 0.0 {
        return f.magnitude().value() <= r;
    }
    let b = 2.0 * f.dot(&d).value();
    let c = f.magnitude_squared().value() - r * r;
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
/// All three positions must be expressed in the same Cartesian center and
/// frame. Only the geometry of the observer-target-occulter triangle matters.
pub fn occultation_fraction<C, F>(
    observer: &Position<C, F, Meter>,
    target_center: &Position<C, F, Meter>,
    target_radius: Meters,
    occulter_center: &Position<C, F, Meter>,
    occulter_radius: Meters,
) -> f64
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let d_t: Displacement<F, Meter> = target_center - observer;
    let d_o: Displacement<F, Meter> = occulter_center - observer;
    let tr = target_radius.value();
    let or_ = occulter_radius.value();
    let nt = d_t.magnitude().value();
    let no = d_o.magnitude().value();
    if nt == 0.0 || no == 0.0 {
        return 0.0;
    }
    let app_t = (tr / nt).atan(); // angular radius of target
    let app_o = (or_ / no).atan(); // angular radius of occulter
    let cos_sep = d_t.dot(&d_o).value() / (nt * no);
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
    use affn::cartesian::Position;
    use affn::DeriveReferenceCenter;
    use affn::DeriveReferenceFrame;
    use qtty::length::{Meter, Meters};
    use qtty::Quantity;

    #[derive(Debug, Clone, Copy, DeriveReferenceCenter)]
    struct C;

    #[derive(Debug, Clone, Copy, DeriveReferenceFrame)]
    struct F;

    type P = Position<C, F, Meter>;

    fn p(x: f64, y: f64, z: f64) -> P {
        P::new(x, y, z)
    }
    fn m(x: f64) -> Meters {
        Quantity::<Meter>::new(x)
    }

    #[test]
    fn los_obstruction_segment_through_origin() {
        let r = m(1.0);
        let a = p(-2.0, 0.0, 0.0);
        let b = p(2.0, 0.0, 0.0);
        assert!(line_of_sight_obstructed(&a, &b, r));
    }

    #[test]
    fn los_obstruction_segment_above_sphere() {
        let r = m(1.0);
        let a = p(-2.0, 0.0, 2.0);
        let b = p(2.0, 0.0, 2.0);
        assert!(!line_of_sight_obstructed(&a, &b, r));
    }

    #[test]
    fn occultation_full_when_occulter_covers_target() {
        let observer = p(0.0, 0.0, 0.0);
        let target = p(10.0, 0.0, 0.0);
        let occ = p(5.0, 0.0, 0.0);
        let f = occultation_fraction(&observer, &target, m(1.0), &occ, m(10.0));
        assert!(f >= 0.99);
    }

    #[test]
    fn occultation_none_when_well_separated() {
        let observer = p(0.0, 0.0, 0.0);
        let target = p(10.0, 0.0, 0.0);
        let occ = p(0.0, 10.0, 0.0);
        let f = occultation_fraction(&observer, &target, m(0.1), &occ, m(0.1));
        assert_eq!(f, 0.0);
    }

    #[test]
    fn occultation_partial_intermediate() {
        let observer = p(0.0, 0.0, 0.0);
        let target = p(10.0, 0.0, 0.0);
        // Occulter offset just at the disk edge — partial overlap.
        // app_t = atan(0.5/10) ≈ 0.04996, app_o = atan(0.5/5) ≈ 0.09967.
        // Need sep > |app_o - app_t| ≈ 0.0497 and < app_o + app_t ≈ 0.1496.
        // Place occulter so sep ≈ 0.1: y_off = 5 * tan(0.1) ≈ 0.501.
        let occ = p(5.0, 0.501, 0.0);
        let f = occultation_fraction(&observer, &target, m(0.5), &occ, m(0.5));
        assert!(f > 0.0 && f < 1.0, "got f = {f}");
    }
}
