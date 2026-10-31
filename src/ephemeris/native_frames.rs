// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Fixed native dynamical ecliptics → canonical IAU 2006 J2000 ecliptic.
//!
//! The two native realizations have separate published FK5 ties. These are
//! fixed axes: the FK5 catalogue spin is not a rotation of these frames.
//! See `doc/frames_and_centers.md` for provenance and precision limitations.

use crate::coordinates::cartesian::{Position, Velocity};
use crate::coordinates::frames::EclipticMeanJ2000;
use crate::coordinates::transform::frames::bias;
use crate::qtty::{AstronomicalUnit, AstronomicalUnits, Quantity};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar};
use affn::{ReferenceCenter, Rotation3};
use spin::Once;

// IMCCE vsop87.doc, REFERENCE SYSTEM: native VSOP87A/E → FK5 J2000.
const VSOP_TO_FK5: Rotation3 = Rotation3::from_matrix_unchecked([
    [1.0, 0.000_000_440_360, -0.000_000_190_919],
    [-0.000_000_479_966, 0.917_482_137_087, -0.397_776_982_902],
    [0.0, 0.397_776_982_902, 0.917_482_137_087],
]);

// Chapront-Touzé & Chapront, ELP2000-82B notice, §8 p.12:
// native (after Laskar) → FK5 J2000. Do not substitute the VSOP matrix.
const ELP_TO_FK5: Rotation3 = Rotation3::from_matrix_unchecked([
    [1.0, 0.000_000_437_913, -0.000_000_189_859],
    [-0.000_000_477_299, 0.917_482_137_607, -0.397_776_981_701],
    [0.0, 0.397_776_981_701, 0.917_482_137_607],
]);

// ERFA 2.0.1 eraFk5hip r5h (J2000 orientation only), from the published
// Mignard & Froeschlé (2000) orientation (-19.9, -9.1, +22.9) mas.
// This is the FK5 catalogue-to-Hipparcos/ICRS tie, NOT the IAU 2006 bias.
const FK5_TO_ICRS: Rotation3 = Rotation3::from_matrix_unchecked([
    [
        0.999_999_999_999_992_9,
        1.110_223_350_983_549_5e-7,
        4.411_803_963_527_301e-8,
    ],
    [
        -1.110_223_308_498_116_4e-7,
        0.999_999_999_999_989_2,
        -9.647_792_498_531_089e-8,
    ],
    [
        -4.411_805_032_666_219_5e-8,
        9.647_792_009_628_368e-8,
        0.999_999_999_999_994_3,
    ],
]);

pub(crate) fn vsop_to_canonical() -> Rotation3 {
    static ROTATION: Once<Rotation3> = Once::new();
    *ROTATION.call_once(|| bias::icrs_to_ecliptic_j2000() * FK5_TO_ICRS * VSOP_TO_FK5)
}

pub(crate) fn elp_to_canonical() -> Rotation3 {
    static ROTATION: Once<Rotation3> = Once::new();
    *ROTATION.call_once(|| bias::icrs_to_ecliptic_j2000() * FK5_TO_ICRS * ELP_TO_FK5)
}

#[inline]
pub(super) fn vsop_position<C: ReferenceCenter<Params = ()>>(
    x: AstronomicalUnits,
    y: AstronomicalUnits,
    z: AstronomicalUnits,
) -> Position<C, EclipticMeanJ2000, AstronomicalUnit> {
    let [x, y, z] = vsop_to_canonical().apply_array([x.value(), y.value(), z.value()]);
    Position::new(x, y, z)
}

#[inline]
pub(super) fn vsop_velocity(
    x: Quantity<super::AuPerDay>,
    y: Quantity<super::AuPerDay>,
    z: Quantity<super::AuPerDay>,
) -> Velocity<EclipticMeanJ2000, super::AuPerDay> {
    let [x, y, z] = vsop_to_canonical().apply_array([x.value(), y.value(), z.value()]);
    Velocity::new(x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_erfa_basis_vectors_and_handedness() {
        // Independent ERFA 2.0.1: ecm06(J2000) @ fk5hip r5h @ IMCCE.
        let cases = [
            (
                vsop_to_canonical(),
                [
                    [
                        0.999_999_999_999_972,
                        5.268_731_538_978_922e-7,
                        -9.253_316_636_593_184e-8,
                    ],
                    [
                        -5.268_733_925_653_201e-7,
                        1.000_000_000_000_128_3,
                        -7.844_799_937_875_412e-8,
                    ],
                    [
                        9.253_355_243_060_187e-8,
                        7.844_796_404_843_995e-8,
                        1.000_000_000_000_184_3,
                    ],
                ],
            ),
            (
                elp_to_canonical(),
                [
                    [
                        0.999_999_999_999_972_1,
                        5.244_261_537_690_754e-7,
                        -9.147_316_625_277_154e-8,
                    ],
                    [
                        -5.244_264_676_224_863e-7,
                        0.999_999_999_999_489_3,
                        -7.713_925_940_459_517e-8,
                    ],
                    [
                        9.147_268_140_916e-8,
                        7.713_922_423_104e-8,
                        0.999_999_999_999_544_8,
                    ],
                ],
            ),
        ];
        for (rotation, expected) in cases {
            for axis in 0..3 {
                let mut basis = [0.; 3];
                basis[axis] = 1.;
                let out = rotation.apply_array(basis);
                for row in 0..3 {
                    assert!((out[row] - expected[row][axis]).abs() < 3e-16);
                }
                // Published IMCCE matrices have only 12 decimals: their
                // transpose is an inverse to ~1e-12, not machine epsilon.
                let back = rotation.inverse().apply_array(out);
                for row in 0..3 {
                    assert!((back[row] - basis[row]).abs() < 2e-12);
                }
            }
            let m = rotation.as_matrix();
            let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            assert!((det - 1.).abs() < 2e-12);
        }
    }
}
