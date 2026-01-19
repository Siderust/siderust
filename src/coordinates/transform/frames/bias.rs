use affn::Rotation3;

/// Frame bias rotation from ICRS to mean equator/equinox of J2000.0.
///
/// Values from IERS Conventions (2003), expressed as a small rotation matrix.
const FRAME_BIAS_ICRS_TO_J2000: Rotation3 = Rotation3::from_matrix([
    [
        0.999_999_999_999_994_2,
        0.000_000_070_782_794_8,
        -0.000_000_080_562_171_5,
    ],
    [
        -0.000_000_070_782_797_4,
        0.999_999_999_999_996_9,
        -0.000_000_033_060_408_8,
    ],
    [
        0.000_000_080_562_169_6,
        0.000_000_033_060_414_5,
        0.999_999_999_999_993_2,
    ],
]);

#[inline]
pub(super) fn frame_bias_icrs_to_j2000() -> Rotation3 {
    FRAME_BIAS_ICRS_TO_J2000
}

#[inline]
pub(super) fn frame_bias_j2000_to_icrs() -> Rotation3 {
    FRAME_BIAS_ICRS_TO_J2000.inverse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bias_inverse_roundtrip_is_identity() {
        let forward = frame_bias_icrs_to_j2000();
        let inverse = frame_bias_j2000_to_icrs();
        let v = [0.1, -0.2, 0.3];

        let rotated = forward.apply_array(v);
        let back = inverse.apply_array(rotated);

        let eps = 1e-10;
        assert!((back[0] - v[0]).abs() < eps);
        assert!((back[1] - v[1]).abs() < eps);
        assert!((back[2] - v[2]).abs() < eps);
    }

    #[test]
    fn bias_matrix_is_close_to_identity() {
        let rot = frame_bias_icrs_to_j2000();
        let mat = rot.as_matrix();
        let eps = 1e-7;
        for (i, row) in mat.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < eps);
            }
        }
    }
}
