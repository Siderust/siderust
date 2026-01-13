use affn::Rotation3;

/// Frame bias rotation from ICRS to mean equator/equinox of J2000.0.
///
/// Values from IERS Conventions (2003), expressed as a small rotation matrix.
const FRAME_BIAS_ICRS_TO_J2000: Rotation3 = Rotation3::from_matrix([
    [0.999_999_999_999_994_2, 0.000_000_070_782_794_8, -0.000_000_080_562_171_5],
    [-0.000_000_070_782_797_4, 0.999_999_999_999_996_9, -0.000_000_033_060_408_8],
    [0.000_000_080_562_169_6, 0.000_000_033_060_414_5, 0.999_999_999_999_993_2],
]);

#[inline]
pub(super) fn frame_bias_icrs_to_j2000() -> Rotation3 {
    FRAME_BIAS_ICRS_TO_J2000
}

#[inline]
pub(super) fn frame_bias_j2000_to_icrs() -> Rotation3 {
    FRAME_BIAS_ICRS_TO_J2000.inverse()
}
