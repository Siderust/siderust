pub trait Kind {
    #[inline(always)]
    fn validate(_x: f64, _y: f64, _z: f64) {}
}

#[derive(Debug, Copy, Clone)]
pub struct PositionKind;
impl Kind for PositionKind {
    // TODO: Check Distance > 0.0
}

#[derive(Debug, Copy, Clone)]
pub struct DirectionKind;
impl Kind for DirectionKind {
    #[inline(always)]
    fn validate(x: f64, y: f64, z: f64) {
        debug_assert!(
            ((x * x + y * y + z * z) - 1.0).abs() < 1.0e-12,
            "Vector must be unitary, got ({}, {}, {})", x, y, z
        );
    }
}
