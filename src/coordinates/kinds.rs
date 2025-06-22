pub trait Kind {
    #[inline(always)]
    fn validate(_x: f64, _y: f64, _z: f64) {}
}

#[derive(Debug, Copy, Clone)]
pub struct PositionKind;
impl Kind for PositionKind {
    #[inline(always)]
    fn validate(x: f64, y: f64, z: f64) {
        debug_assert!(
            (x * x + y * y + z * z).sqrt() !=  std::f64::INFINITY,
            "Vector cannot be infinite, got ({}, {}, {})", x, y, z
        );
        debug_assert!(
            (x * x + y * y + z * z).sqrt() !=  std::f64::NAN,
            "Vector cannot be NaN, got ({}, {}, {})", x, y, z
        );
    }
}

#[derive(Debug, Copy, Clone)]
pub struct VelocityKind;
impl Kind for VelocityKind {
    #[inline(always)]
    fn validate(x: f64, y: f64, z: f64) {
        debug_assert!(
            (x * x + y * y + z * z).sqrt() !=  std::f64::INFINITY,
            "Vector cannot be infinite, got ({}, {}, {})", x, y, z
        );
        debug_assert!(
            (x * x + y * y + z * z).sqrt() !=  std::f64::NAN,
            "Vector cannot be NaN, got ({}, {}, {})", x, y, z
        );
    }
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
