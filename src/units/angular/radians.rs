use std::f64::consts::PI;
use std::fmt;
use super::Degrees;

/// A strongly-typed representation of an angle in radians.
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct Radians(f64);

/// Conversion from radians to degrees using the standard library ratio.
impl From<Radians> for Degrees {
    fn from(rad: Radians) -> Self {
        Degrees::new(rad.as_f64() * 180.0 / PI)
    }
}

/// Methods specific to `Radians`.
impl Radians {

    pub const TAU: Radians = Radians::new(std::f64::consts::TAU);

    /// Creates a new `Radians` from an `f64`.
    #[inline]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner `f64` value for this angle in radians.
    #[inline]
    pub const fn as_f64(&self) -> f64 {
        self.0
    }

    /// Convert radians to degrees (equivalent to `Degrees::from(*self)`).
    #[inline]
    pub fn to_degrees(self) -> Degrees {
        self.into()
    }

    /// Compute the cosine of the angle.
    #[inline]
    pub fn cos(self) -> f64 {
        self.as_f64().cos()
    }

    /// Compute the sine of the angle.
    #[inline]
    pub fn sin(self) -> f64 {
        self.as_f64().sin()
    }

    /// Compute the tangent of the angle.
    #[inline]
    pub fn tan(self) -> f64 {
        self.as_f64().tan()
    }

    pub fn atan2(self, other: Self) -> Self {
        Self(self.0.atan2(other.0))
    }

    /// Returns the signum of the angle (i.e. sign of the inner value).
    #[inline]
    pub fn signum(self) -> f64 {
        self.as_f64().signum()
    }

    /// Absolute value of the radians.
    #[inline]
    pub fn abs(self) -> Self {
        Self(self.as_f64().abs())
    }

    /// Simultaneously computes the sine and cosine of the angle.
    #[inline]
    pub fn sin_cos(self) -> (f64, f64) {
        self.as_f64().sin_cos()
    }

    /// Normalize range (-π, π].
    #[inline]
    pub fn wrap_pi(self) -> Self {
        let x = self.as_f64();
        let y = (x + PI).rem_euclid(2.0 * PI) - PI;
        let norm = if y <= -PI { y + 2.0 * PI } else { y };
        Radians(norm)
    }
}

/// Implement `Display` for `Radians`.
impl fmt::Display for Radians {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.6} rad", self.as_f64())
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(Radians);
