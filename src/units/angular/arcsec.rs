use std::fmt;
use super::Degrees;

/// A strongly-typed representation of an angle in arcseconds.
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct ArcSecond(pub f64);

impl ArcSecond {
    /// Creates a new `ArcSecond` from an `f64`.
    #[inline]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner `f64` value for this angle in arcseconds.
    #[inline]
    pub const fn as_f64(&self) -> f64 {
        self.0
    }

    /// Convert arcseconds to degrees.
    #[inline]
    pub fn to_degrees(self) -> Degrees {
        Degrees::new(self.0 / 3600.0)
    }

    /// Absolute value of the arcseconds.
    #[inline]
    pub fn abs(self) -> Self {
        Self(self.0.abs())
    }
}

/// Conversion from ArcSecond to Degrees.
impl From<ArcSecond> for Degrees {
    fn from(arcsec: ArcSecond) -> Self {
        arcsec.to_degrees()
    }
}

/// Conversion from Degrees to ArcSecond.
impl From<Degrees> for ArcSecond {
    fn from(deg: Degrees) -> Self {
        ArcSecond(deg.as_f64() * 3600.0)
    }
}

/// Implement `Display` for `ArcSecond`.
impl fmt::Display for ArcSecond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.6}\"", self.0)
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(ArcSecond);
