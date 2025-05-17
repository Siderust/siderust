use std::f64::consts::PI;
use super::{Radians, DMS};

pub const DEG: Degrees = Degrees::new(1.0);

/// A strongly-typed representation of an angle in degrees.
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct Degrees(f64);

/// Conversion from degrees to radians using the standard library ratio.
impl From<Degrees> for Radians {
    fn from(deg: Degrees) -> Self {
        deg.to_radians()
    }
}

/// Methods specific to `Degrees`.
impl Degrees {
    /// Creates a new `Degrees` from an `f64`.
    #[inline]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner `f64` value for this angle in degrees.
    pub const fn from_hms(hours: i32, minutes: u32, seconds: f64) -> Degrees {
        let h_deg = hours as f64 * 15.0;
        let m_deg = minutes as f64 * 15.0 / 60.0;
        let s_deg = seconds * 15.0 / 3600.0;
        Self::new(h_deg + m_deg + s_deg)
    }

    /// Returns the inner `f64` value for this angle in degrees.
    #[inline]
    pub const fn as_f64(&self) -> f64 {
        self.0
    }

    /// Convert degrees to radians (equivalent to `Radians::from(*self)`).
    #[inline]
    pub const fn to_radians(self) -> Radians {
        Radians::new(self.as_f64() * PI / 180.0)
    }

    /// Compute the sine of the angle (in degrees), by converting internally to radians.
    #[inline]
    pub fn sin(self) -> f64 {
        self.to_radians().as_f64().sin()
    }

    /// Compute the cosine of the angle (in degrees).
    #[inline]
    pub fn cos(self) -> f64 {
        self.to_radians().as_f64().cos()
    }

    /// Compute the tangent of the angle (in degrees).
    #[inline]
    pub fn tan(self) -> f64 {
        self.to_radians().as_f64().tan()
    }

    /// Normalize an angle in degrees to the range [0, 360).
    #[inline]
    pub fn normalize(self) -> Self {
        Self(self.as_f64().rem_euclid(360.0))
    }

    /// Take the absolute value of the angle.
    #[inline]
    pub fn abs(self) -> Self {
        Self(self.as_f64().abs())
    }

    /// Creates a `Degrees` from a number of arcseconds.
    #[inline]
    pub fn arcsec_to_deg(arcsec: f64) -> Self {
        Self(arcsec / 3600.0)
    }

    /// Normalize to the range [−90, +90] in a symmetrical fashion.
    ///
    /// (Implementation depends on your use case—this is just an example.)
    #[inline]
    pub fn normalize_to_90_range(self) -> Self {
        let y = (self.as_f64() + 90.0).rem_euclid(360.0);
        Self(90.0 - (y - 180.0).abs())
    }
    
    /// Normalize to the range [−180, +180].
    #[inline]
    pub fn normalize_to_180_range(self) -> Self {
        Self((self.as_f64() + 180.0).rem_euclid(360.0) - 180.0)
    }

    pub fn diff_deg(self, other: Degrees) -> Self {
        // diferencia normalizada al intervalo 0 … 360
        let d = (self.as_f64() - other.as_f64()).rem_euclid(360.0);
        if d > 180.0 { Self(360.0 - d) } else { Self(d) }
    }

    pub fn to_dms(self) -> super::DMS {
        DMS::from_degrees(self)
    }
}

/// Implement `Display` for `Degrees`.
impl std::fmt::Display for Degrees {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.6}°", self.as_f64())
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(Degrees);