use super::*;
use std::ops::Mul;

/// A strongly-typed representation of an angle in Degrees, Minutes, and Seconds (DMS),
/// with an explicit sign flag to avoid ambiguity.
#[derive(Debug, Clone, PartialEq)]
pub struct DMS {
    /// Sign of the angle: `true` for negative, `false` for positive or zero.
    pub sign: bool,
    /// Degree component, always non-negative.
    pub degrees: u32,
    /// Minute component, always in 0..60.
    pub minutes: u32,
    /// Second component, always in 0.0..60.0.
    pub seconds: f64,
}

impl DMS {
    pub const POSITIVE: bool = false;
    pub const NEGATIVE: bool = true;

    /// Creates a new `DMS` instance, normalizing all components and extracting the sign.
    pub fn new(sign: bool, degrees: u32, mut minutes: u32, mut seconds: f64) -> Self {
        // Normalize seconds overflow
        minutes += seconds.div_euclid(60.0) as u32;
        seconds = seconds.rem_euclid(60.0);

        DMS{
            sign,
            degrees: degrees + minutes.div_euclid(60),
            minutes: minutes % 60,
            seconds,
        }
    }

    pub fn from_degrees(degrees: Degrees) -> DMS{
        let mut decimal_degrees = degrees.value();
        let sign: bool = decimal_degrees < 0.0;
        decimal_degrees = decimal_degrees.abs();
    
        let degrees = decimal_degrees.floor();
        let minutes_total = (decimal_degrees - degrees) * 60.0;
        let minutes = minutes_total.floor();
        let seconds = (minutes_total - minutes) * 60.0;
    
        DMS{sign, degrees: degrees as u32, minutes: minutes as u32, seconds}
    }

    pub fn from_milliarcseconds(mut milliarcseconds: f64) -> DMS {
        let sign: bool = milliarcseconds < 0.0;
        milliarcseconds = milliarcseconds.abs();
        let arcsec = milliarcseconds / 1000.0;
        let degrees = (arcsec / 3600.0).floor() as u32;
        let minutes = ((arcsec % 3600.0) / 60.0).floor() as u32;
        let seconds = arcsec % 60.0;
        DMS{sign, degrees, minutes, seconds}
    }

    /// Converts the DMS angle to a decimal-degree `f64`, applying the stored sign.
    pub fn to_degrees(&self) -> Degrees {
        let abs_val = (self.degrees as f64)
            + (self.minutes as f64 / 60.0)
            + (self.seconds / 3600.0);
        if self.sign { Degrees::new(-abs_val) } else { Degrees::new(abs_val) }
    }

    /// Converts to a `Radians` type.
    pub fn to_radians(&self) -> Radians {
        self.to_degrees().to::<Radian>()
    }

    pub fn value(&self) -> f64 {
        self.to_degrees().value()
    }
}

/// Implement `Display` for `DMS`.
impl std::fmt::Display for DMS {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = if self.sign { "-" } else { "" };
        write!(f, "{}{}° {}' {:.2}\"", prefix, self.degrees, self.minutes, self.seconds)
    }
}

impl Mul<f64> for DMS {
    type Output = DMS;

    fn mul(self, rhs: f64) -> DMS {
        let result = self.to_degrees() * rhs;
        DMS::from_degrees(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dms_creation_normal() {
        let angle = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
        assert!(angle.sign == DMS::POSITIVE);
        assert_eq!(angle.degrees, 12);
        assert_eq!(angle.minutes, 34);
        assert_eq!(angle.seconds, 56.0);
    }

    #[test]
    fn test_dms_creation_seconds_overflow() {
        let angle = DMS::new(DMS::POSITIVE, 12, 34, 120.0);
        assert!(angle.sign == DMS::POSITIVE);
        assert_eq!(angle.degrees, 12);
        assert_eq!(angle.minutes, 36);
        assert_eq!(angle.seconds, 0.0);
    }

    #[test]
    fn test_dms_creation_minutes_overflow() {
        let angle = DMS::new(DMS::POSITIVE, 12, 120, 30.0);
        assert!(angle.sign == DMS::POSITIVE);
        assert_eq!(angle.degrees, 14);
        assert_eq!(angle.minutes, 0);
        assert_eq!(angle.seconds, 30.0);
    }



    #[test]
    fn test_dms_creation_negative_degrees() {
        let angle = DMS::new(DMS::NEGATIVE, 12, 34, 56.0);
        assert!(angle.sign == DMS::NEGATIVE);
        assert_eq!(angle.degrees, 12);
        assert_eq!(angle.minutes, 34);
        assert_eq!(angle.seconds, 56.0);
    }

    #[test]
    fn test_dms_creation_all_overflow() {
        let angle = DMS::new(DMS::POSITIVE, 12, 120, 120.0);
        assert!(angle.sign == DMS::POSITIVE);
        assert_eq!(angle.degrees, 14);
        assert_eq!(angle.minutes, 2);
        assert_eq!(angle.seconds, 0.0);
    }

    #[test]
    fn test_to_decimal_positive() {
        let angle = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
        let decimal = angle.value();
        assert!((decimal - 12.582222).abs() < 1e-6);
    }

    #[test]
    fn test_to_decimal_negative() {
        let angle = DMS::new(DMS::NEGATIVE, 12, 34, 56.0);
        let decimal = angle.value();
        assert!((decimal + 12.582222).abs() < 1e-6);
    }

    #[test]
    fn test_to_string_positive() {
        let angle = DMS::new(DMS::POSITIVE, 12, 34, 56.0);
        assert_eq!(angle.to_string(), "12° 34' 56.00\"");
    }

    #[test]
    fn test_to_string_negative() {
        let angle = DMS::new(DMS::NEGATIVE, 12, 34, 56.0);
        assert_eq!(angle.to_string(), "-12° 34' 56.00\"");
    }

    #[test]
    fn test_combined_overflow() {
        let angle = DMS::new(DMS::POSITIVE, 0, 0, 3660.0);
        assert!(angle.sign == DMS::POSITIVE);
        assert_eq!(angle.degrees, 1);
        assert_eq!(angle.minutes, 1);
        assert_eq!(angle.seconds, 0.0);
    }

    #[test]
    fn test_dms_mul() {
        let dms = DMS::new(DMS::POSITIVE, 2, 30, 0.0); // 2.5 deg
        let result = dms * -2.0;
        // 2.5 * -2 = -5.0 deg
        assert!(result.sign == DMS::NEGATIVE);
        assert_eq!(result.degrees, 5);
        assert_eq!(result.minutes, 0);
        assert!((result.seconds - 0.0).abs() < 1e-8);
    }
}
