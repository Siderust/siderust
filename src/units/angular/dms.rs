use super::{Degrees, Radians};
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
        let mut decimal_degrees = degrees.as_f64();
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
        self.to_degrees().to_radians()
    }

    pub fn as_f64(&self) -> f64 {
        self.to_degrees().as_f64()
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
    fn test_dms_new() {
        let dms = DMS::new(false, 30, 15, 45.5);
        assert_eq!(dms.degrees, 30);
        assert_eq!(dms.minutes, 15);
        assert_eq!(dms.seconds, 45.5);
    }

    #[test]
    fn test_dms_from_degrees() {
        let dms = DMS::from_degrees(Degrees::new(30.2625));
        assert_eq!(dms.degrees, 30);
        assert_eq!(dms.minutes, 15);
        assert!((dms.seconds - 45.0).abs() < 1e-10);
    }

    #[test]
    fn test_dms_to_degrees() {
        let dms = DMS::new(false, 30, 15, 45.0);
        let degrees = dms.to_degrees();
        assert!((degrees.as_f64() - 30.2625).abs() < 1e-10);
    }

    #[test]
    fn test_dms_round_trip() {
        let original = Degrees::new(45.123456);
        let dms = DMS::from_degrees(original);
        let converted = dms.to_degrees();
        assert!((original.as_f64() - converted.as_f64()).abs() < 1e-6);
    }

    #[test]
    fn test_dms_display() {
        let dms = DMS::new(false, 30, 15, 45.5);
        let display = format!("{}", dms);
        assert_eq!(display, "30° 15' 45.50\"");
    }

    #[test]
    fn test_dms_debug() {
        let dms = DMS::new(false, 30, 15, 45.5);
        let debug = format!("{:?}", dms);
        assert!(debug.contains("DMS"));
        assert!(debug.contains("30"));
        assert!(debug.contains("15"));
        assert!(debug.contains("45.5"));
    }

    #[test]
    fn test_dms_clone() {
        let dms1 = DMS::new(false, 30, 15, 45.5);
        let dms2 = dms1.clone();
        assert_eq!(dms1.degrees, dms2.degrees);
        assert_eq!(dms1.minutes, dms2.minutes);
        assert_eq!(dms1.seconds, dms2.seconds);
    }

    #[test]
    fn test_dms_copy() {
        let dms1 = DMS::new(false, 30, 15, 45.5);
        let dms2 = dms1.clone(); // Clone since DMS doesn't implement Copy
        assert_eq!(dms1.degrees, dms2.degrees);
        assert_eq!(dms1.minutes, dms2.minutes);
        assert_eq!(dms1.seconds, dms2.seconds);
    }

    #[test]
    fn test_dms_edge_cases() {
        // Zero values
        let dms = DMS::new(false, 0, 0, 0.0);
        assert_eq!(dms.degrees, 0);
        assert_eq!(dms.minutes, 0);
        assert_eq!(dms.seconds, 0.0);

        // Large values
        let dms = DMS::new(false, 359, 59, 59.999);
        assert_eq!(dms.degrees, 359);
        assert_eq!(dms.minutes, 59);
        assert!((dms.seconds - 59.999).abs() < 1e-10);

        // Negative seconds (should be normalized)
        let dms = DMS::new(false, 30, 15, -45.5);
        assert_eq!(dms.degrees, 30);
        assert_eq!(dms.minutes, 15);
        assert_eq!(dms.seconds, 14.5); // -45.5 % 60 = 14.5
    }

    #[test]
    fn test_dms_from_degrees_edge_cases() {
        // Zero degrees
        let dms = DMS::from_degrees(Degrees::new(0.0));
        assert_eq!(dms.degrees, 0);
        assert_eq!(dms.minutes, 0);
        assert_eq!(dms.seconds, 0.0);

        // Very small angle
        let dms = DMS::from_degrees(Degrees::new(0.001));
        assert_eq!(dms.degrees, 0);
        assert_eq!(dms.minutes, 0);
        assert!((dms.seconds - 3.6).abs() < 1e-10);

        // Large angle
        let dms = DMS::from_degrees(Degrees::new(359.999999));
        assert_eq!(dms.degrees, 359);
        assert_eq!(dms.minutes, 59);
        assert!((dms.seconds - 59.9964).abs() < 1e-4);
    }

    #[test]
    fn test_dms_to_degrees_edge_cases() {
        // Zero
        let dms = DMS::new(false, 0, 0, 0.0);
        let degrees = dms.to_degrees();
        assert_eq!(degrees.as_f64(), 0.0);

        // Full circle
        let dms = DMS::new(false, 360, 0, 0.0);
        let degrees = dms.to_degrees();
        assert_eq!(degrees.as_f64(), 360.0);

        // Negative seconds (normalized)
        let dms = DMS::new(false, 30, 15, -45.0);
        let degrees = dms.to_degrees();
        assert!((degrees.as_f64() - 30.254166666666666).abs() < 1e-10); // -45.0 % 60 = 15.0, so 30°15'15"
    }
}
