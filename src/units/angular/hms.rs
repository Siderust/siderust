use super::*;

/// A strongly-typed representation of an angle in Hours, Minutes, and Seconds (HMS),
/// with an explicit sign flag to avoid ambiguity.
#[derive(Debug, Clone, PartialEq)]
pub struct HMS {
    /// Sign of the time: `true` for negative, `false for positive or zero.
    pub sign: bool,
    /// Hour component, always non-negative.
    pub hours: u32,
    /// Minute component, always in 0..60.
    pub minutes: u32,
    /// Second component, always in 0.0..60.0.
    pub seconds: f64,
}

impl HMS {
    /// Creates a new `HMS` instance, normalizing all components and extracting the sign.
    pub fn new(mut hours: i32, mut minutes: i32, mut seconds: f64) -> Self {
        // Normalize seconds overflow (positive or negative)
        if seconds >= 60.0 || seconds <= -60.0 {
            let extra_minutes = (seconds / 60.0).trunc() as i32;
            seconds -= extra_minutes as f64 * 60.0;
            minutes += extra_minutes;
        }

        // Handle negative seconds by borrowing from minutes
        if seconds < 0.0 {
            seconds += 60.0;
            minutes -= 1;
        }

        // Normalize minutes overflow (positive or negative)
        if minutes >= 60 || minutes <= -60 {
            let extra_hours = minutes / 60;
            minutes %= 60;
            hours += extra_hours;
        }

        // Handle negative minutes by borrowing from hours
        if minutes < 0 {
            minutes += 60;
            hours -= 1;
        }

        // Determine final sign and make all components non-negative
        let sign = hours < 0;
        let hours_abs = hours.unsigned_abs();
        let minutes_abs = minutes as u32;

        HMS { sign, hours: hours_abs, minutes: minutes_abs, seconds }
    }

    /// Converts the HMS angle to a decimal-hour `f64`, applying the stored sign.
    pub fn to_decimal(&self) -> f64 {
        let abs_val = (self.hours as f64)
            + (self.minutes as f64 / 60.0)
            + (self.seconds / 3600.0);
        if self.sign { -abs_val } else { abs_val }
    }

    /// Converts to a `Degrees` type (1h = 15°).
    pub fn to_degrees(&self) -> Degrees {
        Degrees::new(self.to_decimal() * 15.0)
    }

    /// Converts to a `Radians` type.
    pub fn to_radians(&self) -> Radians {
        self.to_degrees().to::<Radian>()
    }
}

/// Implement `Display` for `HMS`.
impl std::fmt::Display for HMS {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = if self.sign { "-" } else { "" };
        write!(f, "{}{}h {}m {:.2}s", prefix, self.hours, self.minutes, self.seconds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hms_creation_normal() {
        let time = HMS::new(12, 34, 56.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 12);
        assert_eq!(time.minutes, 34);
        assert_eq!(time.seconds, 56.0);
    }

    #[test]
    fn test_hms_creation_seconds_overflow() {
        let time = HMS::new(12, 34, 120.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 12);
        assert_eq!(time.minutes, 36);
        assert_eq!(time.seconds, 0.0);
    }

    #[test]
    fn test_hms_creation_minutes_overflow() {
        let time = HMS::new(12, 120, 30.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 14);
        assert_eq!(time.minutes, 0);
        assert_eq!(time.seconds, 30.0);
    }

    #[test]
    fn test_hms_creation_negative_seconds() {
        let time = HMS::new(12, 34, -30.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 12);
        assert_eq!(time.minutes, 33);
        assert_eq!(time.seconds, 30.0);
    }

    #[test]
    fn test_hms_creation_negative_minutes() {
        let time = HMS::new(12, -30, 30.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 11);
        assert_eq!(time.minutes, 30);
        assert_eq!(time.seconds, 30.0);
    }

    #[test]
    fn test_hms_creation_negative_hours() {
        let time = HMS::new(-12, 34, 56.0);
        assert!(time.sign);
        assert_eq!(time.hours, 12);
        assert_eq!(time.minutes, 34);
        assert_eq!(time.seconds, 56.0);
    }

    #[test]
    fn test_hms_creation_all_overflow() {
        let time = HMS::new(12, 120, 120.0);
        assert!(!time.sign);
        assert_eq!(time.hours, 14);
        assert_eq!(time.minutes, 2);
        assert_eq!(time.seconds, 0.0);
    }

    #[test]
    fn test_to_decimal_positive() {
        let time = HMS::new(12, 34, 56.0);
        let decimal = time.to_decimal();
        assert!((decimal - 12.582222).abs() < 1e-6);
    }

    #[test]
    fn test_to_decimal_negative() {
        let time = HMS::new(-12, 34, 56.0);
        let decimal = time.to_decimal();
        assert!((decimal + 12.582222).abs() < 1e-6);
    }

    #[test]
    fn test_to_string_positive() {
        let time = HMS::new(12, 34, 56.0);
        assert_eq!(time.to_string(), "12h 34m 56.00s");
    }

    #[test]
    fn test_to_string_negative() {
        let time = HMS::new(-12, 34, 56.0);
        assert_eq!(time.to_string(), "-12h 34m 56.00s");
    }
}
