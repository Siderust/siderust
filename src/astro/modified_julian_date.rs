use chrono::{DateTime, Utc};

/// Represents Modified Julian Date (MJD), which is the Julian Date
/// minus 2400000.5, used in various scientific and technical applications.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModifiedJulianDate(f64);

impl ModifiedJulianDate {
    pub fn new(mjd: f64) -> Self {
        ModifiedJulianDate(mjd)
    }

    /// Returns the inner Modified Julian Day value.
    pub const fn value(&self) -> f64 {
        self.0
    }

    pub const fn to_julian_day(&self) -> super::JulianDate {
        super::JulianDate::new(self.value() + 2400000.5)
    }

    pub fn to_utc(&self) -> Option<DateTime<Utc>> {
        let jd = self.to_julian_day().value();
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = (jd - unix_epoch_jd) * 86400.0;
        let secs = seconds_since_epoch.floor() as i64;
        let nanos = ((seconds_since_epoch - secs as f64) * 1e9) as u32;
        DateTime::<Utc>::from_timestamp(secs, nanos)
    }

    pub fn from_utc(datetime: DateTime<Utc>) -> Self {
        let unix_epoch_jd = 2440587.5; // Julian Day for Unix epoch (1970-01-01T00:00:00Z)
        let seconds_since_epoch = datetime.timestamp() as f64;
        let nanos = datetime.timestamp_subsec_nanos() as f64 / 1e9;
        let jd = unix_epoch_jd + (seconds_since_epoch + nanos) / 86400.0;
        ModifiedJulianDate::new(jd - 2400000.5)
    }
}

impl From<ModifiedJulianDate> for super::JulianDate {
    fn from(mjd: ModifiedJulianDate) -> Self {
        mjd.to_julian_day()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::DateTime;

    #[test]
    fn test_modified_julian_day_creation() {
        let mjd = ModifiedJulianDate::new(51544.5);
        assert_eq!(mjd.value(), 51544.5);
    }

    #[test]
    fn test_to_julian_day() {
        let mjd = ModifiedJulianDate::new(51544.5);
        assert_eq!(mjd.to_julian_day().value(), 2451545.0);
    }

    #[test]
    fn test_to_utc() {
        let mjd = ModifiedJulianDate::new(51544.5);
        let datetime = mjd.to_utc();
        assert_eq!(datetime, DateTime::from_timestamp(946728000, 0));
    }

    #[test]
    fn test_from_utc() {
        let datetime = DateTime::from_timestamp(946728000, 0);
        let mjd = ModifiedJulianDate::from_utc(datetime.unwrap());
        assert_eq!(mjd.value(), 51544.5);
    }
}
