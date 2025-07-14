use crate::astro::JulianDate;

/// Represents a culmination event — the moment a celestial body crosses
/// the observer’s meridian.
///
/// - `Upper`: transit across the upper meridian (highest altitude).  
/// - `Lower`: transit across the lower meridian (lowest altitude).
///
/// The `jd` field stores the Julian Day of the event.
#[derive(Debug)]
pub enum Culmination {
    Upper { jd: JulianDate },
    Lower { jd: JulianDate },
}
