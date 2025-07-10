//! Astronomical Unit (AU) type and conversions.
//!
//! Provides a strongly-typed representation of a length in Astronomical Units (AU)
//! and conversions to and from Light Years (LY).

impl super::AstronomicalUnits {

    pub const KM_PER_AU: f64 = 149_597_870.7;

    pub const fn to_light_year(&self) -> super::LightYears {
        super::LightYears::new(self.0 / super::LightYears::AU_PER_LY)
    }

}

impl std::fmt::Display for super::AstronomicalUnits {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} Au", self.value())
    }
}


/// Converts an `AstronomicalUnit` to a `LightYear`.
impl From<super::Quantity<super::AstronomicalUnit>> for super::Quantity<super::LightYear> {
    fn from(au: super::Quantity<super::AstronomicalUnit>) -> Self {
        au.to_light_year()
    }
}
