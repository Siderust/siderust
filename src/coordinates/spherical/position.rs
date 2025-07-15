use super::SphericalCoord;

use crate::coordinates::{
    centers, frames,
    kinds::PositionKind
};

use crate::units::Kilometer;
pub type Position<C, F, U>  = SphericalCoord<C, F, U, PositionKind>;
pub type Ecliptic<U>    = Position<centers::Heliocentric, frames::Ecliptic, U>;   // L (l), B (b), R (r)
pub type Equatorial<U>  = Position<centers::Geocentric,   frames::Equatorial, U>; // Dec (δ), RA (α), LengthUnit (d)
pub type Horizontal<U>  = Position<centers::Topocentric,  frames::Horizontal, U>; // Alt (α), Az (θ), LengthUnit (d)
pub type ICRS<U>        = Position<centers::Barycentric,  frames::ICRS, U>; // Dec (δ), RA (α), LengthUnit (d)
pub type HCRS<U>        = Position<centers::Heliocentric, frames::ICRS, U>; // Dec (δ), RA (α), LengthUnit (d)
pub type GCRS<U>        = Position<centers::Geocentric,   frames::ICRS, U>; // Dec (δ), RA (α), LengthUnit (d)
pub type Geographic     = Position<centers::Geocentric,   frames::ECEF, Kilometer>; //Latitude (φ),Longitude (λ), Altitude (h)

use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::LengthUnit;

impl<C, F, U> Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    /// The zero point (origin) in this coordinate system.
    pub const CENTER: Self = Self::from_degrees(0.0, 0.0, None);

    pub fn direction(&self) -> super::Direction<C, F> {
        super::Direction::new_spherical_coord(self.polar, self.azimuth, None)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{AstronomicalUnit, Degrees, AU};
    
    #[test]
    fn test_spherical_coord_creation() {
        let coord = ICRS::<AstronomicalUnit>::new(Degrees::new(45.0), Degrees::new(90.0), 1.0);
        assert_eq!(coord.ra().value(), 45.0);
        assert_eq!(coord.dec().value(), 90.0);
        assert_eq!(coord.distance.unwrap(), 1.0);
    }
    
    #[test]
    fn test_spherical_coord_to_string() {
        let coord = GCRS::<AstronomicalUnit>::new(Degrees::new(30.0), Degrees::new(60.0), 1000.0);
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 60"));
        assert!(coord_string.contains("φ: 30"));
        assert!(coord_string.contains("r: 1000"));
    }
    
    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = HCRS::<AstronomicalUnit>::new(Degrees::new(0.0), Degrees::new(0.0), 0.0);
        assert_eq!(coord.polar.value(), 0.0);
        assert_eq!(coord.azimuth.value(), 0.0);
        assert_eq!(coord.distance.unwrap(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = ICRS::<AstronomicalUnit>::new(Degrees::new(90.654321), Degrees::new(45.123456), 1234.56789);
        assert!((coord.dec().value() - 45.123456).abs() < 1e-6);
        assert!((coord.ra().value() - 90.654321).abs() < 1e-6);
        assert!((coord.distance.unwrap() - 1234.56789*AU).abs() < 1e-6*AU);
    }
}
