use super::SphericalCoord;
use crate::coordinates::kinds::PositionKind;

pub type Position<C, F, U=f64>  = SphericalCoord<C, F, U, PositionKind>;

use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Unit;

impl<C, F, U> Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    /// The zero point (origin) in this coordinate system.
    pub const CENTER: Self = Self::from_degrees(0.0, 0.0, Some(0.0));

    pub fn direction(&self) -> super::Direction<C, F, U> {
        super::Direction::new_spherical_coord(self.polar, self.azimuth, None)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    
    #[test]
    fn test_spherical_coord_creation() {
        let coord = Position::<Barycentric, ICRS>::new(Degrees::new(45.0), Degrees::new(90.0), 1.0);
        assert_eq!(coord.ra().as_f64(), 45.0);
        assert_eq!(coord.dec().as_f64(), 90.0);
        assert_eq!(coord.distance.unwrap(), 1.0);
    }
    
    #[test]
    fn test_spherical_coord_to_string() {
        let coord = Position::<Geocentric, ICRS>::new(Degrees::new(30.0), Degrees::new(60.0), 1000.0);
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 60"));
        assert!(coord_string.contains("φ: 30"));
        assert!(coord_string.contains("r: 1000"));
    }
    
    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = Position::<Heliocentric, ICRS>::new(Degrees::new(0.0), Degrees::new(0.0), 0.0);
        assert_eq!(coord.polar.as_f64(), 0.0);
        assert_eq!(coord.azimuth.as_f64(), 0.0);
        assert_eq!(coord.distance.unwrap(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = Position::<Barycentric, ICRS>::new(Degrees::new(90.654321), Degrees::new(45.123456), 1234.56789);
        assert!((coord.dec().as_f64() - 45.123456).abs() < 1e-6);
        assert!((coord.ra().as_f64() - 90.654321).abs() < 1e-6);
        assert!((coord.distance.unwrap() - 1234.56789).abs() < 1e-6);
    }
}
