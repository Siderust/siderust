use crate::coordinates::{frames, centers};
use crate::units::{Quantity, LengthUnit, Unitless};
use super::Vector;

pub type Direction<C, F> = Vector<C, F, Unitless>;

impl<C, F> Direction<C, F>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    /// Returns a Position vector in the same direction, scaled by the given magnitude.
    pub fn position<U: LengthUnit>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U> {
        super::Position::new(
            magnitude * self.x().value(),
            magnitude * self.y().value(),
            magnitude * self.z().value()
        )
    }

    pub fn normalize(x: f64, y: f64, z: f64) -> Self {
        let norm = nalgebra::Vector3::<f64>::new(x, y, z).normalize();
        Self::new(
            Quantity::<Unitless>::new(norm.x),
            Quantity::<Unitless>::new(norm.y),
            Quantity::<Unitless>::new(norm.z)
        )
    }
}

pub type Ecliptic<C=centers::Heliocentric>  = Direction<C, frames::Ecliptic>;
pub type Equatorial<C=centers::Geocentric>  = Direction<C, frames::Equatorial>;
pub type Horizontal<C=centers::Topocentric> = Direction<C, frames::Horizontal>;
pub type Geographic<C=centers::Geocentric>  = Direction<C, frames::ECEF>;
pub type ICRS<C=centers::Barycentric>       = Direction<C, frames::ICRS>;
pub type HCRS = Direction<centers::Heliocentric, frames::ICRS>;
pub type GCRS = Direction<centers::Geocentric,   frames::ICRS>;
pub type TCRS = Direction<centers::Topocentric,  frames::ICRS>;
