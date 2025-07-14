use crate::coordinates::{
    frames, centers,
    kinds::DirectionKind,
};
use crate::units::{Quantity, LengthUnit};
use super::Vector;

// TODO: Remove Unit from Direction
pub type Direction<C, F> = Vector<C, F, f64, DirectionKind>;

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
}

pub type Ecliptic<C=centers::Heliocentric>  = Direction<C, frames::Ecliptic>;
pub type Equatorial<C=centers::Geocentric>  = Direction<C, frames::Equatorial>;
pub type Horizontal<C=centers::Topocentric> = Direction<C, frames::Horizontal>;
pub type Geographic<C=centers::Geocentric>  = Direction<C, frames::ECEF>;
pub type ICRS<C=centers::Barycentric>       = Direction<C, frames::ICRS>;
pub type HCRS       = Direction<centers::Heliocentric, frames::ICRS>;
pub type GCRS       = Direction<centers::Geocentric,   frames::ICRS>;
pub type TCRS       = Direction<centers::Topocentric,  frames::ICRS>;
