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

pub type ICRS       = Direction<centers::Barycentric,  frames::ICRS>;
pub type HCRS       = Direction<centers::Heliocentric, frames::ICRS>;
pub type GCRS       = Direction<centers::Geocentric,   frames::ICRS>;
pub type TCRS       = Direction<centers::Topocentric,  frames::ICRS>;
pub type Ecliptic   = Direction<centers::Heliocentric, frames::Ecliptic>;
pub type Equatorial = Direction<centers::Geocentric,   frames::Equatorial>;
pub type Horizontal = Direction<centers::Topocentric,  frames::Horizontal>;
pub type Geographic = Direction<centers::Geocentric,   frames::ECEF>;
