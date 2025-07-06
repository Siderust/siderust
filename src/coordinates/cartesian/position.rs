use crate::coordinates::{
    centers, frames,
    kinds::{PositionKind, VelocityKind}
};
use super::Vector;
use crate::units::{Quantity, LengthUnit};

pub type Position<C, F, U>  = Vector<C, F, U, PositionKind>;
pub type Velocity<C, F, U>  = Vector<C, F, U, VelocityKind>;

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{
    pub const CENTER: Self = Self::new_const(Quantity::<U>::new(0.0), Quantity::<U>::new(0.0), Quantity::<U>::new(0.0));

    pub fn direction(&self) -> super::Direction<C, F> {
        let d = self.distance();
        super::Direction::<C, F>::new(
            Quantity::<f64>::new(self.x() / d),
            Quantity::<f64>::new(self.y() / d),
            Quantity::<f64>::new(self.z() / d)
        )
    }
}

pub type ICRS<U>       = Position<centers::Barycentric,  frames::ICRS,       U>;
pub type HCRS<U>       = Position<centers::Heliocentric, frames::ICRS,       U>;
pub type GCRS<U>       = Position<centers::Geocentric,   frames::ICRS,       U>;
pub type TCRS<U>       = Position<centers::Topocentric,  frames::ICRS,       U>;
pub type Ecliptic<U>   = Position<centers::Heliocentric, frames::Ecliptic,   U>;
pub type Equatorial<U> = Position<centers::Geocentric,   frames::Equatorial, U>;
pub type Horizontal<U> = Position<centers::Topocentric,  frames::Horizontal, U>;
pub type Geographic<U> = Position<centers::Geocentric,   frames::ECEF,       U>;
