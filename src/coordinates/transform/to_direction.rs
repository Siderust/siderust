use crate::coordinates::{cartesian, centers::ReferenceCenter, frames::ReferenceFrame, spherical};
use crate::units::LengthUnit;

impl<C, F, U> From<&spherical::Position<C, F, U>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    fn from(position: &spherical::Position<C, F, U>) -> Self {
        position.direction()
    }
}

impl<C, F, U> From<&cartesian::Position<C, F, U>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    fn from(position: &cartesian::Position<C, F, U>) -> Self {
        position.direction()
    }
}
