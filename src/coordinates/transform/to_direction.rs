use crate::coordinates::{
    cartesian, spherical,
    frames::ReferenceFrame,
    centers::ReferenceCenter,
};
use crate::units::Unit;

impl<C, F, U> From<&spherical::Position<C, F, U>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    fn from(position: &spherical::Position<C, F, U>) -> Self {
        position.direction()
    }
}

impl<C, F, U> From<&cartesian::Position<C, F, U>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    fn from(position: &cartesian::Position<C, F, U>) -> Self {
        position.direction()
    }
}
