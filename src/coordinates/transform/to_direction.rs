use crate::coordinates::{
    cartesian, spherical,
    frames::ReferenceFrame,
    centers::ReferenceCenter,
};


impl<C, F> From<&spherical::Position<C, F,>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(position: &spherical::Position<C, F>) -> Self {
        position.direction()
    }
}

impl<C, F> From<&cartesian::Position<C, F,>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(position: &cartesian::Position<C, F>) -> Self {
        position.direction()
    }
}
