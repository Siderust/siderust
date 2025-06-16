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
        spherical::Direction::new_spherical_coord(position.polar, position.azimuth, None)
    }
}

impl<C, F> From<&cartesian::Position<C, F,>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(position: &cartesian::Position<C, F>) -> Self {
        cartesian::Direction::from_vec3(position.as_vec3().normalize())
    }
}
