use crate::coordinates::{
    frames, centers,
    kinds::DirectionKind,
};
use crate::units::Unit;
use super::Vector;

pub type Direction<C, F, U = f64> = Vector<C, F, U, DirectionKind>;

impl<C, F, U> Direction<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// For a unit direction vector, the "distance" is always 1 in its unit system.
    pub fn distance(&self) -> U {
        U::from(1.0)
    }

    /// Returns a Position vector in the same direction, scaled by the given magnitude.
    pub fn position(&self, magnitude: U) -> super::Position<C, F, U> {
        super::Position::from_vec3(self.as_vec3() * magnitude)
    }
}
