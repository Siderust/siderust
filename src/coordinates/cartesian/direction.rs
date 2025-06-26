use crate::coordinates::{
    frames, centers,
    kinds::DirectionKind,
};
use crate::units::Unit;
use super::Vector;

pub type Direction<C, F> = Vector<C, F, f64, DirectionKind>;

impl<C, F> Direction<C, F>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    /// Returns a Position vector in the same direction, scaled by the given magnitude.
    pub fn position<U: Unit>(&self, magnitude: U) -> super::Position<C, F, U> {
        super::Position::new(
            magnitude * self.x(),
            magnitude * self.y(),
            magnitude * self.z()
        )
    }
}
