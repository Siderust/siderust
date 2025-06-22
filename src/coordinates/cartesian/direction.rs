use crate::coordinates::{
    frames, centers,
    kinds::DirectionKind,
};

use super::CartesianCoord;

pub type Direction<C, F> = CartesianCoord<C, F, DirectionKind>;

impl<C, F> Direction<C, F>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in AU.
    pub const fn distance(&self) -> f64 { 1.0 }

    pub fn position(&self, magnitude: f64) -> super::Position<C, F> {
        super::Position::from_vec3(self.as_vec3() * magnitude)
    }

}