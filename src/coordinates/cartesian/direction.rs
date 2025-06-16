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
    pub const fn distance_from_origin(&self) -> f64 { 1.0 }

}