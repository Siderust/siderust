use crate::coordinates::{
    cartesian::Direction, centers, frames, kinds::{PositionKind, VelocityKind}
};
use super::Vector;
use crate::units::Distance;

pub type Position<C, F, U = f64>  = Vector<C, F, U, PositionKind>;
pub type Velocity<C, F, U = f64>  = Vector<C, F, U, VelocityKind>;

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Distance,
{
    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in units of U.
    pub fn distance(&self) -> U {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        // TODO: Use nalgebra for this?
        let norm = nalgebra::Vector3::<f64>::new(
            self.x().into(),
            self.y().into(),
            self.z().into()
        ).normalize();
        Self::new(
            U::from(norm.x),
            U::from(norm.y),
            U::from(norm.z)
        )
    }

    pub fn direction(&self) -> super::Direction<C, F> {
        // TODO: Use nalgebra for this?
        let norm = nalgebra::Vector3::<f64>::new(
            self.x().into(),
            self.y().into(),
            self.z().into()
        ).normalize();
        Direction::new(
            norm.x,
            norm.y,
            norm.z
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
