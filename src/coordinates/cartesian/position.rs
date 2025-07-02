use crate::coordinates::{
    cartesian::Direction, centers, frames, kinds::{PositionKind, VelocityKind}
};
use super::Vector;
use crate::units::Unit;

pub type Position<C, F, U = f64>  = Vector<C, F, U, PositionKind>;
pub type Velocity<C, F, U = f64>  = Vector<C, F, U, VelocityKind>;

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
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

// === ICRS-based Cartesian coordinate types ===
pub type ICRS<U = f64> = Position<centers::Barycentric,  frames::ICRS, U>;
pub type HCRS<U = f64> = Position<centers::Heliocentric, frames::ICRS, U>;
pub type GCRS<U = f64> = Position<centers::Geocentric,   frames::ICRS, U>;
pub type TCRS<U = f64> = Position<centers::Topocentric,  frames::ICRS, U>;

// === Ecliptic frame ===
pub type EclipticPos<U = f64> = Position<centers::Heliocentric, frames::Ecliptic, U>;

// === Equatorial frame ===
pub type EquatorialPos<U = f64>   = Position<centers::Geocentric,   frames::Equatorial, U>;

// === Horizontal and Earth-fixed frame ===
pub type HorizontalPos<U = f64>  = Position<centers::Topocentric,  frames::Horizontal, U>;

// === Geographic and Earth-fixed frames ===
pub type GeographicPos<U = f64> = Position<centers::Geocentric, frames::ECEF, U>;
