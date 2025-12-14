use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::{centers, frames};
use qtty::{LengthUnit, Quantity};

pub type Direction<C, F> = super::Vector<C, F, DirectionUnit>;

impl<C, F> Direction<C, F>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    /// Returns a Position vector in the same direction, scaled by the given magnitude.
    pub fn position<U: LengthUnit>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U> {
        super::Position::new(
            magnitude * self.x().value(),
            magnitude * self.y().value(),
            magnitude * self.z().value(),
        )
    }

    pub fn normalize(x: f64, y: f64, z: f64) -> Self {
        let norm = nalgebra::Vector3::<f64>::new(x, y, z).normalize();
        Self::new(
            Quantity::<DirectionUnit>::new(norm.x),
            Quantity::<DirectionUnit>::new(norm.y),
            Quantity::<DirectionUnit>::new(norm.z),
        )
    }

    /// Returns a formatted string representation of this direction vector.
    /// This is provided as a method because DirectionUnit doesn't implement Display.
    pub fn display(&self) -> String {
        format!(
            "Center: {}, Frame: {}, X: {:.6}, Y: {:.6}, Z: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.x().value(),
            self.y().value(),
            self.z().value()
        )
    }
}

pub type Ecliptic<C = centers::Heliocentric> = Direction<C, frames::Ecliptic>;
pub type Equatorial<C = centers::Geocentric> = Direction<C, frames::Equatorial>;
pub type Horizontal<C = centers::Topocentric> = Direction<C, frames::Horizontal>;
pub type Geographic<C = centers::Geocentric> = Direction<C, frames::ECEF>;
pub type ICRS<C = centers::Barycentric> = Direction<C, frames::ICRS>;
pub type HCRS = Direction<centers::Heliocentric, frames::ICRS>;
pub type GCRS = Direction<centers::Geocentric, frames::ICRS>;
pub type TCRS = Direction<centers::Topocentric, frames::ICRS>;

// Display is implemented in vector.rs for all Vector<C, F, U> where U: Unit
