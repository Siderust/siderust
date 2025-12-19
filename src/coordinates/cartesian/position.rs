use crate::coordinates::{centers, frames};
use qtty::Simplify;
use qtty::{LengthUnit, Quantity};

// TODO: Bound U to LengthUnit
// see issue #112792 <https://github.com/rust-lang/rust/issues/112792> for more information
//pub type Position<C, F, U: LengthUnit> = super::Vector<C, F, U>;
pub type Position<C, F, U> = super::Vector<C, F, U>;

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{
    /// Returns the direction (unit vector) from the center to this position.
    ///
    /// Note: Directions are frame-only types (no center). This extracts the
    /// normalized direction regardless of the position's center.
    pub fn direction(&self) -> super::Direction<F> {
        use crate::coordinates::spherical::direction::DirectionUnit;
        let d = self.distance();
        super::Direction::<F>::from_vec3(nalgebra::Vector3::new(
            Quantity::<DirectionUnit>::new((self.x() / d).simplify().value()),
            Quantity::<DirectionUnit>::new((self.y() / d).simplify().value()),
            Quantity::<DirectionUnit>::new((self.z() / d).simplify().value()),
        ))
    }
}

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{
    /// The origin of this coordinate system (all coordinates 0). AKA Null Vector.
    pub const CENTER: Self = Self::new_const(
        (),
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
    );
}

pub type Ecliptic<U, C = centers::Heliocentric> = Position<C, frames::Ecliptic, U>;
pub type Equatorial<U, C = centers::Geocentric> = Position<C, frames::Equatorial, U>;
pub type Horizontal<U, C = centers::Topocentric> = Position<C, frames::Horizontal, U>;
pub type Geographic<U, C = centers::Geocentric> = Position<C, frames::ECEF, U>;
pub type ICRS<U, C = centers::Barycentric> = Position<C, frames::ICRS, U>;
pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;
pub type TCRS<U> = Position<centers::Topocentric, frames::ICRS, U>;
