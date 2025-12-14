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
    pub const CENTER: Self = Self::new_const(
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
    );

    pub fn direction(&self) -> super::Direction<C, F> {
        use crate::coordinates::spherical::direction::DirectionUnit;
        let d = self.distance();
        super::Direction::<C, F>::new_const(
            Quantity::<DirectionUnit>::new((self.x() / d).simplify().value()),
            Quantity::<DirectionUnit>::new((self.y() / d).simplify().value()),
            Quantity::<DirectionUnit>::new((self.z() / d).simplify().value()),
        )
    }
}

pub type Ecliptic<U, C = centers::Heliocentric> = Position<C, frames::Ecliptic, U>;
pub type Equatorial<U, C = centers::Geocentric> = Position<C, frames::Equatorial, U>;
pub type Horizontal<U, C = centers::Topocentric> = Position<C, frames::Horizontal, U>;
pub type Geographic<U, C = centers::Geocentric> = Position<C, frames::ECEF, U>;
pub type ICRS<U, C = centers::Barycentric> = Position<C, frames::ICRS, U>;
pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;
pub type TCRS<U> = Position<centers::Topocentric, frames::ICRS, U>;
