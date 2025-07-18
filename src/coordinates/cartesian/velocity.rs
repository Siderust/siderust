use crate::coordinates::frames;
use super::Vector;

// TODO: Bound U to VelocityUnit
// see issue #112792 <https://github.com/rust-lang/rust/issues/112792> for more information
//pub type Velocity<C, F, U: VelocityUnit> = Vector<C, F, U>;
pub type Velocity<F, U> = Vector<(), F, U>;

pub type Ecliptic<U>   = Velocity<frames::Ecliptic,   U>;
pub type Equatorial<U> = Velocity<frames::Equatorial, U>;
pub type Horizontal<U> = Velocity<frames::Horizontal, U>;
pub type ICRS<U>       = Velocity<frames::ICRS,       U>;
