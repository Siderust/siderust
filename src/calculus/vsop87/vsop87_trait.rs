use crate::{coordinates::{
    cartesian::Position,
    centers::{Barycentric, Heliocentric},
    frames::Ecliptic,
}, units::AstronomicalUnit};
use crate::units::JulianDay;
use crate::targets::Target;
use crate::bodies::solar_system::*;

pub trait VSOP87 {
    fn vsop87a(&self, jd: JulianDay) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>>;
    fn vsop87e(&self, jd: JulianDay) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>;
}

macro_rules! impl_vsop87_for_planet {
    ($planet:ident) => {
        impl VSOP87 for $planet {
            fn vsop87a(&self, jd: JulianDay) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>> {
                $planet::vsop87a(jd)
            }

            fn vsop87e(&self, jd: JulianDay) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>> {
                $planet::vsop87e(jd)
            }
        }
    };
}

impl_vsop87_for_planet!(Mercury);
impl_vsop87_for_planet!(Venus);
impl_vsop87_for_planet!(Earth);
impl_vsop87_for_planet!(Mars);
impl_vsop87_for_planet!(Jupiter);
impl_vsop87_for_planet!(Saturn);
impl_vsop87_for_planet!(Uranus);
impl_vsop87_for_planet!(Neptune);
