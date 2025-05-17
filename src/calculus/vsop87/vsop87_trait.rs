use crate::coordinates::{
    CartesianCoord,
    centers::{Heliocentric, Barycentric},
    frames::Ecliptic,
};
use crate::units::JulianDay;
use crate::targets::Target;
use crate::bodies::solar_system::*;

pub trait VSOP87 {
    fn vsop87a(&self, jd: JulianDay) -> Target<CartesianCoord<Heliocentric, Ecliptic>>;
    fn vsop87e(&self, jd: JulianDay) -> Target<CartesianCoord<Barycentric, Ecliptic>>;
}

macro_rules! impl_vsop87_for_planet {
    ($planet:ident) => {
        impl VSOP87 for $planet {
            fn vsop87a(&self, jd: JulianDay) -> Target<CartesianCoord<Heliocentric, Ecliptic>> {
                $planet::vsop87a(jd)
            }

            fn vsop87e(&self, jd: JulianDay) -> Target<CartesianCoord<Barycentric, Ecliptic>> {
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
