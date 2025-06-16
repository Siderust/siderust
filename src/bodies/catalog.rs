//! # Star Catalog Module
//!
//! This module provides a set of constant definitions for well-known bright stars,
//! including their physical properties and equatorial coordinates at epoch J2000.0.
//! These constants can be used for quick access to reference stars in astronomical calculations.

use super::Star;
use crate::coordinates::{centers::Geocentric, frames::Equatorial, spherical::Position};
use crate::targets::Target;
use crate::units::*;

pub const VEGA: Star<'static> = Star::new_const(
    "Vega",
    LightYear::new(25.0),
    SolarMass::new(2.135),
    SolarRadius::new(2.59),
    SolarLuminosity::new(40.12),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(279.2347),
            Degrees::new(38.7837),
            25.0, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const POLARIS: Star<'static> = Star::new_const(
    "Polaris",
    LightYear::new(433.0),
    SolarMass::new(6.5),
    SolarRadius::new(46.0),
    SolarLuminosity::new(2500.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(37.95456067),
            Degrees::new(89.26410897),
            433.0, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const SIRIUS: Star<'static> = Star::new_const(
    "Sirius",
    LightYear::new(8.6),
    SolarMass::new(2.063),
    SolarRadius::new(1.713),
    SolarLuminosity::new(24.7),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(101.28715533),
            Degrees::new(-16.716115867),
            8.6, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const CANOPUS: Star<'static> = Star::new_const(
    "Canopus",
    LightYear::new(310.0),
    SolarMass::new(8.0),
    SolarRadius::new(71.0),
    SolarLuminosity::new(13_600.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(95.98787778),
            Degrees::new(-52.69566111),
            310.0, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const ARCTURUS: Star<'static> = Star::new_const(
    "Arcturus",
    LightYear::new(36.7),
    SolarMass::new(1.1),
    SolarRadius::new(26.0),
    SolarLuminosity::new(170.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(213.9153),
            Degrees::new(19.1825),
            36.7, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const RIGEL: Star<'static> = Star::new_const(
    "Rigel",
    LightYear::new(860.0),
    SolarMass::new(17.0),
    SolarRadius::new(78.9),
    SolarLuminosity::new(120_000.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(78.634467),
            Degrees::new(-8.20163889),
            860.0, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const BETELGEUSE: Star<'static> = Star::new_const(
    "Betelgeuse",
    LightYear::new(548.0),
    SolarMass::new(11.6),
    SolarRadius::new(724.0),
    SolarLuminosity::new(14_000.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(88.792939),
            Degrees::new(7.407064),
            548.0, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const PROCYON: Star<'static> = Star::new_const(
    "Procyon",
    LightYear::new(11.5),
    SolarMass::new(1.499),
    SolarRadius::new(2.048),
    SolarLuminosity::new(6.93),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(114.825493),
            Degrees::new(5.224993),
            11.5, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const ALDEBARAN: Star<'static> = Star::new_const(
    "Aldebaran",
    LightYear::new(65.1),
    SolarMass::new(1.16),
    SolarRadius::new(45.1),
    SolarLuminosity::new(439.0),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(68.980163),
            Degrees::new(16.509302),
            65.1, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);

pub const ALTAIR: Star<'static> = Star::new_const(
    "Altair",
    LightYear::new(16.7),
    SolarMass::new(1.86),
    SolarRadius::new(1.79),
    SolarLuminosity::new(10.6),
    Target::<Position::<Geocentric, Equatorial>>::new_static(
        Position::<Geocentric, Equatorial>::new_const(
            Degrees::new(297.695827),
            Degrees::new(8.868321),
            16.7, // Distance in Light Years
        ),
        JulianDay::J2000,
    ),
);
