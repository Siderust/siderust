// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar System
//!
//! This module exposes a type‑safe, canonical description of our Solar System
//! for astronomical calculations, ephemeris generation, and visualisation
//! engines.  It bundles **all** major bodies and selected Lagrange points into
//! a single aggregate constant: [`SOLAR_SYSTEM`].
//!
//! ## What’s Included
//!
//! | Category                | Bodies / Points                                                                                 |
//! |-------------------------|--------------------------------------------------------------------------------------------------|
//! | **Star**                | Sun                                                                                             |
//! | **Planets**             | Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune                                   |
//! | **Dwarf planets**       | Ceres, Pluto, Haumea, Makemake, Eris                                                            |
//! | **Major moons**         | Moon, Io, Europa, Ganymede, Callisto, Titan, Triton                                             |
//! | **Lagrange points**     | Sun–Earth L₁, Sun–Earth L₂                                                                      |
//!
//! Each body re‑uses the strongly‑typed structures defined elsewhere in the crate
//! (`Star`, `Planet`, `Satellite`, [`Orbit`], …) so values can be consumed
//! directly without conversion.
//!
//! ### Physical & Orbital Parameters
//!
//! Every object comes with the following fields (J2000.0 epoch unless noted):
//!
//! * **Mass** (kg)
//! * **Mean radius** (km)
//! * **Keplerian orbital elements**:
//!   * `a` – semi‑major axis [`AstronomicalUnit`]
//!   * `e` – eccentricity (unitless)
//!   * `i` – inclination [`Degrees`]
//!   * `Ω` – longitude of ascending node [`Degrees`]
//!   * `ω` – argument of perihelion [`Degrees`]
//!   * `M₀` – mean anomaly at epoch [`Degrees`]
//!
//! ```text
//! Abbreviations for orbital elements:
//!   a  — semi‑major axis         Ω — longitude of ascending node
//!   e  — eccentricity            ω — argument of perihelion
//!   i  — inclination             M₀ — mean anomaly at epoch
//! ```
//!
//! ---
//! ## References
//!
//! 1. NASA – Planetary Fact Sheet: <https://nssdc.gsfc.nasa.gov/planetary/factsheet/>
//! 2. Williams, D. R. (2024). *Planetary Fact Sheet – Metric*. NASA Goddard Space Flight Center.

use super::{Planet, Satellite, Star};
use crate::astro::{orbit::Orbit, JulianDate};
use crate::coordinates::spherical::position::{Ecliptic, EquatorialMeanJ2000};
use crate::targets::Target;
use qtty::length::nominal::RSUN;
use qtty::*;

pub struct Sun;
pub struct Mercury;
pub struct Venus;
pub struct Earth;
pub struct Moon;
pub struct Mars;
pub struct Jupiter;
pub struct Saturn;
pub struct Uranus;
pub struct Neptune;

/// **Sun** – the central star of the Solar System.
///
/// | Parameter       | Value                |
/// |-----------------|----------------------|
/// | Mass            | 1.9885 × 10³⁰ kg     |
/// | Radius          | 696,340 km           |
/// | Luminosity      | 3.828 × 10²⁶ W       |
/// | Right Ascension | 18h 44m 48s (J2000)  |
/// | Declination     | −23° 00′ 00″ (J2000) |
/// | LengthUnit      | 1 AU (~0.0000158 ly) |
pub const SUN: super::Star<'static> = super::Star::new_const(
    "Sun",
    LightYears::new(1.58125e-5), // 1 AstronomicalUnits in LightYears
    SolarMasses::new(1.0),
    RSUN,
    L_SUN,
    Target::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            HourAngles::from_hms(-23, 0, 0.0).to_const::<Degree>(), // Dec (polar) Approx at J2000
            HourAngles::from_hms(18, 44, 48.0).to_const::<Degree>(), // RA (azimuth) Approx at J2000
            LightYears::new(1.58125e-5), // 1 AstronomicalUnits in LightYears
        ),
        JulianDate::J2000,
    ),
);

/// **Mercury** – innermost planet of the Solar System.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 3.3011 × 10²³ kg |
/// | Radius    | 2,439.7 km       |
/// | a         | 0.387099 AU      |
/// | e         | 0.20563          |
/// | i         | 7.00487°         |
/// | Ω         | 48.33167°        |
/// | ω         | 29.12478°        |
/// | M₀        | 174.79439°       |
pub const MERCURY: super::Planet = super::Planet {
    mass: Kilograms::new(3.3011e23),
    radius: Kilometers::new(2439.7),
    orbit: Orbit::new(
        AstronomicalUnits::new(0.38709893),
        0.20563069,
        Degrees::new(7.00487),
        Degrees::new(48.33167),
        Degrees::new(29.12478),
        Degrees::new(174.79439),
        JulianDate::J2000,
    ),
};

/// **Venus** – second planet from the Sun.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 4.8675 × 10²⁴ kg |
/// | Radius    | 6,051.8 km       |
/// | a         | 0.723332 AU      |
/// | e         | 0.006773         |
/// | i         | 3.39471°         |
/// | Ω         | 76.68069°        |
/// | ω         | 54.85229°        |
/// | M₀        | 50.44675°        |
pub const VENUS: super::Planet = super::Planet {
    mass: Kilograms::new(4.8675e24),
    radius: Kilometers::new(6051.8),
    orbit: Orbit::new(
        AstronomicalUnits::new(0.72333199),
        0.00677323,
        Degrees::new(3.39471),
        Degrees::new(76.68069),
        Degrees::new(54.85229),
        Degrees::new(50.44675),
        JulianDate::J2000,
    ),
};

/// **Earth** – third planet from the Sun and our home.
///
/// | Parameter | Value             |
/// |-----------|-------------------|
/// | Mass      | 5.97237 × 10²⁴ kg |
/// | Radius    | 6,371.0 km        |
/// | a         | 1.000000 AU       |
/// | e         | 0.016710          |
/// | i         | 0.00005°          |
/// | Ω         | -11.26064°        |
/// | ω         | 114.20783°        |
/// | M₀        | 357.51716°        |
pub const EARTH: super::Planet = super::Planet {
    mass: Kilograms::new(5.97237e24),
    radius: Kilometers::new(6371.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(1.00000011),
        0.01671022,
        Degrees::new(0.00005),
        Degrees::new(-11.26064),
        Degrees::new(114.20783),
        Degrees::new(357.51716),
        JulianDate::J2000,
    ),
};

/// **Moon** – natural satellite of Earth.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 7.346 × 10²² kg  |
/// | Radius    | 1,737.4 km       |
/// | a         | ~0.00257 AU      |
/// | e         | 0.0549           |
/// | i         | 5.145°           |
/// | Ω         | 125.08°          |
/// | ω         | 318.15°          |
/// | M₀        | 135.27°          |
pub const MOON: super::Satellite = super::Satellite::new_const(
    "Moon",
    Kilograms::new(7.346e22),
    Kilometers::new(1_737.4),
    Orbit {
        // 384 400 km → AstronomicalUnits
        semi_major_axis: AstronomicalUnits::new(2.566881e-6),
        eccentricity: 0.054_9,
        inclination: Degrees::new(5.145),
        longitude_of_ascending_node: Degrees::new(125.08),
        argument_of_perihelion: Degrees::new(318.15),
        mean_anomaly_at_epoch: Degrees::new(135.27),
        epoch: JulianDate::J2000,
    },
);

/// **Mars** – the red planet, fourth from the Sun.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 6.4171 × 10²³ kg |
/// | Radius    | 3,389.5 km       |
/// | a         | 1.523662 AU      |
/// | e         | 0.093412         |
/// | i         | 1.85061°         |
/// | Ω         | 49.57854°        |
/// | ω         | 286.46230°       |
/// | M₀        | 19.41248°        |
pub const MARS: super::Planet = super::Planet {
    mass: Kilograms::new(6.4171e23),
    radius: Kilometers::new(3389.5),
    orbit: Orbit::new(
        AstronomicalUnits::new(1.52366231),
        0.09341233,
        Degrees::new(1.85061),
        Degrees::new(49.57854),
        Degrees::new(286.46230),
        Degrees::new(19.41248),
        JulianDate::J2000,
    ),
};

/// **Jupiter** – the largest planet in the Solar System.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 1.8982 × 10²⁷ kg |
/// | Radius    | 69,911.0 km      |
/// | a         | 5.203363 AU      |
/// | e         | 0.048393         |
/// | i         | 1.30530°         |
/// | Ω         | 100.55615°       |
/// | ω         | 274.19770°       |
/// | M₀        | 19.65053°        |
pub const JUPITER: super::Planet = super::Planet {
    mass: Kilograms::new(1.8982e27),
    radius: Kilometers::new(69911.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(5.20336301),
        0.04839266,
        Degrees::new(1.30530),
        Degrees::new(100.55615),
        Degrees::new(274.19770),
        Degrees::new(19.65053),
        JulianDate::J2000,
    ),
};

/// **Saturn** – known for its prominent ring system.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 5.6834 × 10²⁶ kg |
/// | Radius    | 58,232.0 km      |
/// | a         | 9.537070 AU      |
/// | e         | 0.054151         |
/// | i         | 2.48446°         |
/// | Ω         | 113.71504°       |
/// | ω         | 338.71690°       |
/// | M₀        | 317.51238°       |
pub const SATURN: super::Planet = super::Planet {
    mass: Kilograms::new(5.6834e26),
    radius: Kilometers::new(58232.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(9.53707032),
        0.05415060,
        Degrees::new(2.48446),
        Degrees::new(113.71504),
        Degrees::new(338.71690),
        Degrees::new(317.51238),
        JulianDate::J2000,
    ),
};

/// **Uranus** – icy giant with extreme axial tilt.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 8.6810 × 10²⁵ kg |
/// | Radius    | 25,362.0 km      |
/// | a         | 19.191264 AU     |
/// | e         | 0.047168         |
/// | i         | 0.76986°         |
/// | Ω         | 74.22988°        |
/// | ω         | 96.73436°        |
/// | M₀        | 142.26794°       |
pub const URANUS: super::Planet = super::Planet {
    mass: Kilograms::new(8.6810e25),
    radius: Kilometers::new(25362.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(19.19126393),
        0.04716771,
        Degrees::new(0.76986),
        Degrees::new(74.22988),
        Degrees::new(96.73436),
        Degrees::new(142.26794),
        JulianDate::J2000,
    ),
};

/// **Neptune** – outermost giant planet.
///
/// | Parameter | Value             |
/// |-----------|-------------------|
/// | Mass      | 1.02409 × 10²⁶ kg |
/// | Radius    | 24,622.0 km       |
/// | a         | 30.068963 AU      |
/// | e         | 0.008586          |
/// | i         | 1.76917°          |
/// | Ω         | 131.72169°        |
/// | ω         | 273.24966°        |
/// | M₀        | 259.90868°        |
pub const NEPTUNE: super::Planet = super::Planet {
    mass: Kilograms::new(1.02409e26),
    radius: Kilometers::new(24622.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(30.06896348),
        0.00858587,
        Degrees::new(1.76917),
        Degrees::new(131.72169),
        Degrees::new(273.24966),
        Degrees::new(259.90868),
        JulianDate::J2000,
    ),
};

/// **Pluto** – dwarf planet once considered the ninth planet.
///
/// | Parameter | Value            |
/// |-----------|------------------|
/// | Mass      | 1.303 × 10²² kg  |
/// | Radius    | 1,188.3 km       |
/// | a         | 39.481687 AU     |
/// | e         | 0.248808         |
/// | i         | 17.14175°        |
/// | Ω         | 110.30347°       |
/// | ω         | 113.76329°       |
/// | M₀        | 14.86205°        |
pub const PLUTO: super::Planet = super::Planet {
    mass: Kilograms::new(1.303e22),
    radius: Kilometers::new(1188.3),
    orbit: Orbit::new(
        AstronomicalUnits::new(39.48168677),
        0.24880766,
        Degrees::new(17.14175),
        Degrees::new(110.30347),
        Degrees::new(113.76329),
        Degrees::new(14.86205),
        JulianDate::J2000,
    ),
};

// -------------------------------------------------------------------------------------------------
//  Dwarf planets
// -------------------------------------------------------------------------------------------------

pub const CERES: Planet = Planet {
    mass: Kilograms::new(9.393e20),
    radius: Kilometers::new(473.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(2.7675),
        0.0758,
        Degrees::new(10.5941),
        Degrees::new(80.3055),
        Degrees::new(73.5977),
        Degrees::new(95.9892),
        JulianDate::J2000,
    ),
};

pub const HAUMEA: Planet = Planet {
    mass: Kilograms::new(4.006e21),
    radius: Kilometers::new(620.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(43.218),
        0.195,
        Degrees::new(28.19),
        Degrees::new(121.9),
        Degrees::new(239.1),
        Degrees::new(205.7),
        JulianDate::J2000,
    ),
};

pub const MAKEMAKE: Planet = Planet {
    mass: Kilograms::new(3.1e21),
    radius: Kilometers::new(715.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(45.791),
        0.159,
        Degrees::new(29.01),
        Degrees::new(79.3),
        Degrees::new(286.1),
        Degrees::new(159.8),
        JulianDate::J2000,
    ),
};

pub const ERIS: Planet = Planet {
    mass: Kilograms::new(1.66e22),
    radius: Kilometers::new(1163.0),
    orbit: Orbit::new(
        AstronomicalUnits::new(67.864),
        0.44177,
        Degrees::new(44.04),
        Degrees::new(35.95),
        Degrees::new(151.231),
        Degrees::new(204.17),
        JulianDate::J2000,
    ),
};

pub const DWARF_PLANETS: &[&Planet] = &[&super::PLUTO, &CERES, &HAUMEA, &MAKEMAKE, &ERIS];

// -------------------------------------------------------------------------------------------------
//  Major moons (planet‑centric orbital elements approximated; not heliocentric)
// -------------------------------------------------------------------------------------------------

pub const IO: Satellite = Satellite::new_const(
    "Io",
    Kilograms::new(8.9319e22),
    Kilometers::new(1821.6),
    // Semi‑major axis 421 700 km
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.002823),
        eccentricity: 0.0041,
        inclination: Degrees::new(0.036),
        longitude_of_ascending_node: Degrees::new(43.977),
        argument_of_perihelion: Degrees::new(84.129),
        mean_anomaly_at_epoch: Degrees::new(171.016),
        epoch: JulianDate::J2000,
    },
);

pub const EUROPA: Satellite = Satellite::new_const(
    "Europa",
    Kilograms::new(4.7998e22),
    Kilometers::new(1560.8),
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.004485),
        eccentricity: 0.009,
        inclination: Degrees::new(0.465),
        longitude_of_ascending_node: Degrees::new(219.106),
        argument_of_perihelion: Degrees::new(88.970),
        mean_anomaly_at_epoch: Degrees::new(324.528),
        epoch: JulianDate::J2000,
    },
);

pub const GANYMEDE: Satellite = Satellite::new_const(
    "Ganymede",
    Kilograms::new(1.4819e23),
    Kilometers::new(2634.1),
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.007155),
        eccentricity: 0.0013,
        inclination: Degrees::new(0.177),
        longitude_of_ascending_node: Degrees::new(63.552),
        argument_of_perihelion: Degrees::new(192.417),
        mean_anomaly_at_epoch: Degrees::new(317.654),
        epoch: JulianDate::J2000,
    },
);

pub const CALLISTO: Satellite = Satellite::new_const(
    "Callisto",
    Kilograms::new(1.0759e23),
    Kilometers::new(2410.3),
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.012585),
        eccentricity: 0.0074,
        inclination: Degrees::new(0.192),
        longitude_of_ascending_node: Degrees::new(298.848),
        argument_of_perihelion: Degrees::new(52.643),
        mean_anomaly_at_epoch: Degrees::new(51.483),
        epoch: JulianDate::J2000,
    },
);

pub const TITAN: Satellite = Satellite::new_const(
    "Titan",
    Kilograms::new(1.3452e23),
    Kilometers::new(2574.73),
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.008167),
        eccentricity: 0.0288,
        inclination: Degrees::new(0.34854),
        longitude_of_ascending_node: Degrees::new(168.650),
        argument_of_perihelion: Degrees::new(186.585),
        mean_anomaly_at_epoch: Degrees::new(30.744),
        epoch: JulianDate::J2000,
    },
);

pub const TRITON: Satellite = Satellite::new_const(
    "Triton",
    Kilograms::new(2.14e22),
    Kilometers::new(1353.4),
    Orbit {
        semi_major_axis: AstronomicalUnits::new(0.002371),
        eccentricity: 0.000016,
        inclination: Degrees::new(156.865), // retrograde
        longitude_of_ascending_node: Degrees::new(216.732),
        argument_of_perihelion: Degrees::new(185.965),
        mean_anomaly_at_epoch: Degrees::new(144.960),
        epoch: JulianDate::J2000,
    },
);

pub const MAJOR_MOONS: &[&Satellite] = &[
    &super::MOON,
    &IO,
    &EUROPA,
    &GANYMEDE,
    &CALLISTO,
    &TITAN,
    &TRITON,
];

// -------------------------------------------------------------------------------------------------
//  Lagrange‑point helper type
// -------------------------------------------------------------------------------------------------

/// Simple position wrapper for a restricted‑three‑body Lagrange point.
#[derive(Debug, Clone)]
pub struct LagrangePoint {
    /// Designation, e.g. "Sun–Earth L1".
    pub name: &'static str,
    /// Primary–secondary pair the point belongs to (for information only).
    pub parent_system: &'static str,
    /// Heliocentric ecliptic coordinates referenced to J2000.0.
    pub position: Ecliptic<AstronomicalUnit>,
}

// For demonstration purposes we approximate L1/L2 as ±0.01 AU along the Sun–Earth line.
const SUN_EARTH_L1: LagrangePoint = LagrangePoint {
    name: "Sun–Earth L1",
    parent_system: "Sun–Earth",
    position: Ecliptic::<AstronomicalUnit>::new_raw(
        Degrees::new(0.0), // lat (polar)
        Degrees::new(0.0), // lon (azimuth)
        AstronomicalUnits::new(0.99),
    ),
};
const SUN_EARTH_L2: LagrangePoint = LagrangePoint {
    name: "Sun–Earth L2",
    parent_system: "Sun–Earth",
    position: Ecliptic::<AstronomicalUnit>::new_raw(
        Degrees::new(0.0),   // lat (polar)
        Degrees::new(180.0), // lon (azimuth)
        AstronomicalUnits::new(1.01),
    ),
};

pub const LAGRANGE_POINTS: &[&LagrangePoint] = &[&SUN_EARTH_L1, &SUN_EARTH_L2];

// -------------------------------------------------------------------------------------------------
//  Aggregate constant
// -------------------------------------------------------------------------------------------------

/// Bundles all canonical Solar‑System bodies into a single constant for
/// ergonomic iteration.
#[derive(Debug)]
pub struct SolarSystem<'a> {
    pub sun: &'a Star<'a>,
    pub planets: &'a [&'a Planet],
    pub dwarf_planets: &'a [&'a Planet],
    pub moons: &'a [&'a Satellite<'a>],
    pub lagrange_points: &'a [&'a LagrangePoint],
}

pub const PLANETS: &[&Planet] = &[
    &super::MERCURY,
    &super::VENUS,
    &super::EARTH,
    &super::MARS,
    &super::JUPITER,
    &super::SATURN,
    &super::URANUS,
    &super::NEPTUNE,
];

pub const SOLAR_SYSTEM: SolarSystem<'static> = SolarSystem {
    sun: &super::SUN,
    planets: PLANETS,
    dwarf_planets: DWARF_PLANETS,
    moons: MAJOR_MOONS,
    lagrange_points: LAGRANGE_POINTS,
};
