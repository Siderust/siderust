// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar System
//!
//! ## Scientific scope
//!
//! Our Solar System consists of the Sun, eight planets, five IAU-recognised
//! dwarf planets, more than 200 moons, millions of asteroids, and billions
//! of comets. Physical parameters (mass, radius, Bond albedo) and Keplerian
//! heliocentric orbital elements at the J2000.0 epoch are tabulated here.
//! Orbital elements are osculating values from the JPL planetary ephemeris;
//! they capture mean orbital behaviour well for visualisation and coarse
//! ephemeris but should not be used for precision close-approach analysis
//! (use VSOP87 / ELP-2000 via [`crate::calculus`] for that). Bond albedo
//! values are IAU/NASA consensus; they are typed as [`Albedos`] to prevent
//! unit-confusion with geometric albedo.
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - Marker structs: [`Sun`], [`Mercury`], [`Venus`], [`Earth`], [`Moon`],
//!   [`Mars`], [`Jupiter`], [`Saturn`], [`Uranus`], [`Neptune`], [`Pluto`].
//! - IAU rotation parameters via the [`HasIauRotation`] trait + `ROTATION`
//!   associated constants for each body.
//! - Planet constants ([`MERCURY`] … [`NEPTUNE`], [`PLUTO`]) as typed
//!   [`Planet`] values with mass ([`Kilograms`]), radius ([`Kilometers`]),
//!   Bond albedo ([`Albedos`]), and [`KeplerianOrbit`].
//! - Dwarf planet constants ([`CERES`], [`HAUMEA`], [`MAKEMAKE`], [`ERIS`]).
//! - Moon constants ([`MOON`], [`IO`], [`EUROPA`], [`GANYMEDE`],
//!   [`CALLISTO`], [`TITAN`], [`TRITON`]) as [`Satellite`] values.
//! - [`LagrangePoint`] type + [`LAGRANGE_POINTS`] slice.
//! - [`SolarSystem`] aggregate struct + [`SOLAR_SYSTEM`] constant.
//! - [`Trackable`] implementations for every planet and the Moon.
//!
//! ## References
//!
//! - Williams, D. R. (2024). *Planetary Fact Sheet – Metric*. NASA GSFC.
//!   <https://nssdc.gsfc.nasa.gov/planetary/factsheet/>
//! - IAU Working Group on Cartographic Coordinates and Rotational Elements
//!   (2015). *Celestial Mechanics and Dynamical Astronomy* 130, 22.
//!   doi:10.1007/s10569-017-9805-5
//! - Seidelmann, P. K. (Ed.) (1992). *Explanatory Supplement to the
//!   Astronomical Almanac*. University Science Books.

use super::{Planet, Satellite, Star};
use crate::astro::orbit::KeplerianOrbit;
use crate::astro::{HasIauRotation, IauRotationParams};
use crate::coordinates::spherical::position::{EclipticMeanJ2000, EquatorialMeanJ2000};
use crate::qtty::length::nominal::RSUN;
use crate::qtty::*;
use crate::targets::CoordinateWithPM;
use crate::time::JulianDate;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Sun;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Mercury;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Venus;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Earth;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Moon;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Mars;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Jupiter;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Saturn;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Uranus;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Neptune;
#[derive(Debug, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Pluto;

// =============================================================================
// Planetary IAU Rotation Parameters
// =============================================================================

macro_rules! define_body_rotations {
    (
        $(
            $rotation_const:ident => $body:ty {
                alpha0_deg: $alpha0_deg:expr,
                alpha0_rate: $alpha0_rate:expr,
                delta0_deg: $delta0_deg:expr,
                delta0_rate: $delta0_rate:expr,
                w0_deg: $w0_deg:expr,
                w_rate: $w_rate:expr $(,)?
            }
        ),+ $(,)?
    ) => {
        $(
            pub const $rotation_const: IauRotationParams = IauRotationParams {
                alpha0_deg: $alpha0_deg,
                alpha0_rate: $alpha0_rate,
                delta0_deg: $delta0_deg,
                delta0_rate: $delta0_rate,
                w0_deg: $w0_deg,
                w_rate: $w_rate,
            };

            impl HasIauRotation for $body {
                const ROTATION: IauRotationParams = $rotation_const;
            }

            impl $body {
                pub const ROTATION: IauRotationParams = $rotation_const;
            }
        )+
    };
}

define_body_rotations!(
    MERCURY_ROTATION => Mercury {
        alpha0_deg: Degrees::new(281.0103),
        alpha0_rate: Degrees::new(-0.0328),
        delta0_deg: Degrees::new(61.4155),
        delta0_rate: Degrees::new(-0.0049),
        w0_deg: Degrees::new(329.5988),
        w_rate: Degrees::new(6.1385108),
    },
    VENUS_ROTATION => Venus {
        alpha0_deg: Degrees::new(272.76),
        alpha0_rate: Degrees::new(0.0),
        delta0_deg: Degrees::new(67.16),
        delta0_rate: Degrees::new(0.0),
        w0_deg: Degrees::new(160.20),
        w_rate: Degrees::new(-1.4813688),
    },
    MARS_ROTATION => Mars {
        alpha0_deg: Degrees::new(317.269),
        alpha0_rate: Degrees::new(-0.10927),
        delta0_deg: Degrees::new(54.432),
        delta0_rate: Degrees::new(-0.05827),
        w0_deg: Degrees::new(176.049),
        w_rate: Degrees::new(350.891982443),
    },
    MOON_ROTATION => Moon {
        alpha0_deg: Degrees::new(269.9949),
        alpha0_rate: Degrees::new(0.0031),
        delta0_deg: Degrees::new(66.5392),
        delta0_rate: Degrees::new(0.0130),
        w0_deg: Degrees::new(38.3213),
        w_rate: Degrees::new(13.17635815),
    },
    JUPITER_ROTATION => Jupiter {
        alpha0_deg: Degrees::new(268.057),
        alpha0_rate: Degrees::new(-0.006),
        delta0_deg: Degrees::new(64.495),
        delta0_rate: Degrees::new(0.002),
        w0_deg: Degrees::new(284.95),
        w_rate: Degrees::new(870.5360000),
    },
    SATURN_ROTATION => Saturn {
        alpha0_deg: Degrees::new(40.589),
        alpha0_rate: Degrees::new(-0.036),
        delta0_deg: Degrees::new(83.537),
        delta0_rate: Degrees::new(-0.004),
        w0_deg: Degrees::new(38.90),
        w_rate: Degrees::new(810.7939024),
    },
    URANUS_ROTATION => Uranus {
        alpha0_deg: Degrees::new(257.311),
        alpha0_rate: Degrees::new(0.0),
        delta0_deg: Degrees::new(-15.175),
        delta0_rate: Degrees::new(0.0),
        w0_deg: Degrees::new(203.81),
        w_rate: Degrees::new(-501.1600928),
    },
    NEPTUNE_ROTATION => Neptune {
        alpha0_deg: Degrees::new(299.36),
        alpha0_rate: Degrees::new(0.0), // periodic term handled separately
        delta0_deg: Degrees::new(43.46),
        delta0_rate: Degrees::new(0.0), // periodic term handled separately
        w0_deg: Degrees::new(249.978),
        w_rate: Degrees::new(541.1397757),
    },
    PLUTO_ROTATION => Pluto {
        alpha0_deg: Degrees::new(132.993),
        alpha0_rate: Degrees::new(0.0),
        delta0_deg: Degrees::new(-6.163),
        delta0_rate: Degrees::new(0.0),
        w0_deg: Degrees::new(302.695),
        w_rate: Degrees::new(56.3625225),
    },
);

/// **Sun** – the central star of the Solar System.
///
/// | Parameter       | Value                |
/// |-----------------|----------------------|
/// | Mass            | 1.9885 × 10³⁰ kg     |
/// | Radius          | 695,700 km (IAU 2015 nominal) |
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
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_unchecked(
            HourAngles::from_hms(-23, 0, 0.0).to_const::<Degree>(), // Dec (polar) Approx at J2000
            HourAngles::from_hms(18, 44, 48.0).to_const::<Degree>(), // RA (azimuth) Approx at J2000
            LightYears::new(1.58125e-5), // 1 AstronomicalUnits in LightYears
        ),
        crate::J2000,
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(0.38709893),
        0.20563069,
        Degrees::new(7.00487),
        Degrees::new(48.33167),
        Degrees::new(29.12478),
        Degrees::new(174.79439),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.068)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(0.72333199),
        0.00677323,
        Degrees::new(3.39471),
        Degrees::new(76.68069),
        Degrees::new(54.85229),
        Degrees::new(50.44675),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.76)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(1.00000011),
        0.01671022,
        Degrees::new(0.00005),
        Degrees::new(-11.26064),
        Degrees::new(114.20783),
        Degrees::new(357.51716),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.306)),
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
    KeplerianOrbit::new(
        AstronomicalUnits::new(2.566881e-6),
        0.054_9,
        Degrees::new(5.145),
        Degrees::new(125.08),
        Degrees::new(318.15),
        Degrees::new(135.27),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.12));

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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(1.52366231),
        0.09341233,
        Degrees::new(1.85061),
        Degrees::new(49.57854),
        Degrees::new(286.46230),
        Degrees::new(19.41248),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.25)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(5.20336301),
        0.04839266,
        Degrees::new(1.30530),
        Degrees::new(100.55615),
        Degrees::new(274.19770),
        Degrees::new(19.65053),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.503)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(9.53707032),
        0.05415060,
        Degrees::new(2.48446),
        Degrees::new(113.71504),
        Degrees::new(338.71690),
        Degrees::new(317.51238),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.342)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(19.19126393),
        0.04716771,
        Degrees::new(0.76986),
        Degrees::new(74.22988),
        Degrees::new(96.73436),
        Degrees::new(142.26794),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.3)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(30.06896348),
        0.00858587,
        Degrees::new(1.76917),
        Degrees::new(131.72169),
        Degrees::new(273.24966),
        Degrees::new(259.90868),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.29)),
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
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(39.48168677),
        0.24880766,
        Degrees::new(17.14175),
        Degrees::new(110.30347),
        Degrees::new(113.76329),
        Degrees::new(14.86205),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.52)),
};

// -------------------------------------------------------------------------------------------------
//  Dwarf planets
// -------------------------------------------------------------------------------------------------

pub const CERES: Planet = Planet {
    mass: Kilograms::new(9.393e20),
    radius: Kilometers::new(473.0),
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(2.7675),
        0.0758,
        Degrees::new(10.5941),
        Degrees::new(80.3055),
        Degrees::new(73.5977),
        Degrees::new(95.9892),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.09)),
};

pub const HAUMEA: Planet = Planet {
    mass: Kilograms::new(4.006e21),
    radius: Kilometers::new(620.0),
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(43.218),
        0.195,
        Degrees::new(28.19),
        Degrees::new(121.9),
        Degrees::new(239.1),
        Degrees::new(205.7),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.8)),
};

pub const MAKEMAKE: Planet = Planet {
    mass: Kilograms::new(3.1e21),
    radius: Kilometers::new(715.0),
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(45.791),
        0.159,
        Degrees::new(29.01),
        Degrees::new(79.3),
        Degrees::new(286.1),
        Degrees::new(159.8),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.81)),
};

pub const ERIS: Planet = Planet {
    mass: Kilograms::new(1.66e22),
    radius: Kilometers::new(1163.0),
    orbit: KeplerianOrbit::new(
        AstronomicalUnits::new(67.864),
        0.44177,
        Degrees::new(44.04),
        Degrees::new(35.95),
        Degrees::new(151.231),
        Degrees::new(204.17),
        crate::J2000,
    ),
    albedo: Some(Albedos::new(0.96)),
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
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.002823),
        0.0041,
        Degrees::new(0.036),
        Degrees::new(43.977),
        Degrees::new(84.129),
        Degrees::new(171.016),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.62));

pub const EUROPA: Satellite = Satellite::new_const(
    "Europa",
    Kilograms::new(4.7998e22),
    Kilometers::new(1560.8),
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.004485),
        0.009,
        Degrees::new(0.465),
        Degrees::new(219.106),
        Degrees::new(88.970),
        Degrees::new(324.528),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.67));

pub const GANYMEDE: Satellite = Satellite::new_const(
    "Ganymede",
    Kilograms::new(1.4819e23),
    Kilometers::new(2634.1),
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.007155),
        0.0013,
        Degrees::new(0.177),
        Degrees::new(63.552),
        Degrees::new(192.417),
        Degrees::new(317.654),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.43));

pub const CALLISTO: Satellite = Satellite::new_const(
    "Callisto",
    Kilograms::new(1.0759e23),
    Kilometers::new(2410.3),
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.012585),
        0.0074,
        Degrees::new(0.192),
        Degrees::new(298.848),
        Degrees::new(52.643),
        Degrees::new(51.483),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.17));

pub const TITAN: Satellite = Satellite::new_const(
    "Titan",
    Kilograms::new(1.3452e23),
    Kilometers::new(2574.73),
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.008167),
        0.0288,
        Degrees::new(0.34854),
        Degrees::new(168.650),
        Degrees::new(186.585),
        Degrees::new(30.744),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.22));

pub const TRITON: Satellite = Satellite::new_const(
    "Triton",
    Kilograms::new(2.14e22),
    Kilometers::new(1353.4),
    KeplerianOrbit::new(
        AstronomicalUnits::new(0.002371),
        0.000016,
        Degrees::new(156.865),
        Degrees::new(216.732),
        Degrees::new(185.965),
        Degrees::new(144.960),
        crate::J2000,
    ),
)
.with_albedo(Albedos::new(0.76));

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
    pub position: EclipticMeanJ2000<AstronomicalUnit>,
}

// For demonstration purposes we approximate L1/L2 as ±0.01 AU along the Sun–Earth line.
const SUN_EARTH_L1: LagrangePoint = LagrangePoint {
    name: "Sun–Earth L1",
    parent_system: "Sun–Earth",
    position: EclipticMeanJ2000::<AstronomicalUnit>::new_unchecked(
        Degrees::new(0.0), // lat (polar)
        Degrees::new(0.0), // lon (azimuth)
        AstronomicalUnits::new(0.99),
    ),
};
const SUN_EARTH_L2: LagrangePoint = LagrangePoint {
    name: "Sun–Earth L2",
    parent_system: "Sun–Earth",
    position: EclipticMeanJ2000::<AstronomicalUnit>::new_unchecked(
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

#[cfg(test)]
mod tests {
    use super::*;
    use affn::frames::{
        JupiterSystemIII, MarsFixed, MercuryFixed, MoonPrincipalAxes, NeptuneFixed, PlutoFixed,
        ReferenceFrame, SaturnFixed, SphericalNaming, UranusFixed, VenusFixed,
    };

    #[test]
    fn planetary_frame_names() {
        assert_eq!(MercuryFixed::frame_name(), "MercuryFixed");
        assert_eq!(VenusFixed::frame_name(), "VenusFixed");
        assert_eq!(MarsFixed::frame_name(), "MarsFixed");
        assert_eq!(MoonPrincipalAxes::frame_name(), "MoonPrincipalAxes");
        assert_eq!(JupiterSystemIII::frame_name(), "JupiterSystemIII");
        assert_eq!(SaturnFixed::frame_name(), "SaturnFixed");
        assert_eq!(UranusFixed::frame_name(), "UranusFixed");
        assert_eq!(NeptuneFixed::frame_name(), "NeptuneFixed");
        assert_eq!(PlutoFixed::frame_name(), "PlutoFixed");
    }

    #[test]
    fn spherical_naming_is_planetocentric() {
        assert_eq!(MercuryFixed::polar_name(), "lat");
        assert_eq!(MercuryFixed::azimuth_name(), "lon");
        assert_eq!(MercuryFixed::distance_name(), "radius");
        assert_eq!(MarsFixed::polar_name(), "lat");
        assert_eq!(MarsFixed::azimuth_name(), "lon");
    }

    #[test]
    fn rotation_params_at_j2000() {
        let p = &MARS_ROTATION;
        assert!((p.alpha0(crate::J2000).value() - 317.269).abs() < 1e-10);
        assert!((p.delta0(crate::J2000).value() - 54.432).abs() < 1e-10);
        assert!((p.w(crate::J2000).value() - 176.049).abs() < 1e-10);
    }

    #[test]
    fn rotation_params_rate() {
        let p = &MARS_ROTATION;
        let alpha = p.alpha0(crate::time::JulianDate::new((crate::J2000.raw() + Days::new(36_525.0)).value()));
        assert!((alpha.value() - (317.269 - 0.10927)).abs() < 1e-10);
        let w = p.w(crate::time::JulianDate::new((
            crate::J2000.raw() + crate::qtty::Days::new(1.0),
        ).value()));
        assert!((w.value() - (176.049 + 350.891982443)).abs() < 1e-8);
    }

    #[test]
    fn body_associated_rotation_constants_match_catalog() {
        assert_eq!(Mars::ROTATION.alpha0_deg, MARS_ROTATION.alpha0_deg);
        assert_eq!(Jupiter::ROTATION.w_rate, JUPITER_ROTATION.w_rate);
        assert_eq!(Moon::ROTATION.delta0_deg, MOON_ROTATION.delta0_deg);
    }
}

// =============================================================================
// Trackable implementations
//
// These live here (not in targets/trackable.rs) to keep targets free of any
// dependency on bodies, breaking the targets ↔ bodies cycle.
// =============================================================================

use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Geocentric},
    frames::EclipticMeanJ2000 as EclJ2000,
};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::targets::Trackable;

/// Barycentric ecliptic position (AU), the natural output of VSOP87e.
type BaryEclPos = CoordinateWithPM<Position<Barycentric, EclJ2000, AstronomicalUnit>>;

/// Geocentric ecliptic position (km), the natural output of ELP2000.
type GeoEclPos = Position<Geocentric, EclJ2000, Kilometer>;

macro_rules! impl_trackable_vsop87 {
    ($($Planet:ident),+ $(,)?) => {
        $(
            impl Trackable for $Planet {
                type Coords = BaryEclPos;

                #[inline]
                fn track(&self, jd: JulianDate) -> Self::Coords {
                    CoordinateWithPM::new_static($Planet::vsop87e(jd), jd)
                }
            }
        )+
    };
}

impl_trackable_vsop87!(Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune);

impl Trackable for Moon {
    type Coords = GeoEclPos;

    #[inline]
    fn track(&self, jd: JulianDate) -> Self::Coords {
        Moon::get_geo_position::<Kilometer>(jd)
    }
}

#[cfg(test)]
mod trackable_tests {
    use super::*;

    #[test]
    fn star_produces_icrs_direction() {
        use crate::bodies::catalog;
        let sirius = &catalog::SIRIUS;
        let dir = <crate::bodies::Star<'_> as Trackable>::track(sirius, crate::J2000);
        let ra_diff = dir.ra() - Degrees::new(101.287);
        assert!(ra_diff <= Degrees::new(0.01) && ra_diff >= -Degrees::new(0.01));
        let dec_diff = dir.dec() + Degrees::new(16.716);
        assert!(dec_diff <= Degrees::new(0.01) && dec_diff >= -Degrees::new(0.01));
    }

    #[test]
    fn sun_produces_barycentric_position() {
        let pos = Sun.track(crate::J2000);
        let dist = pos.position.distance();
        assert!(
            dist < AstronomicalUnits::new(0.02),
            "Sun should be near SSB, got {dist}"
        );
    }

    #[test]
    fn earth_changes_with_time() {
        let p1 = Earth.track(crate::J2000);
        let p2 = Earth.track(crate::time::JulianDate::new((crate::J2000.raw() + Days::new(182.625)).value()));
        let sep = p1.position.distance_to(&p2.position);
        assert!(
            sep > AstronomicalUnits::new(1.0),
            "Half-year separation should be > 1 AU, got {sep}"
        );
    }

    #[test]
    fn moon_produces_geocentric_position() {
        let pos = Moon.track(crate::J2000);
        let dist = pos.distance();
        assert!(
            dist >= Kilometers::new(350_000.0) && dist <= Kilometers::new(410_000.0),
            "Moon distance should be ~384 400 km, got {dist}"
        );
    }

    #[test]
    fn planets_produce_nonzero_positions() {
        let dist = Mars.track(crate::J2000).position.distance();
        assert!(
            dist > AstronomicalUnits::new(1.0),
            "Mars should be > 1 AU from SSB, got {dist}"
        );
    }
}
