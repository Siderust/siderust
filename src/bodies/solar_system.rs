//! # Solar System
//!
//! This module exposes a set of **public constants** representing the Sun, the eight major planets,
//! Earth's Moon, and Pluto. The physical and orbital parameters are sourced from the
//! _NASA Planetary Fact Sheet_\[1]\[2]. Values are given at the **J2000.0** epoch unless otherwise noted.
//!
//! Each constant can be used for astronomical calculations or as reference data. Each entry includes:
//!
//! * **Mass** (kg)
//! * **Mean radius** (km)
//! * **Keplerian orbital elements** (see [Orbit]):
//!   * `a` – semi-major axis [AstronomicalUnit]
//!   * `e` – eccentricity (unitless)
//!   * `i` – inclination [Degrees]
//!   * `Ω` – longitude of ascending node [Degrees]
//!   * `ω` – argument of perihelion [Degrees]
//!   * `M₀` – mean anomaly at epoch [Degrees]
//!
//! ---
//! ## References
//! 1. NASA – Planetary Fact Sheet: <https://nssdc.gsfc.nasa.gov/planetary/factsheet/>
//! 2. Williams, D. R. (2024). _Planetary Fact Sheet – Metric_. NASA Goddard Space Flight Center.
//!
//! ```text
//! Abbreviated notation used for orbital elements :
//!   a  — semi-major axis         Ω — longitude of ascending node
//!   e  — eccentricity            ω — argument of perihelion
//!   i  — inclination             M₀ — mean anomaly at epoch
//! ```

use crate::units::*;
use crate::astro::orbit::Orbit;
use crate::coordinates::{centers::Geocentric, frames::Equatorial, spherical::Position};
use crate::targets::Target;

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
/// | Radius          | 696,340 km          |
/// | Luminosity      | 3.828 × 10²⁶ W       |
/// | Right Ascension | 18h 44m 48s (J2000) |
/// | Declination     | −23° 00′ 00″ (J2000)|
/// | LengthUnit        | 1 AU (~0.0000158 ly)|
pub const SUN: super::Star<'static> = super::Star::new_const(
    "Sun",
    AU.to_light_year(),
    SOLAR_MASS,
    SOLAR_RADIUS,
    SOLAR_LUMINOSITY,
    Target::<Position::<Geocentric, Equatorial, LightYear>>::new_static(
        Position::<Geocentric, Equatorial, LightYear>::new_const(
            Degrees::from_hms(18, 44, 48.0), // Aprox at J2000
            Degrees::from_hms(-23, 0, 0.0), // Aprox at J2000
            AU.to_light_year(),
        ),
        JulianDay::J2000,
    ),
);

/// **Mercury** – innermost planet of the Solar System.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 3.3011 × 10²³ kg  |
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
        AU::new(0.38709893),
        0.20563069,
        Degrees::new(7.00487),
        Degrees::new(48.33167),
        Degrees::new(29.12478),
        Degrees::new(174.79439),
        JulianDay::J2000,
    ),
};

/// **Venus** – second planet from the Sun.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 4.8675 × 10²⁴ kg  |
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
        AU::new(0.72333199),
        0.00677323,
        Degrees::new(3.39471),
        Degrees::new(76.68069),
        Degrees::new(54.85229),
        Degrees::new(50.44675),
        JulianDay::J2000,
    ),
};

/// **Earth** – third planet from the Sun and our home.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 5.97237 × 10²⁴ kg |
/// | Radius    | 6,371.0 km       |
/// | a         | 1.000000 AU      |
/// | e         | 0.016710         |
/// | i         | 0.00005°         |
/// | Ω         | -11.26064°       |
/// | ω         | 114.20783°       |
/// | M₀        | 357.51716°       |
pub const EARTH: super::Planet = super::Planet {
    mass: Kilograms::new(5.97237e24),
    radius: Kilometers::new(6371.0),
    orbit: Orbit::new(
        AU::new(1.00000011),
        0.01671022,
        Degrees::new(0.00005),
        Degrees::new(-11.26064),
        Degrees::new(114.20783),
        Degrees::new(357.51716),
        JulianDay::J2000,
    ),
};

/// **Moon** – natural satellite of Earth.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 7.346 × 10²² kg   |
/// | Radius    | 1,737.4 km       |
/// | a         | ~0.00257 AU      |
/// | e         | 0.0549           |
/// | i         | 5.145°           |
/// | Ω         | 125.08°          |
/// | ω         | 318.15°          |
/// | M₀        | 135.27°          |
pub const MOON: super::Satelite = super::Satelite::new_const(
    "Moon",
    Kilograms::new(7.346e22),
    Kilometers::new(1_737.4),

    Orbit {
        // 384 400 km → AU
        semi_major_axis: Kilometers::new(384_400.0).to_au(),
        eccentricity: 0.054_9,
        inclination: Degrees::new(5.145),
        longitude_of_ascending_node: Degrees::new(125.08),
        argument_of_perihelion: Degrees::new(318.15),
        mean_anomaly_at_epoch: Degrees::new(135.27),
        epoch: JulianDay::J2000,
    },
);

/// **Mars** – the red planet, fourth from the Sun.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 6.4171 × 10²³ kg  |
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
        AU::new(1.52366231),
        0.09341233,
        Degrees::new(1.85061),
        Degrees::new(49.57854),
        Degrees::new(286.46230),
        Degrees::new(19.41248),
        JulianDay::J2000,
    ),
};

/// **Jupiter** – the largest planet in the Solar System.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 1.8982 × 10²⁷ kg  |
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
        AU::new(5.20336301),
        0.04839266,
        Degrees::new(1.30530),
        Degrees::new(100.55615),
        Degrees::new(274.19770),
        Degrees::new(19.65053),
        JulianDay::J2000,
    ),
};

/// **Saturn** – known for its prominent ring system.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 5.6834 × 10²⁶ kg  |
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
        AU::new(9.53707032),
        0.05415060,
        Degrees::new(2.48446),
        Degrees::new(113.71504),
        Degrees::new(338.71690),
        Degrees::new(317.51238),
        JulianDay::J2000,
    ),
};

/// **Uranus** – icy giant with extreme axial tilt.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 8.6810 × 10²⁵ kg  |
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
        AU::new(19.19126393),
        0.04716771,
        Degrees::new(0.76986),
        Degrees::new(74.22988),
        Degrees::new(96.73436),
        Degrees::new(142.26794),
        JulianDay::J2000,
    ),
};

/// **Neptune** – outermost giant planet.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 1.02409 × 10²⁶ kg |
/// | Radius    | 24,622.0 km      |
/// | a         | 30.068963 AU     |
/// | e         | 0.008586         |
/// | i         | 1.76917°         |
/// | Ω         | 131.72169°       |
/// | ω         | 273.24966°       |
/// | M₀        | 259.90868°       |
pub const NEPTUNE: super::Planet = super::Planet {
    mass: Kilograms::new(1.02409e26),
    radius: Kilometers::new(24622.0),
    orbit: Orbit::new(
        AU::new(30.06896348),
        0.00858587,
        Degrees::new(1.76917),
        Degrees::new(131.72169),
        Degrees::new(273.24966),
        Degrees::new(259.90868),
        JulianDay::J2000,
    ),
};

/// **Pluto** – dwarf planet once considered the ninth planet.
///
/// | Parameter | Value             |
/// |-----------|------------------|
/// | Mass      | 1.303 × 10²² kg   |
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
        AU::new(39.48168677),
        0.24880766,
        Degrees::new(17.14175),
        Degrees::new(110.30347),
        Degrees::new(113.76329),
        Degrees::new(14.86205),
        JulianDay::J2000,
    ),
};
