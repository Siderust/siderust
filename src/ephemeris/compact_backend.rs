// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use super::{
    native_frames, AuPerDay, Ephemeris, EphemerisError, MajorPlanet, PlanetEphemerisError,
    PlanetState,
};
use crate::archive::{elp_compact, vsop_compact};
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::qtty::{AstronomicalUnit, AstronomicalUnits, Kilometer};
use crate::time::JulianDate;
use tempoch::TDB;

type BarycentricAu = Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>;
type HeliocentricAu = Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;
type LunarPosition = Position<Geocentric, EclipticMeanJ2000, Kilometer>;
type VsopEvaluator = fn(f64) -> Result<[AstronomicalUnits; 3], vsop_compact::CompactEphemerisError>;

/// Offline VSOP87E and ELP2000-82B backend with truncated coefficients.
///
/// Use the fallible `Ephemeris::try_*` methods for input epochs. The legacy
/// infallible trait methods panic outside the supported interval.
/// VSOP87E uses the original VSOP87 body definitions. The outer planet
/// series approximate the planet directions at the validated accuracy;
/// their independent DE440s reference uses planetary-system barycenters.
#[derive(Debug, Clone, Copy, Default)]
pub struct CompactEphemeris;

impl CompactEphemeris {
    /// First supported Julian date in TDB (inclusive).
    pub const FIRST_JD_TDB: f64 = vsop_compact::FIRST_JD_TDB;
    /// End of supported coverage in TDB (exclusive).
    pub const END_JD_TDB: f64 = vsop_compact::END_JD_TDB;
}

fn planetary_position(evaluate: VsopEvaluator, tdb: f64) -> Result<BarycentricAu, EphemerisError> {
    let [x, y, z] = evaluate(tdb).map_err(|_| EphemerisError::OutOfRange {
        jd: tdb,
        start_jd: vsop_compact::FIRST_JD_TDB,
        end_jd: vsop_compact::END_JD_TDB,
    })?;
    Ok(native_frames::vsop_position(x, y, z))
}

fn heliocentric(body: BarycentricAu, sun: BarycentricAu) -> HeliocentricAu {
    Position::new(body.x() - sun.x(), body.y() - sun.y(), body.z() - sun.z())
}

impl Ephemeris for CompactEphemeris {
    fn try_planet_barycentric(
        planet: MajorPlanet,
        jd: JulianDate,
    ) -> Result<BarycentricAu, PlanetEphemerisError> {
        let evaluate: VsopEvaluator = match planet {
            MajorPlanet::Mercury => vsop_compact::try_mercury,
            MajorPlanet::Venus => vsop_compact::try_venus,
            MajorPlanet::Mars => vsop_compact::try_mars,
            MajorPlanet::Jupiter => vsop_compact::try_jupiter,
            MajorPlanet::Saturn => vsop_compact::try_saturn,
            MajorPlanet::Uranus => vsop_compact::try_uranus,
            MajorPlanet::Neptune => vsop_compact::try_neptune,
        };
        Ok(planetary_position(evaluate, jd.to_scale::<TDB>().value())?)
    }

    fn try_sun_barycentric(jd: JulianDate) -> Result<BarycentricAu, EphemerisError> {
        planetary_position(vsop_compact::try_sun, jd.to_scale::<TDB>().value())
    }

    fn sun_barycentric(jd: JulianDate) -> BarycentricAu {
        Self::try_sun_barycentric(jd).expect("compact Sun position unavailable")
    }

    fn try_earth_barycentric(jd: JulianDate) -> Result<BarycentricAu, EphemerisError> {
        planetary_position(vsop_compact::try_earth, jd.to_scale::<TDB>().value())
    }

    fn earth_barycentric(jd: JulianDate) -> BarycentricAu {
        Self::try_earth_barycentric(jd).expect("compact Earth position unavailable")
    }

    fn try_earth_heliocentric(jd: JulianDate) -> Result<HeliocentricAu, EphemerisError> {
        Ok(heliocentric(
            Self::try_earth_barycentric(jd)?,
            Self::try_sun_barycentric(jd)?,
        ))
    }

    fn earth_heliocentric(jd: JulianDate) -> HeliocentricAu {
        Self::try_earth_heliocentric(jd).expect("compact Earth position unavailable")
    }

    fn try_earth_barycentric_velocity(
        jd: JulianDate,
    ) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
        let tdb = jd.to_scale::<TDB>().value();
        let earth = |t| planetary_position(vsop_compact::try_earth, t);
        // A binary-exact step of 675 seconds keeps all offsets exact at
        // these Julian dates. Switch to a second-order one-sided stencil
        // when the symmetric stencil would leave the validity interval.
        const STEP: f64 = 1.0 / 128.0;
        let components = |position: BarycentricAu| position.as_array().map(|axis| axis.value());
        let derivative: [f64; 3] =
            if tdb - STEP >= Self::FIRST_JD_TDB && tdb + STEP < Self::END_JD_TDB {
                let a = components(earth(tdb - STEP)?);
                let b = components(earth(tdb + STEP)?);
                core::array::from_fn(|i| (b[i] - a[i]) / (2.0 * STEP))
            } else {
                let step = if tdb - STEP < Self::FIRST_JD_TDB {
                    STEP
                } else {
                    -STEP
                };
                let a = components(earth(tdb)?);
                let b = components(earth(tdb + step)?);
                let c = components(earth(tdb + 2.0 * step)?);
                core::array::from_fn(|i| (4.0 * (b[i] - a[i]) - (c[i] - a[i])) / (2.0 * step))
            };
        let [x, y, z] = derivative;
        Ok(Velocity::new(x, y, z))
    }

    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        Self::try_earth_barycentric_velocity(jd).expect("compact Earth velocity unavailable")
    }

    fn try_moon_geocentric(jd: JulianDate) -> Result<LunarPosition, EphemerisError> {
        let tdb = jd.to_scale::<TDB>().value();
        let xyz =
            elp_compact::try_moon_geocentric(tdb).map_err(|_| EphemerisError::OutOfRange {
                jd: tdb,
                start_jd: elp_compact::FIRST_JD_TDB,
                end_jd: elp_compact::END_JD_TDB,
            })?;
        let [x, y, z] = native_frames::elp_to_canonical().apply_array(xyz.map(|axis| axis.value()));
        Ok(Position::new(x, y, z))
    }

    fn moon_geocentric(jd: JulianDate) -> LunarPosition {
        Self::try_moon_geocentric(jd).expect("compact Moon position unavailable")
    }
}

macro_rules! impl_planet_state {
    ($($body:ident => $evaluate:ident),+ $(,)?) => {$(
        impl PlanetState<crate::bodies::$body> for CompactEphemeris {
            fn try_heliocentric(jd: JulianDate) -> Result<HeliocentricAu, EphemerisError> {
                let body = planetary_position(vsop_compact::$evaluate, jd.to_scale::<TDB>().value())?;
                Ok(heliocentric(body, Self::try_sun_barycentric(jd)?))
            }

            fn heliocentric(jd: JulianDate) -> HeliocentricAu {
                <Self as PlanetState<crate::bodies::$body>>::try_heliocentric(jd).expect("compact planet position unavailable")
            }
        }
    )+};
}

impl_planet_state!(
    Mercury => try_mercury,
    Venus => try_venus,
    Mars => try_mars,
    Jupiter => try_jupiter,
    Saturn => try_saturn,
    Uranus => try_uranus,
    Neptune => try_neptune,
);
