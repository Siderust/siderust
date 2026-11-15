// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(clippy::print_stdout)]

//! Offline compact Sun, Moon, and Mars observations.

use siderust::astro::{eop::NullEop, nutation::Iau2000B};
use siderust::coordinates::{
    cartesian::Position,
    centers::{Geocentric, Geodetic, Topocentric},
    frames::{EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialTrueOfDate, ECEF},
    transform::{
        horizontal::project_position_with, AstroContext, PositionAstroExt, TransformCenter,
    },
};
use siderust::ephemeris::{CompactEphemeris, Ephemeris, EphemerisError, PlanetState};
use siderust::event::lunar::{
    phase::moon_phase_from_topocentric, position_angle::moon_bright_limb_from_topocentric,
};
use siderust::qtty::{AstronomicalUnit, Degrees, Meters};
use siderust::{bodies::Mars, JulianDate};

fn main() -> Result<(), EphemerisError> {
    let jd = JulianDate::new(2_461_355.5); // TT
    let site = Geodetic::<ECEF>::new(
        Degrees::new(-0.1275),
        Degrees::new(51.507222),
        Meters::new(20.0),
    );
    // Explicit offline approximation: UT1 = UTC, zero polar motion.
    let ctx = AstroContext::<CompactEphemeris, NullEop>::default();
    let model = ctx.with_model::<Iau2000B>();

    let moon = CompactEphemeris::try_moon_geocentric(jd)?;
    let earth = CompactEphemeris::try_earth_heliocentric(jd)?;
    let mars = <CompactEphemeris as PlanetState<Mars>>::try_heliocentric(jd)?;
    let sun = Position::<Geocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        -earth.x(),
        -earth.y(),
        -earth.z(),
    );
    let mars = Position::<Geocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        mars.x() - earth.x(),
        mars.y() - earth.y(),
        mars.z() - earth.z(),
    );
    let moon: Position<Topocentric, _, _> = moon
        .to_frame_with::<EquatorialMeanJ2000, _>(&jd, &model)
        .to_center_with(site, jd, &model);
    let sun: Position<Topocentric, _, _> = sun
        .to_frame_with::<EquatorialMeanJ2000, _>(&jd, &model)
        .to_center_with(site, jd, &model);
    let mars: Position<Topocentric, _, _> = mars
        .to_frame_with::<EquatorialMeanJ2000, _>(&jd, &model)
        .to_center_with(site, jd, &model);
    let moon = moon
        .to_frame_with::<EquatorialTrueOfDate, _>(&jd, &model)
        .to_spherical();
    let sun = sun
        .to_frame_with::<EquatorialTrueOfDate, _>(&jd, &model)
        .to_spherical();
    let mars = mars
        .to_frame_with::<EquatorialTrueOfDate, _>(&jd, &model)
        .to_spherical();
    let moon_horizon = project_position_with(&moon, jd, &model);
    let mars_horizon = project_position_with(&mars, jd, &model);
    let phase = moon_phase_from_topocentric(&moon, &sun);
    let limb = moon_bright_limb_from_topocentric(&moon, &sun);
    println!("Moon: {:?}, Mars: {:?}", moon_horizon, mars_horizon);
    println!(
        "Illumination: {}, bright limb: {limb:?}",
        phase.illuminated_fraction
    );
    Ok(())
}
