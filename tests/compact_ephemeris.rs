// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compact ephemeris coverage and pointing regressions.

#![cfg(feature = "embedded-compact")]

use siderust::astro::{eop::NullEop, nutation::Iau2000B};
use siderust::coordinates::{
    cartesian::Position,
    centers::{Geocentric, Geodetic},
    frames::{EclipticMeanJ2000, ECEF},
    transform::{AstroContext, TransformFrame},
};
use siderust::ephemeris::{CompactEphemeris, Ephemeris, EphemerisError, MajorPlanet};
use siderust::event::{
    horizontal::{
        equatorial_to_horizontal_true_of_date_with_ctx,
        geocentric_j2000_to_apparent_topocentric_with_ctx,
    },
    lunar::{
        phase::moon_phase_from_topocentric, position_angle::moon_bright_limb_from_topocentric,
    },
};
use siderust::qtty::{AstronomicalUnit, Degrees, Meters};
use siderust::time::JulianDate;
use tempoch::TDB;

fn tt(tdb: f64) -> JulianDate {
    tempoch::JulianDate::<TDB>::new(tdb).to_scale()
}

#[test]
fn coverage_and_velocity_stencils() {
    let first = CompactEphemeris::FIRST_JD_TDB;
    let end = CompactEphemeris::END_JD_TDB;
    let next_below_end = f64::from_bits(end.to_bits() - 1);
    for tdb in [first, first + 1000.0, next_below_end] {
        let jd = tt(tdb);
        assert!(CompactEphemeris::try_sun_barycentric(jd).is_ok());
        assert!(CompactEphemeris::try_moon_geocentric(jd).is_ok());
        let v = CompactEphemeris::try_earth_barycentric_velocity(jd).unwrap();
        assert!((0.016..0.019).contains(&v.magnitude().value()));
    }

    for tdb in [f64::from_bits(first.to_bits() - 1), end] {
        let jd = tt(tdb);
        assert!(matches!(
            CompactEphemeris::try_moon_geocentric(jd),
            Err(EphemerisError::OutOfRange { .. })
        ));
        assert!(matches!(
            CompactEphemeris::try_earth_barycentric_velocity(jd),
            Err(EphemerisError::OutOfRange { .. })
        ));
    }
}

#[test]
fn offline_topocentric_geometry() {
    let jd = tt(2_461_355.5);
    let site = Geodetic::<ECEF>::new(
        Degrees::new(-0.1275),
        Degrees::new(51.507222),
        Meters::new(1.0),
    );
    let ctx = AstroContext::<CompactEphemeris, NullEop>::default();
    let model = ctx.with_model::<Iau2000B>();
    let moon = CompactEphemeris::try_moon_geocentric(jd).unwrap();
    let earth = CompactEphemeris::try_earth_heliocentric(jd).unwrap();
    let sun = Position::<Geocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
        -earth.x(),
        -earth.y(),
        -earth.z(),
    );
    let moon_topo =
        geocentric_j2000_to_apparent_topocentric_with_ctx(&moon.to_frame(), site, jd, &model);
    let sun_topo =
        geocentric_j2000_to_apparent_topocentric_with_ctx(&sun.to_frame(), site, jd, &model);
    let horizontal = equatorial_to_horizontal_true_of_date_with_ctx(&moon_topo, site, jd, &model);
    assert!(horizontal.az().is_finite() && horizontal.alt().is_finite());
    let phase = moon_phase_from_topocentric(&moon_topo, &sun_topo);
    assert!((0.0..=1.0).contains(&phase.illuminated_fraction.value()));
    assert!(moon_bright_limb_from_topocentric(&moon_topo, &sun_topo)
        .unwrap()
        .is_finite());
    for planet in [
        MajorPlanet::Mercury,
        MajorPlanet::Venus,
        MajorPlanet::Mars,
        MajorPlanet::Jupiter,
        MajorPlanet::Saturn,
        MajorPlanet::Uranus,
        MajorPlanet::Neptune,
    ] {
        let position = CompactEphemeris::try_planet_barycentric(planet, jd).unwrap();
        assert!(position
            .as_array()
            .iter()
            .all(|component| component.is_finite()));
    }
}

#[derive(Default)]
struct ReferenceEop;

impl siderust::astro::eop::EopProvider for ReferenceEop {
    fn eop_at(&self, _: tempoch::JulianDate<tempoch::UTC>) -> siderust::astro::eop::EopValues {
        use siderust::qtty::{Arcseconds, Seconds};
        siderust::astro::eop::EopValues {
            dut1: Seconds::new(0.4),
            xp: Arcseconds::new(0.15),
            yp: Arcseconds::new(-0.25),
            ..Default::default()
        }
    }
}

fn angular_error_arcsec(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    norm(cross).atan2(dot).to_degrees() * 3600.0
}

fn norm(vector: [f64; 3]) -> f64 {
    vector
        .into_iter()
        .map(|axis| axis * axis)
        .sum::<f64>()
        .sqrt()
}

fn horizon_vector(az: f64, alt: f64) -> [f64; 3] {
    let (sa, ca) = az.to_radians().sin_cos();
    let (se, ce) = alt.to_radians().sin_cos();
    [ca * ce, sa * ce, se]
}

#[test]
fn de440s_erfa_pointing_and_velocity() {
    use siderust::astro::aberration::apply_aberration_to_direction_with_velocity;
    use siderust::coordinates::frames::EquatorialMeanJ2000;
    use siderust::qtty::{Kilometer, Kilometers};

    let context = AstroContext::<CompactEphemeris, ReferenceEop>::default();
    let model = context.with_model::<Iau2000B>();
    let null = AstroContext::<CompactEphemeris, NullEop>::default();
    let null_model = null.with_model::<Iau2000B>();
    let mut maxima = [0.0_f64; 12];
    let mut nearest_horizon = 90.0_f64;
    for row in include_str!("fixtures/compact_reference.csv")
        .lines()
        .skip(1)
    {
        let data: Vec<f64> = row.split(',').map(|field| field.parse().unwrap()).collect();
        let jd = tt(data[0]);
        let tdb = jd.to_scale::<TDB>();
        let (tt_hi, tt_lo) = jd.raw_seconds_pair();
        let (tdb_hi, tdb_lo) = tdb.raw_seconds_pair();
        let time_error = ((tdb_hi - tt_hi) + (tdb_lo - tt_lo)).value() - data[18];
        let nut = siderust::astro::nutation::nutation_iau2000b(jd);
        let site = Geodetic::new(
            Degrees::new(data[1]),
            Degrees::new(data[2]),
            Meters::new(data[3]),
        );
        let earth = CompactEphemeris::try_earth_barycentric(jd).unwrap();
        let ecliptic = if data[4] == 301.0 {
            CompactEphemeris::try_moon_geocentric(jd).unwrap()
        } else {
            let body = match data[4] as u8 {
                10 => CompactEphemeris::try_sun_barycentric(jd).unwrap(),
                id => CompactEphemeris::try_planet_barycentric(
                    match id {
                        1 => MajorPlanet::Mercury,
                        2 => MajorPlanet::Venus,
                        4 => MajorPlanet::Mars,
                        5 => MajorPlanet::Jupiter,
                        6 => MajorPlanet::Saturn,
                        7 => MajorPlanet::Uranus,
                        8 => MajorPlanet::Neptune,
                        _ => unreachable!(),
                    },
                    jd,
                )
                .unwrap(),
            };
            let [x, y, z] = core::array::from_fn(|i| {
                (body.as_array()[i] - earth.as_array()[i]).to::<Kilometer>()
            });
            Position::new(x, y, z)
        };
        let eq: Position<Geocentric, EquatorialMeanJ2000, Kilometer> = ecliptic.to_frame();
        let actual = ecliptic.as_array().map(|axis| axis.value());
        let reference = [data[5], data[6], data[7]];
        let velocity = CompactEphemeris::try_earth_barycentric_velocity(jd).unwrap();
        let error_velocity = norm(core::array::from_fn(|i| {
            velocity.as_array()[i].value() - data[10 + i]
        }));
        let topo = geocentric_j2000_to_apparent_topocentric_with_ctx(&eq, site, jd, &model);
        let horizon = equatorial_to_horizontal_true_of_date_with_ctx(&topo, site, jd, &model);
        let null_topo =
            geocentric_j2000_to_apparent_topocentric_with_ctx(&eq, site, jd, &null_model);
        let null_horizon =
            equatorial_to_horizontal_true_of_date_with_ctx(&null_topo, site, jd, &null_model);
        let observed = horizon_vector(horizon.az().value(), horizon.alt().value());
        let no_eop = horizon_vector(null_horizon.az().value(), null_horizon.alt().value());
        let ref_horizon = horizon_vector(data[8], data[9]);
        let az_error = (horizon.az().value() - data[8] + 180.0).rem_euclid(360.0) - 180.0;

        let corrected = apply_aberration_to_direction_with_velocity(
            eq.to_spherical().direction().to_cartesian(),
            &velocity.to_frame(),
        );
        let direction = corrected.as_array();
        let corrected = Position::<Geocentric, EquatorialMeanJ2000, Kilometer>::new(
            Kilometers::new(direction[0] * norm(actual)),
            Kilometers::new(direction[1] * norm(actual)),
            Kilometers::new(direction[2] * norm(actual)),
        );
        let topo = geocentric_j2000_to_apparent_topocentric_with_ctx(&corrected, site, jd, &model);
        let horizon = equatorial_to_horizontal_true_of_date_with_ctx(&topo, site, jd, &model);
        let apparent_error = angular_error_arcsec(
            horizon_vector(horizon.az().value(), horizon.alt().value()),
            horizon_vector(data[16], data[17]),
        );
        let radial_error = if data[4] == 301.0 {
            (norm(actual) - norm(reference)).abs()
        } else {
            0.0
        };
        nearest_horizon = nearest_horizon.min(data[9].abs());
        for (maximum, value) in maxima.iter_mut().zip([
            angular_error_arcsec(actual, reference),
            angular_error_arcsec(observed, ref_horizon),
            az_error.abs() * 3600.0,
            (data[9]
                - equatorial_to_horizontal_true_of_date_with_ctx(
                    &geocentric_j2000_to_apparent_topocentric_with_ctx(&eq, site, jd, &model),
                    site,
                    jd,
                    &model,
                )
                .alt()
                .value())
            .abs()
                * 3600.0,
            radial_error,
            error_velocity * 149_597_870_700.0 / 86400.0,
            error_velocity / 173.144632674240 * 206_264.806,
            apparent_error,
            angular_error_arcsec(observed, no_eop),
            time_error.abs(),
            (nut.dpsi.value() - data[19]).abs().to_degrees() * 3600.0,
            (nut.deps.value() - data[20]).abs().to_degrees() * 3600.0,
        ]) {
            *maximum = maximum.max(value);
        }
    }
    eprintln!("maxima (geo, topo, az, alt arcsec, lunar km, m/s, velocity aberration arcsec, apparent arcsec, EOP arcsec, TT-TDB seconds, nutation longitude and obliquity arcsec): {maxima:?}; nearest horizon {nearest_horizon} deg");
    assert!(maxima[0] < 60.0 && maxima[1] < 60.0 && maxima[7] < 60.0);
    assert!(maxima[4] < 150.0 && maxima[5] < 10.0);
}
