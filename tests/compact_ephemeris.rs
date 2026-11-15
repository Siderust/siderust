// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Compact ephemeris coverage and pointing regressions.

#![cfg(feature = "embedded-compact")]

use siderust::astro::{
    aberration::apply_aberration_to_direction_with_velocity,
    eop::{EopProvider, EopValues},
    nutation::{nutation_iau2000b, Iau2000B},
};
use siderust::bodies::{Jupiter, Mars, Mercury, Neptune, Saturn, Uranus, Venus};
use siderust::coordinates::{
    cartesian::Position,
    centers::{Geocentric, Geodetic, Heliocentric, Topocentric},
    frames::{EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialTrueOfDate, ECEF},
    transform::{horizontal::project_position_with, AstroContext, TransformCenter, TransformFrame},
};
use siderust::ephemeris::{
    CompactEphemeris, Ephemeris, EphemerisError, MajorPlanet, PlanetEphemerisError, PlanetState,
};
use siderust::event::lunar::{
    phase::{moon_phase_from_topocentric, moon_phase_topocentric},
    position_angle::{moon_bright_limb_from_topocentric, moon_bright_limb_position_angle},
};
use siderust::qtty::{Arcseconds, AstronomicalUnit, Degrees, Kilometer, Meters, Seconds};
use siderust::time::JulianDate;
use std::collections::BTreeMap;
use tempoch::{TDB, UTC};

type EarthRelative<U> = Position<Geocentric, EclipticMeanJ2000, U>;
type HeliocentricPosition = Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;
type HeliocentricQuery = fn(JulianDate) -> Result<HeliocentricPosition, EphemerisError>;

const PLANETS: [(f64, MajorPlanet, HeliocentricQuery); 7] = [
    (
        1.0,
        MajorPlanet::Mercury,
        <CompactEphemeris as PlanetState<Mercury>>::try_heliocentric,
    ),
    (
        2.0,
        MajorPlanet::Venus,
        <CompactEphemeris as PlanetState<Venus>>::try_heliocentric,
    ),
    (
        4.0,
        MajorPlanet::Mars,
        <CompactEphemeris as PlanetState<Mars>>::try_heliocentric,
    ),
    (
        5.0,
        MajorPlanet::Jupiter,
        <CompactEphemeris as PlanetState<Jupiter>>::try_heliocentric,
    ),
    (
        6.0,
        MajorPlanet::Saturn,
        <CompactEphemeris as PlanetState<Saturn>>::try_heliocentric,
    ),
    (
        7.0,
        MajorPlanet::Uranus,
        <CompactEphemeris as PlanetState<Uranus>>::try_heliocentric,
    ),
    (
        8.0,
        MajorPlanet::Neptune,
        <CompactEphemeris as PlanetState<Neptune>>::try_heliocentric,
    ),
];
const AU_PER_DAY_IN_M_PER_S: f64 = 149_597_870_700.0 / 86_400.0;

fn tt(tdb: f64) -> JulianDate {
    tempoch::JulianDate::<TDB>::new(tdb).to_scale()
}

fn previous_jd(jd: f64) -> f64 {
    f64::from_bits(jd.to_bits() - 1)
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn angular_error_arcsec(a: [f64; 3], b: [f64; 3]) -> f64 {
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    norm(cross).atan2(dot(a, b)).to_degrees() * 3600.0
}

fn assert_below(case: &str, metric: &str, error: f64, tolerance: f64) {
    assert!(
        error.abs() < tolerance,
        "{case}: {metric} discrepancy {error:.6e} exceeds tolerance {tolerance:.6e}"
    );
}

#[test]
fn coverage_and_velocity_stencils() {
    let first = CompactEphemeris::FIRST_JD_TDB;
    let end = CompactEphemeris::END_JD_TDB;
    // The documented 675 s stencil changes sides at these two epochs.
    let stencil_step = 1.0 / 128.0;
    let transitions = [first + stencil_step, end - stencil_step];
    for tdb in [
        first - 1.0,
        previous_jd(first),
        first,
        first + 1_000.0,
        (first + end) / 2.0,
        previous_jd(end),
        end,
        end + 1.0,
    ] {
        let jd = tt(tdb);
        let check = |body: &str, result: Result<f64, EphemerisError>| {
            if (first..end).contains(&tdb) {
                assert!(
                    matches!(result, Ok(magnitude) if magnitude.is_finite() && magnitude > 0.0),
                    "{body} at JD(TDB) {tdb}: expected finite positive magnitude, got {result:?}"
                );
            } else {
                let expected = EphemerisError::OutOfRange {
                    jd: tdb,
                    start_jd: first,
                    end_jd: end,
                };
                assert_eq!(result, Err(expected), "{body} at JD(TDB) {tdb}");
            }
        };
        check(
            "Sun barycentric",
            CompactEphemeris::try_sun_barycentric(jd).map(|p| p.distance().value()),
        );
        check(
            "Earth barycentric",
            CompactEphemeris::try_earth_barycentric(jd).map(|p| p.distance().value()),
        );
        check(
            "Earth heliocentric",
            CompactEphemeris::try_earth_heliocentric(jd).map(|p| p.distance().value()),
        );
        check(
            "Moon geocentric",
            CompactEphemeris::try_moon_geocentric(jd).map(|p| p.distance().value()),
        );
        check(
            "Earth barycentric speed",
            CompactEphemeris::try_earth_barycentric_velocity(jd).map(|v| v.magnitude().value()),
        );
        for (_, planet, heliocentric) in PLANETS {
            check(
                &format!("{planet:?} heliocentric"),
                heliocentric(jd).map(|p| p.distance().value()),
            );
            let barycentric =
                CompactEphemeris::try_planet_barycentric(planet, jd).map_err(|err| match err {
                    PlanetEphemerisError::Ephemeris(err) => err,
                    _ => panic!("{planet:?} barycentric at JD(TDB) {tdb}: unexpected {err:?}"),
                });
            check(
                &format!("{planet:?} barycentric"),
                barycentric.map(|p| p.distance().value()),
            );
        }
    }

    let earth = |tdb| {
        CompactEphemeris::try_earth_barycentric(tt(tdb))
            .unwrap_or_else(|err| panic!("Earth position at JD(TDB) {tdb}: {err}"))
            .as_array()
            .map(|axis| axis.value())
    };
    let velocity = |tdb| {
        CompactEphemeris::try_earth_barycentric_velocity(tt(tdb))
            .unwrap_or_else(|err| panic!("Earth velocity at JD(TDB) {tdb}: {err}"))
            .as_array()
            .map(|axis| axis.value())
    };
    // A fourth-order one-sided derivative taken into the interval uses a
    // different, exactly representable sample spacing (1/32 day).
    // Orbital jerk is of order 10 m/s/day² (v * (2π/year)², plus lunar
    // reflex motion), putting the 675 s second-order error below 1 mm/s.
    // Allow 5 mm/s for higher-frequency terms and roundoff.
    let tolerance_m_per_s = 0.005;
    for tdb in [
        first,
        previous_jd(transitions[0]),
        transitions[0],
        first + 1_000.0,
        (first + end) / 2.0,
        previous_jd(transitions[1]),
        transitions[1],
        previous_jd(end),
    ] {
        let step = if tdb < (first + end) / 2.0 {
            1.0 / 32.0
        } else {
            -1.0 / 32.0
        };
        let origin = earth(tdb);
        let [a, b, c, d] = [1.0, 2.0, 3.0, 4.0].map(|k| difference(earth(tdb + k * step), origin));
        let reference = core::array::from_fn(|i| {
            (48.0 * a[i] - 36.0 * b[i] + 16.0 * c[i] - 3.0 * d[i]) / (12.0 * step)
        });
        assert_below(
            &format!("Earth at JD(TDB) {tdb}"),
            "velocity vs position derivative (m/s)",
            norm(difference(velocity(tdb), reference)) * AU_PER_DAY_IN_M_PER_S,
            tolerance_m_per_s,
        );
    }
    for tdb in transitions {
        assert_below(
            &format!("Earth at JD(TDB) {tdb}"),
            "velocity jump at stencil transition (m/s)",
            norm(difference(velocity(tdb), velocity(previous_jd(tdb)))) * AU_PER_DAY_IN_M_PER_S,
            tolerance_m_per_s,
        );
    }
}

#[test]
fn offline_topocentric_geometry() {
    let jd = tt(2_461_355.5);
    let sites = [
        Geodetic::<ECEF>::new(
            Degrees::new(-0.1275),
            Degrees::new(51.507222),
            Meters::new(0.0),
        ),
        Geodetic::<ECEF>::new(
            Degrees::new(179.8725),
            Degrees::new(-51.507222),
            Meters::new(0.0),
        ),
    ];
    // Match the convenience functions' default context without loading an
    // EOP table. The default provider supplies its built-in UT1 fallback.
    let ctx = AstroContext::<CompactEphemeris>::default();
    let moon = CompactEphemeris::try_moon_geocentric(jd).unwrap();
    let moon_eq: Position<Geocentric, EquatorialMeanJ2000, Kilometer> = moon.to_frame();
    let earth = CompactEphemeris::try_earth_heliocentric(jd).unwrap();
    let sun = EarthRelative::<AstronomicalUnit>::new(-earth.x(), -earth.y(), -earth.z());
    let sun_eq: Position<Geocentric, EquatorialMeanJ2000, _> = sun.to_frame();

    let [(moon_a, sun_a), (moon_b, sun_b)] = sites.map(|site| {
        let case = format!("lunar geometry at {site:?}, JD(TT) {}", jd.value());
        let check = |metric, error, tolerance| assert_below(&case, metric, error, tolerance);
        let moon_relative: Position<Topocentric, _, _> = moon_eq.to_center_with(site, jd, &ctx);
        let moon_topo = siderust::coordinates::transform::PositionAstroExt::to_frame_with::<
            EquatorialTrueOfDate,
            _,
        >(&moon_relative, &jd, &ctx)
        .to_spherical();
        let sun_relative: Position<Topocentric, _, _> = sun_eq.to_center_with(site, jd, &ctx);
        let sun_topo = siderust::coordinates::transform::PositionAstroExt::to_frame_with::<
            EquatorialTrueOfDate,
            _,
        >(&sun_relative, &jd, &ctx)
        .to_spherical();
        let moon_cart = moon_topo.to_cartesian();
        let moon_km = moon_cart.as_array().map(|axis| axis.value());
        let sun_km = sun_topo
            .to_cartesian()
            .to_unit::<Kilometer>()
            .as_array()
            .map(|axis| axis.value());

        let topo_ecl: Position<Topocentric, EclipticMeanJ2000, Kilometer> =
            moon.to_center_with(site, jd, &ctx);
        let topo_eq: Position<Topocentric, EquatorialMeanJ2000, Kilometer> =
            moon_eq.to_center_with(site, jd, &ctx);
        let rotated: Position<Topocentric, EquatorialMeanJ2000, Kilometer> = topo_ecl.to_frame();
        let recovered: EarthRelative<Kilometer> = topo_ecl.to_center_with((), jd, &ctx);
        let relative_au: Position<Topocentric, _, _> = moon_eq
            .to_unit::<AstronomicalUnit>()
            .to_center_with(site, jd, &ctx);
        let moon_in_au = siderust::coordinates::transform::PositionAstroExt::to_frame_with::<
            EquatorialTrueOfDate,
            _,
        >(&relative_au, &jd, &ctx)
        .to_spherical();
        // Equivalence to within 10 µm, allowing f64 roundoff through the rotations.
        for (metric, error) in [
            (
                "frame/center order (km)",
                rotated.try_distance_to(&topo_eq).unwrap(),
            ),
            (
                "recovered geocentric Moon (km)",
                recovered.distance_to(&moon),
            ),
            (
                "Moon AU/km conversion (km)",
                moon_in_au
                    .to_cartesian()
                    .to_unit::<Kilometer>()
                    .try_distance_to(&moon_cart)
                    .unwrap(),
            ),
        ] {
            check(metric, error.value(), 1e-8);
        }

        // The illumination angle is at the Moon, between Moon→observer and
        // Moon→Sun. Build the triangle in a common distance unit.
        let phase = moon_phase_from_topocentric(&moon_topo, &sun_topo);
        let moon_to_sun = difference(sun_km, moon_km);
        let cosine = dot(moon_km.map(|x| -x), moon_to_sun) / (norm(moon_km) * norm(moon_to_sun));
        check(
            "phase angle (rad)",
            phase.phase_angle.value() - cosine.acos(),
            1e-12,
        );
        check(
            "illuminated fraction",
            phase.illuminated_fraction.value() - (1.0 + cosine) / 2.0,
            1e-12,
        );

        let from_epoch = moon_phase_topocentric::<CompactEphemeris>(jd, site);
        check(
            "precomputed vs direct elongation (rad)",
            phase.elongation.value() - from_epoch.elongation.value(),
            1e-12,
        );
        // The direct path uses geocentric distances for the phase triangle.
        // A one-Earth-radius distance change divided by the Sun distance
        // is about 4.4e-5 rad.
        check(
            "precomputed vs direct phase angle (rad)",
            phase.phase_angle.value() - from_epoch.phase_angle.value(),
            5e-5,
        );

        // Project the Sun onto the Moon's local east and north tangent axes.
        let (sin_ra, cos_ra) = moon_topo.ra().value().to_radians().sin_cos();
        let (sin_dec, cos_dec) = moon_topo.dec().value().to_radians().sin_cos();
        let east = [-sin_ra, cos_ra, 0.0];
        let north = [-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec];
        let sun_east = dot(sun_km, east);
        let waxing = sun_east < 0.0;
        let separation = (angular_error_arcsec(moon_km, sun_km) / 3600.0).to_radians();
        let signed_elongation = if waxing {
            separation
        } else {
            std::f64::consts::TAU - separation
        };
        check(
            "signed elongation (rad)",
            phase.elongation.value() - signed_elongation,
            1e-12,
        );
        assert_eq!(phase.waxing, waxing, "{case}: waxing flag");
        let expected_limb = sun_east.atan2(dot(sun_km, north));
        let limb = moon_bright_limb_from_topocentric(&moon_topo, &sun_topo)
            .unwrap()
            .value();
        let direct_limb =
            moon_bright_limb_position_angle::<CompactEphemeris, Topocentric>((jd, site))
                .unwrap()
                .value();
        for (metric, expected) in [
            ("north-through-east bright limb (rad)", expected_limb),
            ("precomputed vs direct bright limb (rad)", direct_limb),
        ] {
            let delta = limb - expected;
            check(metric, delta.sin().atan2(delta.cos()), 1e-12);
        }
        (moon_km, sun_km)
    });

    // Antipodal observers shift the common vectors by twice the same
    // WGS84 site radius. Direction changes scale inversely with distance.
    let site_xyz: Position<Geocentric, ECEF, Kilometer> = sites[0].to_cartesian();
    for (body, a, b) in [("Moon", moon_a, moon_b), ("Sun", sun_a, sun_b)] {
        let baseline = norm(difference(a, b));
        assert_below(
            body,
            "antipodal baseline (km)",
            baseline - 2.0 * site_xyz.distance().value(),
            1e-6,
        );
    }
    let lunar_parallax = angular_error_arcsec(moon_a, moon_b);
    let solar_parallax = angular_error_arcsec(sun_a, sun_b);
    assert!(
        (900.0..7600.0).contains(&lunar_parallax),
        "antipodal lunar parallax: {lunar_parallax} arcsec"
    );
    assert!(
        (1.0..19.0).contains(&solar_parallax),
        "antipodal solar parallax: {solar_parallax} arcsec"
    );
}

#[derive(Default)]
struct ReferenceEop;

impl EopProvider for ReferenceEop {
    fn eop_at(&self, _: tempoch::JulianDate<UTC>) -> EopValues {
        EopValues {
            dut1: Seconds::new(0.4),
            xp: Arcseconds::new(0.15),
            yp: Arcseconds::new(-0.25),
            ..Default::default()
        }
    }
}

fn horizon_vector(az: f64, alt: f64) -> [f64; 3] {
    let (sa, ca) = az.to_radians().sin_cos();
    let (se, ce) = alt.to_radians().sin_cos();
    [ca * ce, sa * ce, se]
}

#[test]
fn de440s_erfa_pointing_and_velocity() {
    let context = AstroContext::<CompactEphemeris, ReferenceEop>::default();
    let model = context.with_model::<Iau2000B>();
    let first = CompactEphemeris::FIRST_JD_TDB;
    let last = previous_jd(CompactEphemeris::END_JD_TDB);
    let mut coverage = BTreeMap::new();
    let mut just_below_horizon = false;
    let mut just_above_horizon = false;
    for row in include_str!("fixtures/compact_reference.csv")
        .lines()
        .skip(1)
    {
        let fields: Vec<f64> = row.split(',').map(|field| field.parse().unwrap()).collect();
        let &[tdb_jd, lon_deg, lat_deg, height_m, body_id, geo_x_km, geo_y_km, geo_z_km, az_deg, alt_deg, vx_au_day, vy_au_day, vz_au_day, dut1_s, xp_arcsec, yp_arcsec, aberrated_az_deg, aberrated_alt_deg, tdb_minus_tt_s, dpsi_rad, deps_rad] =
            fields.as_slice()
        else {
            panic!("Invalid DE440s/ERFA fixture row: {row}");
        };

        let jd = tt(tdb_jd);
        let site = Geodetic::new(
            Degrees::new(lon_deg),
            Degrees::new(lat_deg),
            Meters::new(height_m),
        );
        let (body, ecliptic) = if body_id == 301.0 {
            (
                "Moon".to_owned(),
                CompactEphemeris::try_moon_geocentric(jd).unwrap(),
            )
        } else {
            let (name, barycentric) = if body_id == 10.0 {
                (
                    "Sun".to_owned(),
                    CompactEphemeris::try_sun_barycentric(jd).unwrap(),
                )
            } else {
                let (_, planet, _) = PLANETS
                    .into_iter()
                    .find(|(id, _, _)| *id == body_id)
                    .unwrap_or_else(|| {
                        panic!("Unknown fixture body {body_id} at JD(TDB) {tdb_jd}")
                    });
                (
                    format!("{planet:?}"),
                    CompactEphemeris::try_planet_barycentric(planet, jd).unwrap(),
                )
            };
            let earth = CompactEphemeris::try_earth_barycentric(jd).unwrap();
            let [x, y, z] = core::array::from_fn(|i| {
                (barycentric.as_array()[i] - earth.as_array()[i]).to::<Kilometer>()
            });
            (name, Position::new(x, y, z))
        };
        let case =
            format!("{body} at JD(TDB) {tdb_jd}, site ({lon_deg}°, {lat_deg}°, {height_m} m)");
        let check = |metric, error, tolerance| assert_below(&case, metric, error, tolerance);
        let (at_start, in_middle, at_end) = coverage.entry(body).or_insert((false, false, false));
        *at_start |= tdb_jd == first;
        *in_middle |= (first + 3650.0..last - 3650.0).contains(&tdb_jd);
        *at_end |= tdb_jd == last;
        just_below_horizon |= (-1.0..0.0).contains(&alt_deg);
        just_above_horizon |= alt_deg > 0.0 && alt_deg < 1.0;

        let eop = context.eop_at_tt(jd);
        assert_eq!(
            [eop.dut1.value(), eop.xp.value(), eop.yp.value()],
            [dut1_s, xp_arcsec, yp_arcsec],
            "{case}: fixture EOP"
        );
        let actual = ecliptic.as_array().map(|axis| axis.value());
        let reference = [geo_x_km, geo_y_km, geo_z_km];
        let eq = ecliptic.to_frame();
        let horizontal = |position: &Position<Geocentric, EquatorialMeanJ2000, Kilometer>| {
            let topo: Position<Topocentric, _, _> = position.to_center_with(site, jd, &model);
            let of_date = siderust::coordinates::transform::PositionAstroExt::to_frame_with::<
                EquatorialTrueOfDate,
                _,
            >(&topo, &jd, &model);
            let angles = project_position_with(&of_date.to_spherical(), jd, &model);
            horizon_vector(angles.az().value(), angles.alt().value())
        };
        let observed = horizontal(&eq);
        let ref_horizon = horizon_vector(az_deg, alt_deg);
        let velocity = CompactEphemeris::try_earth_barycentric_velocity(jd).unwrap();
        let direction = apply_aberration_to_direction_with_velocity(
            eq.to_spherical().direction().to_cartesian(),
            &velocity.to_frame(),
        );
        let range = ecliptic.distance();
        assert!(
            range.is_finite() && range.value() > 0.0,
            "{case}: invalid geocentric range {range}"
        );
        let [x, y, z] = direction.as_array().map(|axis| range * axis);
        let aberrated = horizontal(&Position::new(x, y, z));
        let ref_aberrated = horizon_vector(aberrated_az_deg, aberrated_alt_deg);

        // The one-arcminute pointing requirement applies to each body and
        // sight line. Great-circle error stays meaningful at the azimuth
        // wrap and zenith and also bounds absolute altitude error.
        for (metric, observed, expected) in [
            ("geocentric pointing (arcsec)", actual, reference),
            (
                "geometric topocentric pointing (arcsec)",
                observed,
                ref_horizon,
            ),
            (
                "aberrated topocentric pointing (arcsec)",
                aberrated,
                ref_aberrated,
            ),
        ] {
            check(metric, angular_error_arcsec(observed, expected), 60.0);
        }
        if body_id == 301.0 {
            // At lunar perigee, a 150 km range error changes surface
            // parallax by less than 2 arcsec (Earth radius * Δrange / range²).
            check("lunar range (km)", range.value() - norm(reference), 150.0);
        }
        // A 10 m/s vector error contributes at most about 7 mas (Δv/c).
        // Retain this allowance for the analytic series and derivative.
        let velocity_error = difference(
            velocity.as_array().map(|axis| axis.value()),
            [vx_au_day, vy_au_day, vz_au_day],
        );
        check(
            "Earth barycentric velocity (m/s)",
            norm(velocity_error) * AU_PER_DAY_IN_M_PER_S,
            10.0,
        );

        // Check the applied correction itself: omitting it can still meet
        // the 60 arcsec pointing requirement. The 0.1 arcsec allowance
        // includes the velocity budget, small frame/direction differences,
        // and ERFA's additional gravitational term in the aberration model.
        let shift_error = difference(
            difference(aberrated, observed),
            difference(ref_aberrated, ref_horizon),
        );
        check(
            "annual aberration shift (arcsec)",
            norm(shift_error).to_degrees() * 3600.0,
            0.1,
        );

        let tdb = jd.to_scale::<TDB>();
        let (tt_hi, tt_lo) = jd.raw_seconds_pair();
        let (tdb_hi, tdb_lo) = tdb.raw_seconds_pair();
        let time_error = ((tdb_hi - tt_hi) + (tdb_lo - tt_lo)).value() - tdb_minus_tt_s;
        // Tempoch's seven-term TDB-TT approximation is accurate to about
        // 10 µs. Allow 20 µs against ERFA's full geocentric series.
        check("TDB minus TT (s)", time_error, 20e-6);

        let nutation = nutation_iau2000b(jd);
        // IAU 2000B is a milliarcsecond-class approximation. Allow 10 mas
        // per angle against IAU 2006/2000A, including the precession
        // adjustment and the difference between longitude and sky angle.
        for (metric, error) in [
            (
                "nutation in longitude (arcsec)",
                nutation.dpsi.value() - dpsi_rad,
            ),
            (
                "nutation in obliquity (arcsec)",
                nutation.deps.value() - deps_rad,
            ),
        ] {
            check(metric, error.to_degrees() * 3600.0, 0.010);
        }
    }
    assert_eq!(
        coverage.len(),
        PLANETS.len() + 2,
        "fixture must cover each planet, Moon and Sun"
    );
    for (body, epochs) in coverage {
        assert_eq!(
            epochs,
            (true, true, true),
            "{body}: fixture must sample the start, interior and last supported JD(TDB)"
        );
    }
    assert!(
        just_below_horizon && just_above_horizon,
        "fixture must sample within 1° on both sides of the geometric horizon"
    );
}
