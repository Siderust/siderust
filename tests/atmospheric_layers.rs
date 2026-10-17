// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon
//! Analytic column, physical-invariant and strict-domain regressions.
#![cfg(feature = "atmosphere")]

use siderust::atmosphere::{
    rayleigh_optical_depth_bodhaine99, AtmosphereProfile, AtmosphericOptics, LayeredAtmosphere,
    LayeredAtmosphereError as Error,
};
use siderust::optica::medium::{InverseLength, OpticalCoefficients};
use siderust::qtty::{
    unit::Meter, Albedos, Hectopascals, Kilometers, Meters, Nanometers, Quantity,
};

fn c(v: f64) -> Quantity<InverseLength<Meter>> {
    Quantity::new(v)
}
fn domain() -> [Nanometers; 2] {
    [Nanometers::new(360.), Nanometers::new(570.)]
}
fn edges(v: &[f64]) -> Vec<Meters> {
    v.iter().copied().map(Meters::new).collect()
}
fn optics(r: f64, a: f64, ssa: f64, m: f64) -> AtmosphericOptics {
    AtmosphericOptics::from_extinction(c(r), c(a), Albedos::new(ssa), c(m)).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-14 * b.abs().max(1.), "{a} != {b}");
}

#[test]
fn analytic_partial_and_two_layer_columns_and_boundaries() {
    let column = LayeredAtmosphere::try_new(
        &edges(&[100., 200., 400.]),
        domain(),
        vec![
            optics(0.01, 0.02, 0.25, 0.03),
            optics(0.04, 0.05, 0.8, 0.06),
        ],
    )
    .unwrap();
    assert_eq!(column.edges(), edges(&[100., 200., 400.]));
    assert_eq!(column.wavelength_domain(), domain());
    assert!(column.profile().is_none());
    for (lo, hi, expected) in [
        (100., 400., 36.),
        (150., 300., 18.),
        (125., 175., 3.),
        (100., 100., 0.),
        (400., 400., 0.),
    ] {
        let tau = column
            .vertical_optical_depth(Nanometers::new(500.), Meters::new(lo), Meters::new(hi))
            .unwrap();
        close(tau.extinction().value(), expected);
        close(
            tau.extinction().value(),
            tau.scattering().value() + tau.absorption().value(),
        );
    }
    let tau = column
        .vertical_optical_depth(Nanometers::new(500.), Meters::new(150.), Meters::new(300.))
        .unwrap();
    close(tau.rayleigh.value(), 4.5);
    close(tau.aerosol_scattering.value(), 4.25);
    close(tau.aerosol_absorption.value(), 1.75);
    close(tau.molecular_absorption.value(), 7.5);
    for h in [200., 400.] {
        close(
            column
                .coefficients(Nanometers::new(360.), Meters::new(h))
                .unwrap()
                .rayleigh()
                .sigma_s
                .value(),
            0.04,
        );
    }
    close(
        column
            .coefficients(Nanometers::new(570.), Meters::new(100.))
            .unwrap()
            .molecular()
            .sigma_a
            .value(),
        0.03,
    );
    assert!(format!("{column:?}").len() > 10);
}

#[test]
fn split_merge_invariance_and_physical_properties() {
    for r in [0., 1e-8, 0.01] {
        for a in [0., 1e-8, 0.02] {
            for m in [0., 1e-8, 0.03] {
                for ssa in [0., 0.2, 1.] {
                    let op = optics(r, a, ssa, m);
                    let total = op.total().unwrap();
                    close(total.sigma_t.value(), r + a + m);
                    close(total.sigma_s.value(), r + a * ssa);
                    close(total.sigma_a.value(), m + a * (1. - ssa));
                    if total.sigma_t.value() == 0. {
                        assert_eq!(total.ssa.value(), 0.);
                    }
                    close(op.aerosol().sigma_s.value(), a * ssa);
                    let one = LayeredAtmosphere::try_new(&edges(&[0., 1000.]), domain(), vec![op])
                        .unwrap();
                    let split = LayeredAtmosphere::try_new(
                        &edges(&[0., 50., 250., 1000.]),
                        domain(),
                        vec![op; 3],
                    )
                    .unwrap();
                    let query = |column: &LayeredAtmosphere| {
                        column
                            .vertical_optical_depth(
                                Nanometers::new(400.),
                                Meters::new(10.),
                                Meters::new(900.),
                            )
                            .unwrap()
                    };
                    close(query(&one).extinction().value(), 890. * (r + a + m));
                    close(
                        query(&one).extinction().value(),
                        query(&split).extinction().value(),
                    );
                    assert_eq!(query(&split), query(&split));
                }
            }
        }
    }
}

#[test]
fn coefficient_state_validation_including_mutation_and_overflow() {
    let zero = OpticalCoefficients::<Meter>::transparent();
    let good = OpticalCoefficients::<Meter>::try_new(0.5, 0.25).unwrap();
    for v in [-1., f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(AtmosphericOptics::from_extinction(c(v), c(0.), Albedos::new(0.5), c(0.)).is_err());
        assert!(AtmosphericOptics::from_extinction(c(0.), c(v), Albedos::new(0.), c(0.)).is_err());
        assert!(AtmosphericOptics::from_extinction(c(0.), c(0.), Albedos::new(0.5), c(v)).is_err());
    }
    for ssa in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            AtmosphericOptics::from_extinction(c(0.), c(0.), Albedos::new(ssa), c(0.)).unwrap_err(),
            Error::Albedo
        );
    }
    for field in 0..4 {
        for v in [-1., f64::NAN, f64::INFINITY, 0.1] {
            let mut bad = good;
            match field {
                0 => bad.sigma_a = c(v),
                1 => bad.sigma_s = c(v),
                2 => bad.sigma_t = c(v),
                _ => bad.ssa = Albedos::new(v),
            }
            assert_eq!(
                AtmosphericOptics::try_new(zero, bad, zero).unwrap_err(),
                Error::Coefficients
            );
        }
    }
    assert!(AtmosphericOptics::try_new(good, zero, zero).is_err());
    assert!(AtmosphericOptics::try_new(zero, zero, good).is_err());
    let overflow = OpticalCoefficients::<Meter>::try_new(f64::MAX, f64::MAX).unwrap();
    assert!(AtmosphericOptics::try_new(zero, overflow, zero).is_err());
    assert!(
        AtmosphericOptics::from_extinction(c(f64::MAX), c(f64::MAX), Albedos::new(1.), c(0.))
            .is_err()
    );
    assert!(
        AtmosphericOptics::from_extinction(c(0.), c(f64::MAX), Albedos::new(0.), c(f64::MAX))
            .is_err()
    );
    assert_eq!(
        Error::Coefficients.to_string(),
        "optical coefficients must be finite, nonnegative and internally consistent"
    );
}

#[test]
fn construction_domains_and_query_failures() {
    let op = optics(0.01, 0., 0., 0.);
    for e in [
        vec![],
        vec![0.],
        vec![0., 0.],
        vec![1., 0.],
        vec![0., f64::NAN],
        vec![f64::NEG_INFINITY, 0.],
        vec![-f64::MAX, f64::MAX],
        vec![0., 2., 1.],
    ] {
        assert!(LayeredAtmosphere::try_new(&edges(&e), domain(), vec![op]).is_err());
    }
    assert!(LayeredAtmosphere::try_new(&edges(&[0., 1., 2.]), domain(), vec![op]).is_err());
    for (lo, hi) in [
        (0., 570.),
        (-1., 570.),
        (500., 500.),
        (570., 360.),
        (f64::NAN, 570.),
        (360., f64::INFINITY),
    ] {
        assert!(LayeredAtmosphere::try_new(
            &edges(&[0., 1.]),
            [Nanometers::new(lo), Nanometers::new(hi)],
            vec![op]
        )
        .is_err());
    }
    let column = LayeredAtmosphere::try_new(&edges(&[100., 200.]), domain(), vec![op]).unwrap();
    for w in [359., 571., f64::NAN, f64::INFINITY] {
        assert_eq!(
            column
                .coefficients(Nanometers::new(w), Meters::new(100.))
                .unwrap_err(),
            Error::Wavelength
        );
        assert!(column
            .vertical_optical_depth(Nanometers::new(w), Meters::new(100.), Meters::new(200.))
            .is_err());
    }
    for h in [99., 201., f64::NAN, f64::INFINITY] {
        assert_eq!(
            column
                .coefficients(Nanometers::new(400.), Meters::new(h))
                .unwrap_err(),
            Error::Height
        );
        for (lo, hi) in [(h, 150.), (150., h)] {
            assert!(column
                .vertical_optical_depth(Nanometers::new(400.), Meters::new(lo), Meters::new(hi))
                .is_err());
        }
    }
    assert!(column
        .vertical_optical_depth(Nanometers::new(400.), Meters::new(175.), Meters::new(125.))
        .is_err());
}

fn from_profile(
    p: AtmosphereProfile,
    e: &[f64],
    rf: &[f64],
    af: &[f64],
) -> Result<LayeredAtmosphere, Error> {
    LayeredAtmosphere::from_profile(p, &edges(e), domain(), rf, af, Albedos::new(0.7))
}

#[test]
fn column_fraction_conservation_local_pressure_and_no_double_height() {
    let mut p = AtmosphereProfile::EL_PARANAL;
    p.surface_pressure = Hectopascals::new(1013.25 * (-2.635_f64 / 8.).exp());
    let col = from_profile(p, &[2635., 3000., 10000.], &[0.2, 0.8], &[0.9, 0.1]).unwrap();
    assert_eq!(col.profile(), Some(&p));
    for w in [360., 375., 500., 570.] {
        let wavelength = Nanometers::new(w);
        let tau = col
            .vertical_optical_depth(wavelength, Meters::new(2635.), Meters::new(10000.))
            .unwrap();
        let reference = rayleigh_optical_depth_bodhaine99(
            wavelength,
            Hectopascals::new(1013.25),
            Kilometers::new(2.635),
            Kilometers::new(8.),
        )
        .value();
        close(tau.rayleigh.value(), reference);
        close(
            tau.extinction().value(),
            p.optical_depth(wavelength).value(),
        );
        let local = col.coefficients(wavelength, Meters::new(2700.)).unwrap();
        close(local.rayleigh().sigma_s.value() * 365., reference * 0.2);
        close(
            tau.aerosol_scattering.value()
                / (tau.aerosol_scattering.value() + tau.aerosol_absorption.value()),
            0.7,
        );
    }
    let one = from_profile(p, &[2635., 10000.], &[1.], &[1.]).unwrap();
    close(
        one.vertical_optical_depth(
            Nanometers::new(500.),
            Meters::new(2635.),
            Meters::new(10000.),
        )
        .unwrap()
        .extinction()
        .value(),
        p.optical_depth(Nanometers::new(500.)).value(),
    );
}

#[test]
fn fractions_and_profile_parameters_are_checked() {
    let p = AtmosphereProfile::EL_PARANAL;
    let e = [2635., 4000., 10000.];
    for fractions in [
        vec![],
        vec![1.],
        vec![0.1, 0.1],
        vec![-0.1, 1.1],
        vec![f64::NAN, 1.],
        vec![f64::MAX, f64::MAX],
    ] {
        assert!(from_profile(p, &e, &fractions, &[0.5, 0.5]).is_err());
        assert!(from_profile(p, &e, &[0.5, 0.5], &fractions).is_err());
    }
    assert!(from_profile(p, &e, &[0., 1.], &[1., 0.]).is_ok());
    assert!(from_profile(p, &e, &[0.5, 0.5 + f64::EPSILON], &[1., 0.]).is_ok());
    assert!(from_profile(p, &[0., 10000.], &[1.], &[1.]).is_err());
    for field in 0..7 {
        for value in [f64::NAN, f64::INFINITY, -1.] {
            let mut bad = p;
            match field {
                0 => bad.surface_pressure = Hectopascals::new(value),
                1 => bad.rayleigh_scale_height = Kilometers::new(value),
                2 => bad.observer_altitude = Kilometers::new(value),
                3 => bad.mie_params.tau0 = siderust::qtty::OpticalDepths::new(value),
                4 => bad.mie_params.alpha = value,
                5 => bad.mie_params.lambda_ref = Nanometers::new(value),
                _ => bad.rayleigh_scale_height = Kilometers::new(0.),
            }
            // A finite negative aerosol exponent is physically allowed.
            if field == 4 && value == -1. {
                assert!(from_profile(bad, &e, &[0.5, 0.5], &[0.5, 0.5]).is_ok());
            } else {
                assert!(from_profile(bad, &e, &[0.5, 0.5], &[0.5, 0.5]).is_err());
            }
        }
    }
    let mut bad = p;
    bad.mie_params.lambda_ref = Nanometers::new(0.);
    assert!(from_profile(bad, &e, &[0.5, 0.5], &[0.5, 0.5]).is_err());
    for d in [
        [Nanometers::new(350.), Nanometers::new(570.)],
        [Nanometers::new(360.), Nanometers::new(600.)],
    ] {
        assert!(LayeredAtmosphere::from_profile(
            p,
            &edges(&e),
            d,
            &[0.5, 0.5],
            &[0.5, 0.5],
            Albedos::new(0.5)
        )
        .is_err());
    }
    assert!(LayeredAtmosphere::from_profile(
        p,
        &edges(&e),
        domain(),
        &[0.5, 0.5],
        &[0.5, 0.5],
        Albedos::new(2.)
    )
    .is_err());
    // No silent overflowing conversion of a fraction into inverse thickness.
    let mut sea = p;
    sea.observer_altitude = Kilometers::new(0.);
    assert!(from_profile(sea, &[0., f64::from_bits(1)], &[1.], &[1.]).is_err());
}

#[test]
fn finite_coefficients_do_not_allow_overflowing_column_integrals() {
    let op = optics(f64::MAX / 4., 0., 0., 0.);
    let column = LayeredAtmosphere::try_new(&edges(&[0., 16.]), domain(), vec![op]).unwrap();
    assert!(column
        .vertical_optical_depth(Nanometers::new(500.), Meters::new(0.), Meters::new(16.))
        .is_err());
    let op = optics(f64::MAX / 4., f64::MAX / 4., 0., 0.);
    let column = LayeredAtmosphere::try_new(&edges(&[0., 3.]), domain(), vec![op]).unwrap();
    assert!(column
        .vertical_optical_depth(Nanometers::new(500.), Meters::new(0.), Meters::new(3.))
        .is_err());
    let mut p = AtmosphereProfile::EL_PARANAL;
    p.mie_params.alpha = f64::MAX;
    assert!(from_profile(p, &[2635., 5000.], &[1.], &[1.]).is_err());
}

#[test]
fn fallible_optica_medium_binds_origin_and_datum_explicitly() {
    use siderust::coordinates::{
        cartesian::Position,
        centers::{Geodetic, Topocentric},
        frames::{Horizontal, ECEF},
    };
    use siderust::optica::medium::TryMedium;
    let column = LayeredAtmosphere::try_new(
        &edges(&[100., 200.]),
        domain(),
        vec![optics(0.01, 0., 0., 0.)],
    )
    .unwrap();
    let site = Geodetic::<ECEF>::new(
        siderust::qtty::Degrees::new(0.),
        siderust::qtty::Degrees::new(0.),
        Meters::new(25.),
    );
    // Explicit absolute datum: the 25-m ellipsoidal site height is NOT substituted
    // for the caller's 100-m column origin. Local z=50 means absolute z=150.
    let medium = column
        .at_origin::<Topocentric>(site, Meters::new(100.))
        .unwrap();
    let pos = |x, y, z| Position::<Topocentric, Horizontal, Meter>::new_with_params(site, x, y, z);
    close(
        medium
            .try_coefficients(pos(0., 0., 50.), Nanometers::new(500.))
            .unwrap()
            .sigma_s
            .value(),
        0.01,
    );
    assert!(column
        .at_origin::<Topocentric>(site, Meters::new(f64::NAN))
        .is_err());
    for p in [
        pos(f64::NAN, 0., 50.),
        pos(0., f64::INFINITY, 50.),
        pos(0., 0., 101.),
        pos(0., 0., f64::MAX),
    ] {
        assert!(medium.try_coefficients(p, Nanometers::new(500.)).is_err());
    }
    assert!(medium
        .try_coefficients(pos(0., 0., 50.), Nanometers::new(600.))
        .is_err());
    let other = Geodetic::<ECEF>::new(
        siderust::qtty::Degrees::new(1.),
        siderust::qtty::Degrees::new(0.),
        Meters::new(25.),
    );
    let p = Position::new_with_params(other, 0., 0., 50.);
    assert_eq!(
        medium
            .try_coefficients(p, Nanometers::new(500.))
            .unwrap_err(),
        Error::Position
    );
}
