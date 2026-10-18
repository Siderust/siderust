//! Physical and numerical regression tests for strict normalized phase tables.

#![cfg(feature = "atmosphere")]

use core::f64::consts::{PI, TAU};
use siderust::atmosphere::{
    NormalizedPhaseError as Error, NormalizedPhaseTable, TabulatedPhaseFunction,
};
use siderust::optica::phase::{
    DoubleHenyeyGreensteinPhaseFunction, HenyeyGreensteinPhaseFunction as Hg, PhaseFunction,
    RayleighPhaseFunction, ScatteringFactor,
};
use siderust::qtty::{Nanometers, Quantity, Radians};

fn table(xs: &[f64], ys: &[f64], data: &[f64]) -> Result<NormalizedPhaseTable, Error> {
    NormalizedPhaseTable::new(
        &xs.iter().map(|&x| Nanometers::new(x)).collect::<Vec<_>>(),
        &ys.iter().map(|&y| Radians::new(y)).collect::<Vec<_>>(),
        &data
            .iter()
            .map(|&p| Quantity::<ScatteringFactor>::new(p))
            .collect::<Vec<_>>(),
    )
}

fn at(t: &NormalizedPhaseTable, x: f64, y: f64) -> Result<f64, Error> {
    t.try_phase(Nanometers::new(x), Radians::new(y))
        .map(|p| p.value())
}

// Independent composite Simpson quadrature; never uses production integrals.
fn integrate(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let step = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(a + i as f64 * step);
    }
    sum * step / 3.0
}

fn linear_mass(ys: &[f64], ps: &[f64]) -> f64 {
    TAU * ys
        .windows(2)
        .zip(ps.windows(2))
        .map(|(y, p)| {
            integrate(
                |theta| {
                    let fraction = (theta - y[0]) / (y[1] - y[0]);
                    (p[0] * (1.0 - fraction) + p[1] * fraction) * theta.sin()
                },
                y[0],
                y[1],
                2000,
            )
        })
        .sum::<f64>()
}

#[test]
fn analytic_models_have_unity_mass_and_hg_moment() {
    let wavelength = Nanometers::new(550.0);
    let rayleigh = |theta| {
        RayleighPhaseFunction
            .phase(wavelength, Radians::new(theta))
            .value()
    };
    assert!((TAU * integrate(|t| rayleigh(t) * t.sin(), 0.0, PI, 10000) - 1.0).abs() < 1e-12);
    for g in [-0.99, -0.65, 0.0, 0.65, 0.99] {
        let model = Hg::try_new(g).unwrap();
        let phase = |t| model.phase(wavelength, Radians::new(t)).value();
        let mass = TAU * integrate(|t| phase(t) * t.sin(), 0.0, PI, 200000);
        let moment = TAU * integrate(|t| phase(t) * t.sin() * t.cos(), 0.0, PI, 200000);
        assert!((mass - 1.0).abs() < 2e-10, "g={g}, mass={mass}");
        assert!((moment - g).abs() < 2e-10, "g={g}, moment={moment}");
        assert!(phase(0.0).is_finite() && phase(PI).is_finite());
        if g > 0.0 {
            assert!(phase(0.0) > phase(PI));
        }
        if g < 0.0 {
            assert!(phase(0.0) < phase(PI));
        }
        if g == 0.0 {
            for t in [0.0, 0.3, 2.0, PI] {
                assert!((phase(t) - 1.0 / (4.0 * PI)).abs() < 1e-15);
            }
        }
    }
    let double = DoubleHenyeyGreensteinPhaseFunction::try_new(0.65, -0.3, 0.8).unwrap();
    let mass = TAU
        * integrate(
            |t| double.phase(wavelength, Radians::new(t)).value() * t.sin(),
            0.0,
            PI,
            10000,
        );
    assert!((mass - 1.0).abs() < 1e-12);
}

#[test]
fn invalid_analytic_parameters_are_rejected_by_optica() {
    for g in [
        -2.0,
        -1.0,
        1.0,
        2.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(Hg::try_new(g).is_err());
        assert!(DoubleHenyeyGreensteinPhaseFunction::try_new(g, 0.0, 0.5).is_err());
    }
}

#[test]
fn isotropic_endpoints_strict_queries_and_legacy_clamping() {
    let p = 1.0 / (4.0 * PI);
    let t = table(&[400.0, 700.0], &[0.0, PI], &[p; 4]).unwrap();
    let cloned = t.clone();
    assert!(format!("{cloned:?}").contains("NormalizedPhaseTable"));
    for x in [400.0, 550.0, 700.0] {
        for y in [0.0, 0.7, PI] {
            assert!((at(&t, x, y).unwrap() - p).abs() < 1e-15);
        }
    }
    for (x, y) in [
        (399.0, 1.0),
        (701.0, 1.0),
        (550.0, -f64::EPSILON),
        (550.0, PI + 1e-15),
    ] {
        let err = at(&t, x, y).unwrap_err();
        assert!(matches!(err, Error::Grid(_)));
        assert!(!err.to_string().is_empty());
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(at(&t, bad, 1.0), Err(Error::NonFiniteQuery)));
        assert!(matches!(at(&t, 550.0, bad), Err(Error::NonFiniteQuery)));
    }
    let legacy = TabulatedPhaseFunction::from_raw_row_major(
        &[400.0, 700.0],
        &[0.0, 180.0],
        &[1.0, 2.0, 3.0, 4.0],
    )
    .unwrap();
    assert_eq!(
        legacy
            .interp_at(Nanometers::new(300.0), Radians::new(-1.0))
            .value(),
        1.0
    );
    assert_eq!(
        legacy
            .interp_at(Nanometers::new(800.0), Radians::new(4.0))
            .value(),
        4.0
    );
}

#[test]
fn rejects_bad_axes_shapes_values_and_normalization() {
    let p = 1.0 / (4.0 * PI);
    for xs in [
        vec![],
        vec![500.0],
        vec![500.0, 500.0],
        vec![500.0, f64::NAN],
        vec![500.0, f64::INFINITY],
    ] {
        assert!(matches!(
            table(&xs, &[0.0, PI], &[p; 4]),
            Err(Error::Grid(_))
        ));
    }
    for xs in [[0.0, 500.0], [-1.0, 500.0], [700.0, 400.0]] {
        assert!(matches!(
            table(&xs, &[0.0, PI], &[p; 4]),
            Err(Error::WavelengthDomain)
        ));
    }
    for ys in [
        [0.1, PI],
        [0.0, PI - 0.1],
        [PI, 0.0],
        [-0.1, PI],
        [0.0, PI + 0.1],
    ] {
        assert!(matches!(
            table(&[400.0, 700.0], &ys, &[p; 4]),
            Err(Error::AngularDomain)
        ));
    }
    for ys in [
        vec![],
        vec![0.0],
        vec![0.0, 0.0],
        vec![0.0, f64::NAN],
        vec![0.0, f64::INFINITY],
        vec![0.0, 2.0, 1.0, PI],
    ] {
        assert!(matches!(
            table(&[400.0, 700.0], &ys, &[p; 4]),
            Err(Error::Grid(_))
        ));
    }
    assert!(matches!(
        table(&[400.0, 700.0], &[0.0, PI], &[p; 3]),
        Err(Error::Grid(_))
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let err = table(&[400.0, 700.0], &[0.0, PI], &[bad, p, p, p]).unwrap_err();
        assert!(matches!(err, Error::InvalidValue { index: 0 }));
        assert!(err.to_string().contains("index 0"));
    }
    for bad in [0.0, 1.0, f64::MAX] {
        let err = table(&[400.0, 700.0], &[0.0, PI], &[bad; 4]).unwrap_err();
        assert!(matches!(err, Error::NotNormalized { index: 0 }));
    }
    assert!(matches!(
        table(&[400.0, 700.0], &[0.0, PI], &[p, 1.0, p, 1.0]),
        Err(Error::NotNormalized { index: 1 })
    ));
}

#[test]
fn nonuniform_interpolant_and_wavelength_normalization() {
    let ys = [0.0, 0.03, 0.2, 1.2, 2.0, PI];
    let mut a = [100.0, 20.0, 0.0, 2.0, 0.0, 1.0];
    let mut b = [0.0, 0.0, 1.0, 2.0, 4.0, 3.0];
    for ps in [&mut a, &mut b] {
        let mass = linear_mass(&ys, ps);
        for p in ps {
            *p /= mass;
        }
    }
    let data: Vec<_> = a.iter().zip(b).flat_map(|(&a, b)| [a, b]).collect();
    let t = table(&[400.0, 700.0], &ys, &data).unwrap();
    // A hand-computed bilinear query with distinct angular/wavelength fractions.
    let expected = 0.75 * (a[2] * 0.7 + a[3] * 0.3) + 0.25 * (b[2] * 0.7 + b[3] * 0.3);
    assert!((at(&t, 475.0, 0.5).unwrap() - expected).abs() < 1e-15);
    for x in [400.0, 475.0, 550.0, 700.0] {
        let ps: Vec<_> = ys.iter().map(|&y| at(&t, x, y).unwrap()).collect();
        assert!((linear_mass(&ys, &ps) - 1.0).abs() < 1e-12);
    }
    // A coarse Rayleigh sample is not normalized as a *linear* interpolant.
    let ps: Vec<_> = [0.0, PI / 2.0, PI]
        .iter()
        .flat_map(|&y| {
            [RayleighPhaseFunction
                .phase(Nanometers::new(500.0), Radians::new(y))
                .value(); 2]
        })
        .collect();
    assert!(matches!(
        table(&[400.0, 700.0], &[0.0, PI / 2.0, PI], &ps),
        Err(Error::NotNormalized { .. })
    ));
}

#[test]
fn well_resolved_analytic_tables() {
    for model in [
        &RayleighPhaseFunction as &dyn PhaseFunction,
        &Hg::try_new(0.65).unwrap(),
    ] {
        let ys: Vec<_> = (0..=20000).map(|i| PI * i as f64 / 20000.0).collect();
        let data: Vec<_> = ys
            .iter()
            .flat_map(|&y| [model.phase(Nanometers::new(500.0), Radians::new(y)).value(); 2])
            .collect();
        let t = table(&[400.0, 700.0], &ys, &data).unwrap();
        assert!(at(&t, 550.0, 0.01).unwrap().is_finite());
    }
}

#[test]
fn narrow_segments_peaks_zero_samples_and_tolerance() {
    let p = 1.0 / (4.0 * PI);
    for ys in [
        vec![0.0, 1e-150, 0.001, 1.0, PI - 1e-14, PI],
        vec![0.0, 1e-12, PI],
    ] {
        table(&[f64::MIN_POSITIVE, f64::MAX], &ys, &vec![p; 2 * ys.len()]).unwrap();
    }
    // A triangular lobe of width d at the pole has integral d²/6 + O(d⁴).
    // Its height is chosen independently from this limiting analytical mass.
    for d in [1e-6, 1e-50, 1e-150] {
        let peak = (3.0 / PI / d) / d;
        let t = table(
            &[400.0, 700.0],
            &[0.0, d, PI],
            &[peak, peak, 0.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(at(&t, 550.0, 0.0).unwrap(), peak);
        assert_eq!(at(&t, 550.0, PI).unwrap(), 0.0);
    }
    // Resolve reflected endpoint weights even at the last representable interval.
    for d in [1e-6, 1e-10, f64::EPSILON * 2.0] {
        let start = PI - d;
        let width = PI - start;
        let peak = (3.0 / PI / width) / width;
        let t = table(
            &[400.0, 700.0],
            &[0.0, start, PI],
            &[0.0, 0.0, 0.0, 0.0, peak, peak],
        )
        .unwrap();
        assert_eq!(at(&t, 550.0, PI).unwrap(), peak);
    }
    for deviation in [-0.999e-6, 0.999e-6] {
        table(&[400.0, 700.0], &[0.0, PI], &[p * (1.0 + deviation); 4]).unwrap();
    }
    for deviation in [-1.001e-6, 1.001e-6] {
        assert!(matches!(
            table(&[400.0, 700.0], &[0.0, PI], &[p * (1.0 + deviation); 4]),
            Err(Error::NotNormalized { .. })
        ));
    }
}
