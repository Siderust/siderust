// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

#[test]
fn apparent_magnitude_accepts_signed_finite_values_and_rejects_nonfinite() {
    for input in [-100.0, -0.0, 0.0, 10.0, 1.0e6] {
        let magnitude = ApparentMagnitude::new(input).expect("finite magnitude");
        assert_eq!(magnitude.value(), input);
    }

    for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            ApparentMagnitude::new(input),
            Err(StellarMapError::InvalidMagnitude(value)) if !value.is_finite()
        ));
    }
}

#[test]
fn s10_flux_obeys_the_astronomical_magnitude_law() {
    for (magnitude, expected) in [
        (0.0, 10_000.0),
        (5.0, 100.0),
        (10.0, 1.0),
        (15.0, 0.01),
        (20.0, 0.0001),
    ] {
        let actual = flux_10mag_units(ApparentMagnitude::new(magnitude).expect("magnitude"));
        assert!(
            (actual / expected - 1.0).abs() < 1e-12,
            "magnitude {magnitude}: expected {expected}, got {actual}"
        );
    }

    let faint = flux_10mag_units(ApparentMagnitude::new(12.5).expect("faint"));
    let bright = flux_10mag_units(ApparentMagnitude::new(10.0).expect("bright"));
    assert!((bright / faint - 10.0).abs() < 1e-12);
}

#[test]
fn brightness_validation_checks_every_component() {
    assert!(StellarSurfaceBrightness::zero().is_valid());

    let valid = StellarSurfaceBrightness {
        integrated_ph_cm2_ns_sr: 1.0,
        b_s10: 2.0,
        v_s10: 3.0,
    };
    assert!(valid.is_valid());

    for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(!StellarSurfaceBrightness {
            integrated_ph_cm2_ns_sr: invalid,
            ..valid
        }
        .is_valid());
        assert!(!StellarSurfaceBrightness {
            b_s10: invalid,
            ..valid
        }
        .is_valid());
        assert!(!StellarSurfaceBrightness {
            v_s10: invalid,
            ..valid
        }
        .is_valid());
    }
}
