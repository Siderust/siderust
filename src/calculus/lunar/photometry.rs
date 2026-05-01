// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lunar photometry helpers used by scattered-moonlight models.

use crate::qtty::radiometry::{
    WattPerSquareMeterSteradianNanometer, WattsPerSquareMeterSteradianNanometer,
};
use crate::qtty::{Kilometers, Nanometers, Radians};

/// Mean lunar radius, in kilometres.
pub const MEAN_MOON_RADIUS_KM: f64 = 1_737.4;

/// Mean geocentric lunar distance, in kilometres.
pub const MEAN_MOON_DISTANCE_KM: f64 = 384_400.0;

/// Approximate wavelength-dependent full-Moon geometric albedo.
///
/// The values are intentionally smooth and conservative: optical lunar
/// reflectance rises from blue to red, and this helper is meant to provide the
/// reusable radiometric backbone for domain models that apply their own
/// calibration and scattering terms.
pub fn lunar_full_moon_albedo_jones2013(wavelength: Nanometers) -> f64 {
    let lambda = wavelength.value().clamp(300.0, 1_100.0);
    let t = (lambda - 300.0) / 800.0;
    0.075 + 0.065 * t
}

/// Phase attenuation used by Jones/Noll-style scattered moonlight models.
///
/// This follows the Krisciunas-Schaefer phase-angle polynomial also used by
/// the Jones et al. workflow as the low-order lunar phase brightness law, but
/// normalized to 1 at full Moon by excluding the constant full-Moon offset.
pub fn lunar_phase_attenuation_jones2013(phase_angle: Radians) -> f64 {
    let a = phase_angle.value().abs().to_degrees();
    10f64.powf(-0.4 * (0.026 * a + 4.0e-9 * a.powi(4)))
}

/// Wavelength- and phase-dependent lunar albedo factor.
pub fn lunar_albedo_jones2013(phase_angle: Radians, wavelength: Nanometers) -> f64 {
    lunar_full_moon_albedo_jones2013(wavelength) * lunar_phase_attenuation_jones2013(phase_angle)
}

/// Top-of-atmosphere lunar spectral radiance from incident solar irradiance.
///
/// `solar_irradiance_w_m2_nm` is the solar spectral irradiance at Earth in
/// `W·m⁻²·nm⁻¹`. The output is a lunar spectral radiance in
/// `W·m⁻²·sr⁻¹·nm⁻¹`, scaled by lunar solid angle and actual Moon distance.
pub fn reflected_lunar_spectral_radiance_jones2013(
    solar_irradiance_w_m2_nm: f64,
    wavelength: Nanometers,
    phase_angle: Radians,
    moon_distance: Kilometers,
) -> WattsPerSquareMeterSteradianNanometer {
    let distance = moon_distance.value();
    if !solar_irradiance_w_m2_nm.is_finite() || !distance.is_finite() || distance <= 0.0 {
        return WattsPerSquareMeterSteradianNanometer::new(f64::NAN);
    }
    let albedo = lunar_albedo_jones2013(phase_angle, wavelength);
    let omega_over_pi = (MEAN_MOON_RADIUS_KM / distance).powi(2);
    let distance_scale = (MEAN_MOON_DISTANCE_KM / distance).powi(2);
    WattsPerSquareMeterSteradianNanometer::new(
        solar_irradiance_w_m2_nm * omega_over_pi * albedo * distance_scale,
    )
}

/// Unit marker for the output of [`reflected_lunar_spectral_radiance_jones2013`].
pub type LunarSpectralRadianceUnit = WattPerSquareMeterSteradianNanometer;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{unit::Radian, Degrees};

    #[test]
    fn albedo_increases_from_blue_to_red() {
        let blue = lunar_full_moon_albedo_jones2013(Nanometers::new(400.0));
        let red = lunar_full_moon_albedo_jones2013(Nanometers::new(800.0));
        assert!(red > blue);
        assert!(blue > 0.0);
    }

    #[test]
    fn phase_attenuation_is_one_at_full_and_small_near_new() {
        let full = lunar_phase_attenuation_jones2013(Degrees::new(0.0).to::<Radian>());
        let new = lunar_phase_attenuation_jones2013(Degrees::new(180.0).to::<Radian>());
        assert!((full - 1.0).abs() < 1.0e-12);
        assert!(new < 1.0e-3);
    }

    #[test]
    fn reflected_radiance_is_positive_for_full_moon() {
        let l = reflected_lunar_spectral_radiance_jones2013(
            1.8,
            Nanometers::new(550.0),
            Degrees::new(0.0).to::<Radian>(),
            Kilometers::new(MEAN_MOON_DISTANCE_KM),
        );
        assert!(l.value() > 0.0);
    }
}
