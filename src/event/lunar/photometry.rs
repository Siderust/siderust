// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Photometry
//!
//! ## Scientific scope
//!
//! Reusable photometric helpers feeding scattered‑moonlight sky‑brightness
//! models (Krisciunas–Schaefer 1991 phase law and Jones et al. 2013
//! workflow). Provides a smooth wavelength‑dependent full‑Moon geometric
//! albedo, the phase‑angle attenuation normalized to 1 at full Moon, and
//! a top‑of‑atmosphere reflected lunar spectral radiance computed from
//! incident solar spectral irradiance and the geometric Moon distance.
//!
//! ## Technical scope
//!
//! Pure functions over `qtty` types; no time dependency. Constants
//! [`MEAN_MOON_RADIUS`] and [`MEAN_MOON_DISTANCE`] anchor the
//! solid‑angle / distance scaling. Output unit alias
//! [`LunarSpectralRadianceUnit`] is exported for downstream pipelines.
//! Implementations are deliberately conservative — domain models are
//! expected to apply their own calibration and atmospheric transfer.
//!
//! ## References
//! - Krisciunas, K. & Schaefer, B. E. (1991). "A model of the brightness
//!   of moonlight". *PASP* 103, 1033. doi:10.1086/132921
//! - Jones, A., Noll, S., Kausch, W., Szyszka, C., & Kimeswenger, S.
//!   (2013). "An advanced scattered moonlight model for Cerro Paranal".
//!   *A&A* 560, A91. doi:10.1051/0004-6361/201322433

use crate::qtty::radiometry::{
    WattPerSquareMeterSteradianNanometer, WattsPerSquareMeterSteradianNanometer,
};
use crate::qtty::{Albedos, IlluminationFractions};
use crate::qtty::{Kilometers, Nanometers, Radians};

/// Mean lunar radius.
pub const MEAN_MOON_RADIUS: Kilometers = Kilometers::new(1_737.4);

/// Mean geocentric lunar distance.
pub const MEAN_MOON_DISTANCE: Kilometers = Kilometers::new(384_400.0);

/// Approximate wavelength-dependent full-Moon geometric albedo.
///
/// The values are intentionally smooth and conservative: optical lunar
/// reflectance rises from blue to red, and this helper is meant to provide the
/// reusable radiometric backbone for domain models that apply their own
/// calibration and scattering terms.
///
/// # Arguments
///
/// * `wavelength`, photon wavelength; clamped internally to `[300, 1100] nm`.
///
/// # Returns
///
/// Dimensionless full‑Moon geometric albedo at `wavelength`.
pub fn lunar_full_moon_albedo_jones2013(wavelength: Nanometers) -> Albedos {
    let lambda = wavelength
        .clamp(Nanometers::new(300.0), Nanometers::new(1_100.0))
        .value();
    let t = (lambda - 300.0) / 800.0;
    Albedos::new(0.075 + 0.065 * t)
}

/// Phase attenuation used by Jones/Noll-style scattered moonlight models.
///
/// This follows the Krisciunas-Schaefer phase-angle polynomial also used by
/// the Jones et al. workflow as the low-order lunar phase brightness law, but
/// normalized to 1 at full Moon by excluding the constant full-Moon offset.
///
/// # Arguments
///
/// * `phase_angle`, Sun–Moon–observer phase angle; only `|α|` is used.
///
/// # Returns
///
/// Dimensionless attenuation factor (typed [`IlluminationFractions`]),
/// equal to `1` at full Moon and vanishing near the new Moon.
pub fn lunar_phase_attenuation_jones2013(phase_angle: Radians) -> IlluminationFractions {
    let a = phase_angle.abs().value().to_degrees();
    IlluminationFractions::new(10f64.powf(-0.4 * (0.026 * a + 4.0e-9 * a.powi(4))))
}

/// Wavelength- and phase-dependent lunar albedo factor.
///
/// # Arguments
///
/// * `phase_angle`, Sun–Moon–observer phase angle.
/// * `wavelength`, photon wavelength.
///
/// # Returns
///
/// Product of [`lunar_full_moon_albedo_jones2013`] and
/// [`lunar_phase_attenuation_jones2013`].
pub fn lunar_albedo_jones2013(phase_angle: Radians, wavelength: Nanometers) -> Albedos {
    let a = lunar_full_moon_albedo_jones2013(wavelength).value();
    let p = lunar_phase_attenuation_jones2013(phase_angle).value();
    Albedos::new(a * p)
}

/// Top-of-atmosphere lunar spectral radiance from incident solar irradiance.
///
/// `solar_irradiance_w_m2_nm` is the solar spectral irradiance at Earth in
/// `W·m⁻²·nm⁻¹` (kept as a raw `f64` because the matching spectral-irradiance
/// unit is not yet exposed by `qtty`). The output is a typed lunar spectral
/// radiance in `W·m⁻²·sr⁻¹·nm⁻¹`, scaled by lunar solid angle and actual
/// Moon distance.
///
/// # Arguments
///
/// * `solar_irradiance_w_m2_nm`, solar spectral irradiance at Earth, `W·m⁻²·nm⁻¹`.
/// * `wavelength`, photon wavelength.
/// * `phase_angle`, Sun–Moon–observer phase angle.
/// * `moon_distance`, geometric Earth–Moon distance.
///
/// # Returns
///
/// Reflected lunar spectral radiance at the top of the atmosphere; `NaN`
/// if any input is non‑finite or `moon_distance ≤ 0`.
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
    let albedo = lunar_albedo_jones2013(phase_angle, wavelength).value();
    let radius_ratio: f64 = MEAN_MOON_RADIUS / moon_distance;
    let omega_over_pi = radius_ratio.powi(2);
    let distance_ratio: f64 = MEAN_MOON_DISTANCE / moon_distance;
    let distance_scale = distance_ratio.powi(2);
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
        assert!(red.value() > blue.value());
        assert!(blue.value() > 0.0);
    }

    #[test]
    fn phase_attenuation_is_one_at_full_and_small_near_new() {
        let full = lunar_phase_attenuation_jones2013(Degrees::new(0.0).to::<Radian>());
        let new = lunar_phase_attenuation_jones2013(Degrees::new(180.0).to::<Radian>());
        assert!((full.value() - 1.0).abs() < 1.0e-12);
        assert!(new.value() < 1.0e-3);
    }

    #[test]
    fn reflected_radiance_is_positive_for_full_moon() {
        let l = reflected_lunar_spectral_radiance_jones2013(
            1.8,
            Nanometers::new(550.0),
            Degrees::new(0.0).to::<Radian>(),
            MEAN_MOON_DISTANCE,
        );
        assert!(l.value() > 0.0);
    }
}
