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
//! a top-of-atmosphere disc-integrated lunar spectral irradiance computed from
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
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};

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

/// Spectral irradiance unit, W m⁻² nm⁻¹, expressed using qtty unit algebra.
pub type LunarSpectralIrradianceUnit = crate::qtty::Per<
    crate::qtty::Per<crate::qtty::unit::Watt, crate::qtty::unit::SquareMeter>,
    crate::qtty::unit::Nanometer,
>;

/// Disc-integrated top-of-atmosphere lunar spectral irradiance.
///
/// Jones et al. (2013), Eq. (1), uses a *fixed reference* lunar solid angle
/// and one distance correction. Equivalently, with Ω(D) ≃ π(R/D)²,
/// `E_moon = E_sun A(α,λ) Ω(D)/π = E_sun A(α,λ) (R/D)²`.
/// The disc-average radiance `E_sun A/π` is distance invariant; integrating
/// it over the apparent disc introduces the inverse-square dilution once.
/// This quantity illuminates the scattering column as a point-source beam.
/// A downstream phase function per steradian converts it to sky radiance;
/// do not multiply by the lunar solid angle a second time.
///
/// The smooth albedo and phase approximations remain those of this module,
/// rather than the full ROLO fit. Solar irradiance is held fixed at the Moon
/// (Earth's value is the approximation); solar-distance variation is external.
///
/// # Arguments
/// * `solar_irradiance_w_m2_nm` — incident solar irradiance in W m⁻² nm⁻¹.
/// * `wavelength` — photon wavelength.
/// * `phase_angle` — Sun–Moon–observer phase angle.
/// * `moon_distance` — observer–Moon distance, positive.
///
/// # Returns
/// Irradiance in W m⁻² nm⁻¹; NaN for nonfinite inputs or nonpositive distance.
pub fn reflected_lunar_spectral_irradiance_jones2013(
    solar_irradiance_w_m2_nm: f64,
    wavelength: Nanometers,
    phase_angle: Radians,
    moon_distance: Kilometers,
) -> crate::qtty::Quantity<LunarSpectralIrradianceUnit> {
    let distance = moon_distance.value();
    if !solar_irradiance_w_m2_nm.is_finite()
        || !wavelength.is_finite()
        || !phase_angle.is_finite()
        || !distance.is_finite()
        || distance <= 0.0
    {
        return crate::qtty::Quantity::new(f64::NAN);
    }
    let albedo = lunar_albedo_jones2013(phase_angle, wavelength).value();
    let radius_ratio: f64 = MEAN_MOON_RADIUS / moon_distance;
    crate::qtty::Quantity::new(solar_irradiance_w_m2_nm * radius_ratio.powi(2) * albedo)
}

/// Legacy compatibility wrapper for lunar beam irradiance.
///
/// **Historical unit error:** despite its name and radiance unit marker,
/// this helper's numeric value is disc-integrated irradiance in W m⁻² nm⁻¹.
/// Use [`reflected_lunar_spectral_irradiance_jones2013`] for a dimensionally
/// correct typed API. The old signature is retained for source compatibility;
/// it does not describe the surface radiance or an additional per-sr source.
/// The erroneous second inverse-square dilution has been removed.
/// Arguments and invalid-input behavior follow the irradiance helper.
pub fn reflected_lunar_spectral_radiance_jones2013(
    solar_irradiance_w_m2_nm: f64,
    wavelength: Nanometers,
    phase_angle: Radians,
    moon_distance: Kilometers,
) -> WattsPerSquareMeterSteradianNanometer {
    WattsPerSquareMeterSteradianNanometer::new(
        reflected_lunar_spectral_irradiance_jones2013(
            solar_irradiance_w_m2_nm,
            wavelength,
            phase_angle,
            moon_distance,
        )
        .value(),
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
    #[test]
    fn beam_irradiance_obeys_disc_geometry_and_inverse_square_law() {
        let wave = Nanometers::new(551.0);
        let phase = Degrees::new(30.0).to::<Radian>();
        let albedo = lunar_albedo_jones2013(phase, wave).value();
        let reference =
            reflected_lunar_spectral_irradiance_jones2013(1.8, wave, phase, MEAN_MOON_DISTANCE);
        for d in [350_000.0, 384_400.0, 405_000.0, 768_800.0] {
            let actual =
                reflected_lunar_spectral_irradiance_jones2013(1.8, wave, phase, Kilometers::new(d));
            // Independent disc-radiance × solid-angle construction.
            let disc_radiance = 1.8 * albedo / core::f64::consts::PI;
            let omega = core::f64::consts::PI * (1_737.4 / d).powi(2);
            assert!((actual.value() / (disc_radiance * omega) - 1.0).abs() < 1e-14);
            assert!((actual.value() / reference.value() - (384_400.0 / d).powi(2)).abs() < 1e-14);
            let legacy =
                reflected_lunar_spectral_radiance_jones2013(1.8, wave, phase, Kilometers::new(d));
            assert_eq!(legacy.value(), actual.value());
        }
        // Compile-time dimensional check: irradiance integrated over wavelength
        // is a power per area, not a radiance per steradian.
        let band = reference * Nanometers::new(10.0);
        assert!((band.value() / reference.value() - 10.0).abs() < 1e-14);
    }

    #[test]
    fn lunar_beam_zero_illumination_and_invalid_inputs() {
        let wave = Nanometers::new(550.0);
        let phase = Radians::new(0.0);
        assert_eq!(
            reflected_lunar_spectral_irradiance_jones2013(0.0, wave, phase, MEAN_MOON_DISTANCE,)
                .value(),
            0.0
        );
        for distance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(reflected_lunar_spectral_irradiance_jones2013(
                1.0,
                wave,
                phase,
                Kilometers::new(distance),
            )
            .value()
            .is_nan());
        }
        assert!(reflected_lunar_spectral_irradiance_jones2013(
            1.0,
            Nanometers::new(f64::NAN),
            phase,
            MEAN_MOON_DISTANCE,
        )
        .value()
        .is_nan());
        assert!(reflected_lunar_spectral_irradiance_jones2013(
            1.0,
            wave,
            Radians::new(f64::NAN),
            MEAN_MOON_DISTANCE,
        )
        .value()
        .is_nan());
    }
}
