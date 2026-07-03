// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Passband-aware Gaia DR3 photometry primitives.
//!
//! Gaia XP sampled spectra are represented as spectral flux densities in
//! `W m^-2 nm^-1`. Photon-flux integration converts to SI internally and
//! returns integrated photon flux in `photons m^-2 s^-1`.

use super::error::{GaiaDr3Error, Result};
use crate::qtty::velocity::C;
use crate::qtty::Nanometers;

const PLANCK_J_S: f64 = 6.626_070_15e-34;

fn wavelength_nm(field: &'static str, value: f64) -> Result<Nanometers> {
    if !value.is_finite() {
        return Err(GaiaDr3Error::NonFinite { field, value });
    }
    if value <= 0.0 {
        return Err(GaiaDr3Error::NotPositive { field, value });
    }
    Ok(Nanometers::new(value))
}

/// Spectral flux density in `W m^-2 nm^-1`.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SpectralFluxDensity(f64);

impl SpectralFluxDensity {
    /// Build a finite non-negative spectral flux density in `W m^-2 nm^-1`.
    pub fn new_w_m2_nm(value: f64) -> Result<Self> {
        if !value.is_finite() {
            return Err(GaiaDr3Error::NonFinite {
                field: "flux_density_w_m2_nm",
                value,
            });
        }
        if value < 0.0 {
            return Err(GaiaDr3Error::Negative {
                field: "flux_density_w_m2_nm",
                value,
            });
        }
        Ok(Self(value))
    }

    /// Return the value in `W m^-2 nm^-1`.
    pub const fn w_m2_nm(self) -> f64 {
        self.0
    }

    /// Return the value in `W m^-2 m^-1`.
    pub fn w_m2_m(self) -> f64 {
        self.0 * 1.0e9
    }
}

/// Integrated photon flux in `photons m^-2 s^-1`.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct PhotonFlux(f64);

impl PhotonFlux {
    /// Build a finite non-negative photon flux in `photons m^-2 s^-1`.
    pub fn new_photons_m2_s(value: f64) -> Result<Self> {
        if !value.is_finite() {
            return Err(GaiaDr3Error::NonFinite {
                field: "photon_flux_ph_m2_s",
                value,
            });
        }
        if value < 0.0 {
            return Err(GaiaDr3Error::Negative {
                field: "photon_flux_ph_m2_s",
                value,
            });
        }
        Ok(Self(value))
    }

    /// Return the value in `photons m^-2 s^-1`.
    pub const fn photons_m2_s(self) -> f64 {
        self.0
    }
}

/// Closed wavelength interval used for passband integration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralBand {
    /// Stable band name.
    pub name: &'static str,
    /// Minimum wavelength.
    pub min_wavelength: Nanometers,
    /// Maximum wavelength.
    pub max_wavelength: Nanometers,
}

impl SpectralBand {
    /// Build a spectral band from nanometer bounds.
    pub fn new(name: &'static str, min_wavelength_nm: f64, max_wavelength_nm: f64) -> Result<Self> {
        let band = Self {
            name,
            min_wavelength: wavelength_nm("min_wavelength_nm", min_wavelength_nm)?,
            max_wavelength: wavelength_nm("max_wavelength_nm", max_wavelength_nm)?,
        };
        band.validate()?;
        Ok(band)
    }

    fn validate(self) -> Result<()> {
        if self.min_wavelength.value() >= self.max_wavelength.value() {
            return Err(GaiaDr3Error::InvalidSpectralBand { name: self.name });
        }
        Ok(())
    }
}

/// Gaia XP starlight production band, 330-650 nm.
pub const GAIA_XP_STELLAR_RADIANCE_330_650_NM: SpectralBand = SpectralBand {
    name: "Gaia XP stellar radiance 330-650 nm",
    min_wavelength: Nanometers::new(330.0),
    max_wavelength: Nanometers::new(650.0),
};

/// One sampled spectral-flux-density point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralFluxSample {
    /// Sample wavelength.
    pub wavelength: Nanometers,
    /// Spectral flux density at the sample wavelength.
    pub flux_density: SpectralFluxDensity,
}

impl SpectralFluxSample {
    /// Build a validated sample from `nm` and `W m^-2 nm^-1`.
    pub fn new(wavelength_nm_value: f64, flux_density_w_m2_nm: f64) -> Result<Self> {
        Ok(Self {
            wavelength: wavelength_nm("wavelength_nm", wavelength_nm_value)?,
            flux_density: SpectralFluxDensity::new_w_m2_nm(flux_density_w_m2_nm)?,
        })
    }
}

/// Gaia XP sampled spectrum with strictly increasing wavelengths.
#[derive(Debug, Clone, PartialEq)]
pub struct GaiaXpSampledSpectrum {
    /// Validated samples.
    pub samples: Vec<SpectralFluxSample>,
}

impl GaiaXpSampledSpectrum {
    /// Build and validate a Gaia XP sampled spectrum.
    pub fn new(samples: Vec<SpectralFluxSample>) -> Result<Self> {
        if samples.is_empty() {
            return Err(GaiaDr3Error::EmptySpectrum);
        }
        for pair in samples.windows(2) {
            if pair[0].wavelength.value() >= pair[1].wavelength.value() {
                return Err(GaiaDr3Error::NonMonotonicSpectrum);
            }
        }
        Ok(Self { samples })
    }
}

/// Photometry model metadata for starlight source generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StellarPhotometryModel {
    /// Gaia DR3 XP spectrum integrated as photon flux over 330-650 nm.
    GaiaDr3XpPhotonRadiance330650NmV1,
    /// Gaia DR3 G/BP/RP SED-derived proxy over 330-650 nm.
    GaiaDr3GbprpSedPhotonRadiance330650NmV1,
    /// Tycho BT/VT to Johnson proxy model.
    TychoBtVtJohnsonProxyV1,
}

/// Integrate photon flux over a wavelength band.
///
/// The calculation linearly interpolates the sampled `F_lambda` values and
/// applies trapezoidal integration to
/// `F_lambda(lambda) * lambda / (h * c) d lambda`, using SI units internally.
/// Input `F_lambda` is `W m^-2 nm^-1`; wavelengths are converted to meters.
pub fn integrate_photon_flux(
    spectrum: &GaiaXpSampledSpectrum,
    band: SpectralBand,
) -> Result<PhotonFlux> {
    band.validate()?;

    let c_m_s = C.value();

    let min_nm = band.min_wavelength.value();
    let max_nm = band.max_wavelength.value();
    let mut wavelengths = Vec::with_capacity(spectrum.samples.len() + 2);
    wavelengths.push(min_nm);
    for sample in &spectrum.samples {
        let wavelength = sample.wavelength.value();
        if wavelength > min_nm && wavelength < max_nm {
            wavelengths.push(wavelength);
        }
    }
    wavelengths.push(max_nm);
    wavelengths.sort_by(f64::total_cmp);
    wavelengths.dedup_by(|a, b| (*a - *b).abs() <= f64::EPSILON);

    if wavelengths.len() < 2
        || interpolate_flux_w_m2_nm(spectrum, min_nm).is_none()
        || interpolate_flux_w_m2_nm(spectrum, max_nm).is_none()
    {
        return Err(GaiaDr3Error::InsufficientBandCoverage { name: band.name });
    }

    let mut total = 0.0;
    for pair in wavelengths.windows(2) {
        let l0_nm = pair[0];
        let l1_nm = pair[1];
        let f0 = interpolate_flux_w_m2_nm(spectrum, l0_nm)
            .ok_or(GaiaDr3Error::InsufficientBandCoverage { name: band.name })?;
        let f1 = interpolate_flux_w_m2_nm(spectrum, l1_nm)
            .ok_or(GaiaDr3Error::InsufficientBandCoverage { name: band.name })?;
        let l0_m = l0_nm * 1.0e-9;
        let l1_m = l1_nm * 1.0e-9;
        let y0 = (f0 * 1.0e9) * l0_m / (PLANCK_J_S * c_m_s);
        let y1 = (f1 * 1.0e9) * l1_m / (PLANCK_J_S * c_m_s);
        total += 0.5 * (y0 + y1) * (l1_m - l0_m);
    }

    PhotonFlux::new_photons_m2_s(total)
}

fn interpolate_flux_w_m2_nm(spectrum: &GaiaXpSampledSpectrum, wavelength_nm: f64) -> Option<f64> {
    let first = spectrum.samples.first()?;
    let last = spectrum.samples.last()?;
    if wavelength_nm < first.wavelength.value() || wavelength_nm > last.wavelength.value() {
        return None;
    }
    if (wavelength_nm - first.wavelength.value()).abs() <= f64::EPSILON {
        return Some(first.flux_density.w_m2_nm());
    }
    for pair in spectrum.samples.windows(2) {
        let x0 = pair[0].wavelength.value();
        let x1 = pair[1].wavelength.value();
        if wavelength_nm <= x1 {
            let y0 = pair[0].flux_density.w_m2_nm();
            let y1 = pair[1].flux_density.w_m2_nm();
            let t = (wavelength_nm - x0) / (x1 - x0);
            return Some(y0 + t * (y1 - y0));
        }
    }
    Some(last.flux_density.w_m2_nm())
}
