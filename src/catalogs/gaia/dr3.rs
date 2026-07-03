// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed Gaia DR3 source domain model.

use super::error::{GaiaDr3Error, Result};
use super::photometry::{PhotonFlux, StellarPhotometryModel};
use super::quality::GaiaDr3QualityFlags;
use super::raw::GaiaDr3RawSourceRow;
use crate::astro::proper_motion::ProperMotion;
use crate::qtty::{
    angular_rate::AngularRate, unit, Degrees, JulianYears, KmPerSeconds, MilliArcsecond,
    MilliArcseconds,
};

/// Gaia proper-motion rate in milliarcseconds per Julian year.
type GaiaProperMotionRate = AngularRate<MilliArcsecond, unit::JulianYear>;

/// Gaia source identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GaiaSourceId(u64);

impl GaiaSourceId {
    /// Build a Gaia source identifier.
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Return the raw Gaia source identifier.
    pub const fn value(self) -> u64 {
        self.0
    }
}

/// Validated Gaia DR3 astrometry.
#[derive(Debug, Clone)]
pub struct GaiaDr3Astrometry {
    /// ICRS spherical direction at [`Self::epoch`].
    pub direction: crate::coordinates::spherical::direction::ICRS,
    /// Gaia reference epoch in Julian years (`ref_epoch_jyr`).
    pub epoch: JulianYears,
    /// Proper motion, when both Gaia components are present.
    ///
    /// Gaia publishes `pmra` as `µα⋆ = µα cos(δ)` and `pmdec` as `µδ`, both as
    /// angular rates per Julian year.
    pub proper_motion: Option<ProperMotion>,
    /// Annual trigonometric parallax.
    pub parallax: Option<MilliArcseconds>,
    /// Radial velocity.
    pub radial_velocity: Option<KmPerSeconds>,
}

/// Gaia DR3 broad-band photometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaiaDr3Photometry {
    /// Gaia G mean magnitude.
    pub phot_g_mean_mag: Option<f64>,
    /// Gaia BP mean magnitude.
    pub phot_bp_mean_mag: Option<f64>,
    /// Gaia RP mean magnitude.
    pub phot_rp_mean_mag: Option<f64>,
}

/// Validated Gaia DR3 source.
#[derive(Debug, Clone)]
pub struct GaiaDr3Source {
    /// Gaia source identifier.
    pub source_id: GaiaSourceId,
    /// Validated typed astrometry.
    pub astrometry: GaiaDr3Astrometry,
    /// Gaia DR3 broad-band photometry.
    pub photometry: GaiaDr3Photometry,
    /// Quality metadata.
    pub quality: GaiaDr3QualityFlags,
}

impl TryFrom<GaiaDr3RawSourceRow> for GaiaDr3Source {
    type Error = GaiaDr3Error;

    fn try_from(row: GaiaDr3RawSourceRow) -> Result<Self> {
        ensure_finite("ra_deg", row.ra_deg)?;
        ensure_finite("dec_deg", row.dec_deg)?;
        if !(0.0..360.0).contains(&row.ra_deg) {
            return Err(GaiaDr3Error::RightAscensionOutOfRange { value: row.ra_deg });
        }
        if !(-90.0..=90.0).contains(&row.dec_deg) {
            return Err(GaiaDr3Error::DeclinationOutOfRange { value: row.dec_deg });
        }

        let proper_motion = match (row.pmra_mas_per_yr, row.pmdec_mas_per_yr) {
            (Some(pmra), Some(pmdec)) => Some(ProperMotion::from_mu_alpha_star(
                gaia_proper_motion_rate("pmra", pmra)?,
                gaia_proper_motion_rate("pmdec", pmdec)?,
            )),
            (None, None) => None,
            (Some(_), None) => return Err(GaiaDr3Error::MissingRequiredField { field: "pmdec" }),
            (None, Some(_)) => return Err(GaiaDr3Error::MissingRequiredField { field: "pmra" }),
        };

        let photometry = GaiaDr3Photometry {
            phot_g_mean_mag: validate_optional_finite("phot_g_mean_mag", row.phot_g_mean_mag)?,
            phot_bp_mean_mag: validate_optional_finite("phot_bp_mean_mag", row.phot_bp_mean_mag)?,
            phot_rp_mean_mag: validate_optional_finite("phot_rp_mean_mag", row.phot_rp_mean_mag)?,
        };

        Ok(Self {
            source_id: GaiaSourceId::new(row.source_id),
            astrometry: GaiaDr3Astrometry {
                direction: crate::coordinates::spherical::Direction::<
                    crate::coordinates::frames::ICRS,
                >::new(
                    Degrees::new(row.ra_deg), Degrees::new(row.dec_deg)
                ),
                epoch: julian_year_epoch("ref_epoch_jyr", row.ref_epoch_jyr)?,
                proper_motion,
                parallax: validate_optional_finite("parallax_mas", row.parallax_mas)?
                    .map(MilliArcseconds::new),
                radial_velocity: validate_optional_finite(
                    "radial_velocity_km_s",
                    row.radial_velocity_km_s,
                )?
                .map(KmPerSeconds::new),
            },
            photometry,
            quality: row.quality,
        })
    }
}

/// Passband-integrated source record suitable for later sky accumulation.
#[derive(Debug, Clone)]
pub struct PassbandIntegratedStellarSource {
    /// Catalogue source identifier.
    pub source_id: u64,
    /// ICRS direction of the source.
    pub direction: crate::coordinates::spherical::direction::ICRS,
    /// Source reference epoch in Julian years.
    pub epoch: JulianYears,
    /// Passband-integrated photon flux.
    pub photon_flux: PhotonFlux,
    /// Photometric model used to produce [`Self::photon_flux`].
    pub photometry_model: StellarPhotometryModel,
}

impl PassbandIntegratedStellarSource {
    /// Build a passband-integrated source from typed Gaia DR3 astrometry.
    pub fn from_gaia_source(
        source: &GaiaDr3Source,
        photon_flux: PhotonFlux,
        photometry_model: StellarPhotometryModel,
    ) -> Self {
        Self {
            source_id: source.source_id.value(),
            direction: source.astrometry.direction,
            epoch: source.astrometry.epoch,
            photon_flux,
            photometry_model,
        }
    }
}

fn julian_year_epoch(field: &'static str, value: f64) -> Result<JulianYears> {
    ensure_finite(field, value)?;
    Ok(JulianYears::new(value))
}

fn gaia_proper_motion_rate(field: &'static str, value: f64) -> Result<GaiaProperMotionRate> {
    ensure_finite(field, value)?;
    Ok(GaiaProperMotionRate::new(value))
}

fn validate_optional_finite(field: &'static str, value: Option<f64>) -> Result<Option<f64>> {
    if let Some(value) = value {
        ensure_finite(field, value)?;
    }
    Ok(value)
}

fn ensure_finite(field: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(GaiaDr3Error::NonFinite { field, value })
    }
}
