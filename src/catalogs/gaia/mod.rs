// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Gaia catalogue support.
//!
//! Raw Gaia DR3 types in [`raw`] preserve primitive catalogue units at
//! CSV/ADQL boundaries. Domain types in [`dr3`] validate those rows and expose
//! typed quantities plus a single typed ICRS direction field.

pub mod dr3;
pub mod error;
pub mod photometry;
pub mod quality;
pub mod raw;

pub use dr3::{
    GaiaDr3Astrometry, GaiaDr3Photometry, GaiaDr3Source, GaiaSourceId,
    PassbandIntegratedStellarSource,
};
pub use error::{GaiaDr3Error, Result};
pub use photometry::{
    integrate_photon_flux, GaiaXpSampledSpectrum, Nanometers, PhotonFlux, SpectralBand,
    SpectralFluxDensity, SpectralFluxSample, StellarPhotometryModel,
    GAIA_XP_STELLAR_RADIANCE_330_650_NM,
};
pub use quality::GaiaDr3QualityFlags;
pub use raw::{parse_gaia_dr3_csv_chunk, GaiaDr3RawSourceRow};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::proper_motion::{ProperMotion, RaProperMotionConvention};
    use crate::qtty::{angular_rate::AngularRate, unit, MilliArcsecond};

    type GaiaProperMotionRate = AngularRate<MilliArcsecond, unit::JulianYear>;
    type DegreesPerJulianYear = AngularRate<crate::qtty::unit::Degree, unit::JulianYear>;

    fn raw_row() -> GaiaDr3RawSourceRow {
        GaiaDr3RawSourceRow {
            source_id: 42,
            ra_deg: 120.0,
            dec_deg: -30.0,
            ref_epoch_jyr: 2016.0,
            pmra_mas_per_yr: Some(1.5),
            pmdec_mas_per_yr: Some(-2.5),
            parallax_mas: Some(3.0),
            radial_velocity_km_s: Some(20.0),
            phot_g_mean_mag: Some(12.0),
            phot_bp_mean_mag: Some(12.5),
            phot_rp_mean_mag: Some(11.5),
            quality: GaiaDr3QualityFlags::default(),
        }
    }

    #[test]
    fn raw_gaia_row_becomes_typed_source() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");

        assert_eq!(source.source_id.value(), 42);
        assert_eq!(source.astrometry.epoch.value(), 2016.0);
        let pm = source.astrometry.proper_motion.expect("proper motion");
        assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
        let pmra_mas: GaiaProperMotionRate = pm.pm_ra.to();
        let pmdec_mas: GaiaProperMotionRate = pm.pm_dec.to();
        assert!((pmra_mas.value() - 1.5).abs() < 1e-12);
        assert!((pmdec_mas.value() + 2.5).abs() < 1e-12);
        assert_eq!(source.astrometry.parallax.expect("parallax").value(), 3.0);
        assert_eq!(
            source
                .astrometry
                .radial_velocity
                .expect("radial velocity")
                .value(),
            20.0
        );
    }

    #[test]
    fn invalid_ra_dec_are_rejected() {
        let mut row = raw_row();
        row.ra_deg = 360.0;
        assert!(matches!(
            GaiaDr3Source::try_from(row),
            Err(GaiaDr3Error::RightAscensionOutOfRange { .. })
        ));

        let mut row = raw_row();
        row.dec_deg = -91.0;
        assert!(matches!(
            GaiaDr3Source::try_from(row),
            Err(GaiaDr3Error::DeclinationOutOfRange { .. })
        ));
    }

    #[test]
    fn degrees_are_converted_into_spherical_icrs_direction() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");

        assert!((source.astrometry.direction.azimuth.value() - 120.0).abs() < 1e-12);
        assert!((source.astrometry.direction.polar.value() + 30.0).abs() < 1e-12);
    }

    #[test]
    fn gaia_domain_source_has_single_direction_field() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");

        let _: crate::coordinates::spherical::direction::ICRS = source.astrometry.direction;
    }

    #[test]
    fn ref_epoch_jyr_converts_to_julian_years() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");

        assert_eq!(source.astrometry.epoch.value(), 2016.0);
    }

    #[test]
    fn non_finite_ref_epoch_jyr_is_rejected() {
        let mut row = raw_row();
        row.ref_epoch_jyr = f64::NAN;
        assert!(matches!(
            GaiaDr3Source::try_from(row),
            Err(GaiaDr3Error::NonFinite {
                field: "ref_epoch_jyr",
                ..
            })
        ));
    }

    #[test]
    fn proper_motion_uses_mu_alpha_star_convention() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");
        let pm = source.astrometry.proper_motion.expect("proper motion");

        assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
    }

    #[test]
    fn missing_one_proper_motion_component_is_rejected() {
        let mut row = raw_row();
        row.pmdec_mas_per_yr = None;
        assert!(matches!(
            GaiaDr3Source::try_from(row),
            Err(GaiaDr3Error::MissingRequiredField { field: "pmdec" })
        ));

        let mut row = raw_row();
        row.pmra_mas_per_yr = None;
        assert!(matches!(
            GaiaDr3Source::try_from(row),
            Err(GaiaDr3Error::MissingRequiredField { field: "pmra" })
        ));
    }

    #[test]
    fn gaia_proper_motion_rate_converts_to_degrees_per_julian_year() {
        let rate: GaiaProperMotionRate = GaiaProperMotionRate::new(3_600_000.0);
        let deg_per_jy: DegreesPerJulianYear = rate.to();
        assert!((deg_per_jy.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn finite_pmra_pmdec_build_domain_proper_motion() {
        let pm = ProperMotion::from_mu_alpha_star(
            GaiaProperMotionRate::new(1.5),
            GaiaProperMotionRate::new(-2.5),
        );
        assert_eq!(pm.ra_convention, RaProperMotionConvention::MuAlphaStar);
        let pmra_mas: GaiaProperMotionRate = pm.pm_ra.to();
        let pmdec_mas: GaiaProperMotionRate = pm.pm_dec.to();
        assert!((pmra_mas.value() - 1.5).abs() < 1e-12);
        assert!((pmdec_mas.value() + 2.5).abs() < 1e-12);
    }

    #[test]
    fn epoch_is_preserved() {
        let source = GaiaDr3Source::try_from(raw_row()).expect("typed source");

        assert_eq!(source.astrometry.epoch.value(), 2016.0);
    }

    #[test]
    fn missing_optional_photometry_stays_none() {
        let mut row = raw_row();
        row.phot_g_mean_mag = None;
        row.phot_bp_mean_mag = None;
        row.phot_rp_mean_mag = None;
        let source = GaiaDr3Source::try_from(row).expect("typed source");

        assert_eq!(source.photometry.phot_g_mean_mag, None);
        assert_eq!(source.photometry.phot_bp_mean_mag, None);
        assert_eq!(source.photometry.phot_rp_mean_mag, None);
    }

    #[test]
    fn gaia_xp_sampled_spectrum_validates_samples() {
        let spectrum = GaiaXpSampledSpectrum::new(vec![
            SpectralFluxSample::new(330.0, 1.0e-12).expect("sample"),
            SpectralFluxSample::new(650.0, 1.0e-12).expect("sample"),
        ])
        .expect("spectrum");

        assert_eq!(spectrum.samples.len(), 2);
    }

    #[test]
    fn photon_flux_integration_matches_linear_synthetic_spectrum() {
        let spectrum = GaiaXpSampledSpectrum::new(vec![
            SpectralFluxSample::new(400.0, 1.0e-12).expect("sample"),
            SpectralFluxSample::new(500.0, 1.0e-12).expect("sample"),
        ])
        .expect("spectrum");
        let band = SpectralBand::new("test", 400.0, 500.0).expect("band");

        let flux = integrate_photon_flux(&spectrum, band).expect("flux");
        let h = 6.626_070_15e-34;
        let c = 299_792_458.0;
        let expected =
            1.0e-12 * 1.0e9 / (h * c) * 0.5 * ((500.0e-9 * 500.0e-9) - (400.0e-9 * 400.0e-9));

        assert!((flux.photons_m2_s() - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn malformed_spectral_values_are_rejected() {
        assert!(SpectralFluxSample::new(f64::NAN, 1.0).is_err());
        assert!(SpectralFluxSample::new(400.0, f64::INFINITY).is_err());
        assert!(SpectralFluxSample::new(400.0, -1.0).is_err());

        let err = GaiaXpSampledSpectrum::new(vec![
            SpectralFluxSample::new(500.0, 1.0).expect("sample"),
            SpectralFluxSample::new(400.0, 1.0).expect("sample"),
        ])
        .expect_err("non-monotonic");
        assert_eq!(err, GaiaDr3Error::NonMonotonicSpectrum);
    }

    #[test]
    fn metadata_model_names_include_gaia_xp_production_candidate() {
        let model = StellarPhotometryModel::GaiaDr3XpPhotonRadiance330650NmV1;

        assert_eq!(format!("{model:?}"), "GaiaDr3XpPhotonRadiance330650NmV1");
    }

    #[test]
    fn parser_accepts_documented_columns() {
        let input = "\
source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,radial_velocity,phot_g_mean_mag,quality_ok,duplicated_source
7,10.0,20.0,2016.0,,,,,13.2,true,false";
        let rows = parse_gaia_dr3_csv_chunk(input).expect("rows");

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].source_id, 7);
        assert_eq!(rows[0].phot_g_mean_mag, Some(13.2));
        assert!(rows[0].quality.quality_ok);
        assert_eq!(rows[0].quality.duplicated_source, Some(false));
    }
}
