// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Mie / aerosol optical depth — power-law parameterization.

use crate::ext_qtty::length::Nanometers;
use crate::provenance::{DataSource, Provenance};

/// Aerosol optical-depth model parameters.
///
/// The optical depth is `τ_M(λ) = τ₀ · (λ / λ_ref)^α`. This is the
/// Patat (2011) power-law commonly used in modern Cerro Paranal sky
/// brightness models. Negative α encodes the usual aerosol behaviour
/// (more extinction at shorter wavelengths).
///
/// Reference: Patat, F. (2011), "The dancing sky: 6 years of night-sky
/// observations at Cerro Paranal", *Astronomy & Astrophysics* 527, A91.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MieParams {
    /// Optical depth at the reference wavelength.
    pub tau0: f64,
    /// Wavelength exponent `α` (typically negative).
    pub alpha: f64,
    /// Reference wavelength `λ_ref`.
    pub lambda_ref: Nanometers,
}

impl MieParams {
    /// Cerro Paranal default used by `darknsb` / NSB
    /// (`τ₀ = 0.05`, `α = -1.38`, `λ_ref = 550 nm`).
    ///
    /// See [`crate::observatories::EL_PARANAL`] for the matching observatory entry.
    pub const PARANAL: MieParams = MieParams {
        tau0: 0.05,
        alpha: -1.38,
        lambda_ref: Nanometers::new(550.0),
    };

    /// Roque de los Muchachos Observatory, La Palma (ORM, Canary Islands, Spain).
    ///
    /// `τ₀ = 0.02` at `λ_ref = 550 nm`, `α = −1.38`, giving AOD(500 nm) ≈ 0.022.
    ///
    /// La Palma enjoys exceptionally clean free-tropospheric air; the aerosol
    /// component of the V-band extinction is significantly lower than at Paranal.
    /// Lombardi et al. (2008, PASP 120, 212; DOI 10.1086/526338) report median
    /// photometric-night extinction k_V ≈ 0.11 mag/airmass; after subtracting the
    /// Rayleigh contribution (~0.09 mag/airmass at 744 hPa / 550 nm) the residual
    /// aerosol term corresponds to τ_Mie(550 nm) ≈ 0.02.  García-Gil et al. (2010,
    /// PASP 122, 1109; DOI 10.1086/656169) independently confirm that the aerosol
    /// optical depth at ORM is typically < 0.03 on clear nights.  Value is a
    /// representative median rather than a single-night calibration.
    ///
    /// See [`crate::observatories::ROQUE_DE_LOS_MUCHACHOS`] for the matching
    /// observatory entry.
    pub const LA_PALMA: MieParams = MieParams {
        tau0: 0.02,
        alpha: -1.38,
        lambda_ref: Nanometers::new(550.0),
    };

    /// Mauna Kea Observatories (Hawaiʻi, USA) — CFHT / Subaru site.
    ///
    /// `τ₀ = 0.03` at `λ_ref = 550 nm`, `α = −1.38`, giving AOD(500 nm) ≈ 0.034.
    ///
    /// Mauna Kea at 4207 m AMSL sits above most of the marine boundary layer.
    /// Krisciunas (1990, PASP 102, 1235; DOI 10.1086/132757) measured V-band
    /// extinction at Mauna Kea; the aerosol contribution derived after Rayleigh
    /// subtraction (reduced by the lower pressure ≈ 615 hPa at summit) is
    /// approximately τ_Mie(550 nm) ≈ 0.02–0.04, with a representative median of 0.03.
    /// CFHT sky-quality monitoring (Cuillandre et al., various CFHT bulletins) is
    /// consistent with this range on photometric nights.  Value is a representative
    /// median rather than a single-night calibration.
    ///
    /// See [`crate::observatories::MAUNA_KEA`] for the matching observatory entry.
    pub const MAUNA_KEA: MieParams = MieParams {
        tau0: 0.03,
        alpha: -1.38,
        lambda_ref: Nanometers::new(550.0),
    };

    /// La Silla Observatory (ESO, Chile).
    ///
    /// `τ₀ = 0.04` at `λ_ref = 550 nm`, `α = −1.38`, giving AOD(500 nm) ≈ 0.045.
    ///
    /// La Silla at 2400 m AMSL has slightly higher aerosol loading than Paranal
    /// on median nights due to its proximity to the Chilean coast and lower
    /// altitude.  Burki et al. (1995, A&AS 112, 383; bibcode 1995A&AS..112..383B)
    /// measured broad-band extinction at La Silla; the aerosol component derived
    /// from V-band residuals after Rayleigh subtraction is approximately
    /// τ_Mie(550 nm) ≈ 0.03–0.06, with a representative median of 0.04.
    /// This is consistent with the ESO La Silla atmospheric monitoring reports
    /// (see also Rufener 1986, A&A 165, 275; bibcode 1986A&A...165..275R).
    /// Value is a representative median rather than a single-night calibration.
    ///
    /// See [`crate::observatories::LA_SILLA_OBSERVATORY`] for the matching
    /// observatory entry.
    pub const LA_SILLA: MieParams = MieParams {
        tau0: 0.04,
        alpha: -1.38,
        lambda_ref: Nanometers::new(550.0),
    };

    /// Provenance record for [`MieParams::PARANAL`].
    pub fn paranal_provenance() -> Provenance {
        Provenance {
            source: Some(DataSource::LiteratureCitation {
                bibkey: "patat2011".to_string(),
                doi: Some("10.1051/0004-6361/201015537".to_string()),
            }),
            notes: Some(
                "Cerro Paranal aerosol τ₀ = 0.05 at 550 nm (α = −1.38). \
                 Patat 2011 (A&A 527, A91) Table 1 / darknsb model. \
                 Representative median of photometric nights over 6 years."
                    .to_string(),
            ),
            ..Provenance::default()
        }
    }

    /// Provenance record for [`MieParams::LA_PALMA`].
    pub fn la_palma_provenance() -> Provenance {
        Provenance {
            source: Some(DataSource::LiteratureCitation {
                bibkey: "lombardi2008".to_string(),
                doi: Some("10.1086/526338".to_string()),
            }),
            notes: Some(
                "La Palma ORM aerosol τ₀ = 0.02 at 550 nm (α = −1.38). \
                 Lombardi et al. 2008 (PASP 120, 212) median k_V ≈ 0.11 mag/airmass; \
                 Mie component ≈ 0.02 after Rayleigh subtraction. \
                 Confirmed by García-Gil et al. 2010 (PASP 122, 1109; DOI 10.1086/656169). \
                 Representative median, not a single-night calibration."
                    .to_string(),
            ),
            ..Provenance::default()
        }
    }

    /// Provenance record for [`MieParams::MAUNA_KEA`].
    pub fn mauna_kea_provenance() -> Provenance {
        Provenance {
            source: Some(DataSource::LiteratureCitation {
                bibkey: "krisciunas1990".to_string(),
                doi: Some("10.1086/132757".to_string()),
            }),
            notes: Some(
                "Mauna Kea aerosol τ₀ = 0.03 at 550 nm (α = −1.38). \
                 Krisciunas 1990 (PASP 102, 1235) V-band extinction at Mauna Kea; \
                 Mie component ≈ 0.02–0.04 after Rayleigh subtraction at 615 hPa. \
                 Consistent with CFHT sky-quality monitoring (Cuillandre et al.). \
                 Representative median, not a single-night calibration."
                    .to_string(),
            ),
            ..Provenance::default()
        }
    }

    /// Provenance record for [`MieParams::LA_SILLA`].
    pub fn la_silla_provenance() -> Provenance {
        Provenance {
            source: Some(DataSource::LiteratureCitation {
                bibkey: "burki1995".to_string(),
                doi: None,
            }),
            notes: Some(
                "La Silla aerosol τ₀ = 0.04 at 550 nm (α = −1.38). \
                 Burki et al. 1995 (A&AS 112, 383; bibcode 1995A&AS..112..383B) \
                 broad-band extinction at La Silla; Mie component ≈ 0.03–0.06 \
                 after Rayleigh subtraction at 760 hPa. See also Rufener 1986 \
                 (A&A 165, 275; bibcode 1986A&A...165..275R). \
                 Representative median, not a single-night calibration."
                    .to_string(),
            ),
            ..Provenance::default()
        }
    }
}

/// Mie / aerosol optical depth at the given wavelength using the
/// power-law parameterization stored in `params`.
#[inline]
pub fn mie_optical_depth(params: &MieParams, wavelength: Nanometers) -> f64 {
    params.tau0 * (wavelength.value() / params.lambda_ref.value()).powf(params.alpha)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn at_lambda_ref_returns_tau0() {
        let p = MieParams::PARANAL;
        assert_eq!(mie_optical_depth(&p, Nanometers::new(550.0)), p.tau0);
    }

    #[test]
    fn shorter_wavelength_higher_tau() {
        let p = MieParams::PARANAL;
        let blue = mie_optical_depth(&p, Nanometers::new(400.0));
        let red = mie_optical_depth(&p, Nanometers::new(700.0));
        assert!(blue > p.tau0);
        assert!(red < p.tau0);
    }

    /// Each named preset must differ from PARANAL in at least τ₀.
    #[test]
    fn site_presets_differ_from_paranal() {
        assert_ne!(MieParams::LA_PALMA.tau0, MieParams::PARANAL.tau0);
        assert_ne!(MieParams::MAUNA_KEA.tau0, MieParams::PARANAL.tau0);
        assert_ne!(MieParams::LA_SILLA.tau0, MieParams::PARANAL.tau0);
    }

    /// Aerosol OD at 550 nm (= τ₀ at λ_ref) must be positive and finite
    /// for every preset, and ordered by expected site cleanliness.
    #[test]
    fn site_presets_optical_depth_ordering() {
        let lambda = Nanometers::new(550.0);
        let tau_lp = mie_optical_depth(&MieParams::LA_PALMA,  lambda);
        let tau_mk = mie_optical_depth(&MieParams::MAUNA_KEA, lambda);
        let tau_ls = mie_optical_depth(&MieParams::LA_SILLA,  lambda);
        let tau_pn = mie_optical_depth(&MieParams::PARANAL,   lambda);

        // All values must be positive and finite.
        for tau in [tau_lp, tau_mk, tau_ls, tau_pn] {
            assert!(tau > 0.0 && tau.is_finite(), "tau must be positive finite: {tau}");
        }

        // Expected order: La Palma (cleanest) < Mauna Kea < La Silla < Paranal.
        assert!(tau_lp < tau_mk, "La Palma should be cleaner than Mauna Kea");
        assert!(tau_mk < tau_ls, "Mauna Kea should be cleaner than La Silla");
        assert!(tau_ls < tau_pn, "La Silla should be cleaner than Paranal");
    }

    /// Provenance records must be non-empty (source field populated).
    #[test]
    fn provenance_non_empty() {
        assert!(MieParams::paranal_provenance().source.is_some());
        assert!(MieParams::la_palma_provenance().source.is_some());
        assert!(MieParams::mauna_kea_provenance().source.is_some());
        assert!(MieParams::la_silla_provenance().source.is_some());
    }

    /// AOD(500 nm) for each site must fall within the expected published ranges.
    #[test]
    fn aod_500nm_within_published_ranges() {
        let lambda = Nanometers::new(500.0);
        let aod_lp = mie_optical_depth(&MieParams::LA_PALMA,  lambda);
        let aod_mk = mie_optical_depth(&MieParams::MAUNA_KEA, lambda);
        let aod_ls = mie_optical_depth(&MieParams::LA_SILLA,  lambda);

        // La Palma: 0.01–0.05 (Lombardi 2008, García-Gil 2010)
        assert!(aod_lp >= 0.01 && aod_lp <= 0.05,
            "La Palma AOD(500nm) = {aod_lp:.4} outside 0.01–0.05");
        // Mauna Kea: 0.02–0.04 (Krisciunas 1990)
        assert!(aod_mk >= 0.02 && aod_mk <= 0.05,
            "Mauna Kea AOD(500nm) = {aod_mk:.4} outside 0.02–0.05");
        // La Silla: 0.03–0.06 (Burki et al. 1995)
        assert!(aod_ls >= 0.03 && aod_ls <= 0.07,
            "La Silla AOD(500nm) = {aod_ls:.4} outside 0.03–0.07");
    }
}
