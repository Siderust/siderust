// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observatory catalog
//!
//! [`Observatory`] is Siderust's single scientific model for an observing
//! site. [`ObservatoryCatalog`] returns that same type whether records come
//! from Siderust's bundled catalog or a user-supplied TOML file.
//! Each position is a WGS84 geodetic longitude, latitude, and ellipsoidal
//! height; it is not a spherical or geocentric position. Pressure and optional
//! temperature use Siderust's typed quantities, while relative humidity is a
//! dimensionless fraction in `[0.0, 1.0]`.
//!
//! The bundled records originate in `data/observatories.toml`, which is
//! validated at build time and embedded into the crate. The legacy named
//! constants are generated from that file, so installed applications do not
//! need to locate catalog data on the filesystem.
//!
//! With the default `serde` feature, custom catalogs use the same direct
//! record representation:
//!
//! ```toml
//! [[observatory]]
//! name = "Example Observatory"
//! longitude_deg = 12.5
//! latitude_deg = 41.9
//! height_m = 800.0
//! reference_pressure_hpa = 920.0
//! reference_temperature_k = 283.0
//! reference_relative_humidity = 0.4
//! ```
//!
//! ```rust
//! # #[cfg(feature = "serde")] {
//! use siderust::catalogs::observatories::{Observatory, ObservatoryCatalog};
//!
//! let catalog = ObservatoryCatalog::from_toml(r#"
//!     [[observatory]]
//!     name = "Example Observatory"
//!     longitude_deg = 12.5
//!     latitude_deg = 41.9
//!     height_m = 800.0
//!     reference_pressure_hpa = 920.0
//! "#).unwrap();
//! let site: &Observatory = catalog.get("Example Observatory").unwrap();
//! assert_eq!(site.name, "Example Observatory");
//! # }
//! ```
//!
//! Downstream crates should consume this model and catalog rather than keep
//! duplicate coordinate registries. Domain-specific calibrations and
//! atmospheric models remain separate concerns.
//!
//! ## References
//!
//! - National Imagery and Mapping Agency (2000). *Department of Defense World
//!   Geodetic System 1984*, NIMA Technical Report TR8350.2, 3rd ed.
//! - Patat, F., Moehler, S., O'Brien, K., et al. (2011). "Optical atmospheric
//!   extinction over Cerro Paranal". *Astronomy & Astrophysics* **527**, A91.
//! - Instituto de Astrofísica de Canarias. *Site characterization of the
//!   Observatorio del Roque de los Muchachos*.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::{Degrees, Hectopascals, Kelvins, Meters};
use std::borrow::Cow;

#[cfg(feature = "serde")]
use std::path::{Path, PathBuf};

#[cfg(feature = "serde")]
#[path = "observatory_schema.rs"]
mod observatory_schema;

#[cfg(all(test, feature = "serde"))]
const BUNDLED_CATALOG_TOML: &str = include_str!("../../data/observatories.toml");

/// An astronomical observatory with a WGS84 location and reference atmosphere.
///
/// Bundled and runtime-loaded sites use this same type. The borrowed-or-owned
/// name preserves the existing compile-time constants while allowing custom
/// catalog records to own names loaded at runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct Observatory {
    /// Human-readable name.
    pub name: Cow<'static, str>,
    /// Geodetic position (WGS84 ellipsoidal longitude, latitude, and height).
    pub geodetic: Geodetic<ECEF>,
    /// Reference atmospheric pressure in hectopascals.
    pub reference_pressure: Hectopascals,
    /// Optional reference air temperature on the absolute kelvin scale.
    pub reference_temperature: Option<Kelvins>,
    /// Optional relative humidity as a dimensionless fraction in `[0.0, 1.0]`.
    pub reference_relative_humidity: Option<f64>,
}

impl Observatory {
    /// Returns the geodetic position of this observatory.
    #[inline]
    pub const fn geodetic(&self) -> Geodetic<ECEF> {
        self.geodetic
    }

    #[cfg(feature = "serde")]
    fn from_dto(record: observatory_schema::ObservatoryDto) -> Self {
        Self {
            name: Cow::Owned(record.name),
            geodetic: Geodetic::new_raw(
                Degrees::new(record.longitude_deg),
                Degrees::new(record.latitude_deg),
                Meters::new(record.height_m),
            ),
            reference_pressure: Hectopascals::new(record.reference_pressure_hpa),
            reference_temperature: record.reference_temperature_k.map(Kelvins::new),
            reference_relative_humidity: record.reference_relative_humidity,
        }
    }
}

include!(concat!(env!("OUT_DIR"), "/observatories_generated.rs"));

/// An ordered collection of [`Observatory`] records.
///
/// Iteration preserves TOML record order. Name lookup is exact and
/// case-sensitive; no stable-ID or alias semantics are imposed.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObservatoryCatalog {
    observatories: Vec<Observatory>,
}

impl ObservatoryCatalog {
    /// Returns the catalog bundled with Siderust.
    ///
    /// Its records are generated from and validated against the canonical
    /// `data/observatories.toml` file during the build.
    pub fn builtin() -> Self {
        Self {
            observatories: generated_builtin_observatories(),
        }
    }

    /// Returns the number of observatories in the catalog.
    pub fn len(&self) -> usize {
        self.observatories.len()
    }

    /// Returns whether the catalog contains no observatories.
    pub fn is_empty(&self) -> bool {
        self.observatories.is_empty()
    }

    /// Iterates over observatories in their TOML record order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &Observatory> {
        self.observatories.iter()
    }

    /// Returns the observatories as an ordered slice.
    pub fn as_slice(&self) -> &[Observatory] {
        &self.observatories
    }

    /// Finds an observatory by exact, case-sensitive name.
    pub fn get(&self, name: &str) -> Option<&Observatory> {
        self.observatories
            .iter()
            .find(|observatory| observatory.name == name)
    }

    /// Parses and validates a catalog using `[[observatory]]` TOML records.
    ///
    /// This API is available with the `serde` feature, which is enabled by
    /// default.
    #[cfg(feature = "serde")]
    pub fn from_toml(input: &str) -> Result<Self, ObservatoryCatalogError> {
        let catalog = observatory_schema::parse_catalog(input)?;
        observatory_schema::validate_catalog(&catalog)
            .map_err(ObservatoryCatalogError::from_validation)?;
        Ok(Self {
            observatories: catalog
                .observatory
                .into_iter()
                .map(Observatory::from_dto)
                .collect(),
        })
    }

    /// Reads, parses, and validates an external TOML catalog.
    ///
    /// This API is available with the `serde` feature, which is enabled by
    /// default.
    #[cfg(feature = "serde")]
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, ObservatoryCatalogError> {
        let path = path.as_ref();
        let input =
            std::fs::read_to_string(path).map_err(|source| ObservatoryCatalogError::Read {
                path: path.to_owned(),
                source,
            })?;
        Self::from_toml(&input)
    }
}

impl<'a> IntoIterator for &'a ObservatoryCatalog {
    type Item = &'a Observatory;
    type IntoIter = std::slice::Iter<'a, Observatory>;

    fn into_iter(self) -> Self::IntoIter {
        self.observatories.iter()
    }
}

/// Errors produced while reading, parsing, or validating an observatory catalog.
#[cfg(feature = "serde")]
#[derive(Debug, thiserror::Error)]
pub enum ObservatoryCatalogError {
    /// The external catalog file could not be read.
    #[error("failed to read observatory catalog `{path}`: {source}")]
    Read {
        /// Path that could not be read.
        path: PathBuf,
        /// Underlying filesystem error.
        #[source]
        source: std::io::Error,
    },
    /// The input is not a valid catalog in the canonical TOML representation.
    #[error("invalid observatory catalog TOML: {0}")]
    Toml(#[from] toml::de::Error),
    /// A parsed record contains a scientifically invalid value.
    #[error("invalid observatory record {record} (`{name}`), field `{field}`: {reason}")]
    InvalidField {
        /// One-based record number.
        record: usize,
        /// Observatory name, or a placeholder when the name itself is invalid.
        name: String,
        /// Invalid TOML field.
        field: &'static str,
        /// Human-readable validation requirement.
        reason: String,
    },
}

#[cfg(feature = "serde")]
impl ObservatoryCatalogError {
    fn from_validation(error: observatory_schema::ValidationError) -> Self {
        Self::InvalidField {
            record: error.record,
            name: error.name,
            field: error.field,
            reason: error.reason,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "serde")]
    use std::error::Error as _;

    #[test]
    fn builtin_catalog_is_ordered_and_uses_public_model() {
        let catalog = ObservatoryCatalog::builtin();
        let sites: &[Observatory] = catalog.as_slice();
        assert!(sites.len() >= 4);
        for expected in [
            EL_PARANAL,
            ROQUE_DE_LOS_MUCHACHOS,
            MAUNA_KEA,
            LA_SILLA_OBSERVATORY,
        ] {
            assert_eq!(catalog.get(expected.name.as_ref()), Some(&expected));
        }
    }

    #[cfg(feature = "serde")]
    fn one_record(overrides: &str) -> String {
        format!(
            "[[observatory]]\nname = \"Test Site\"\nlongitude_deg = 1.0\nlatitude_deg = 2.0\nheight_m = 3.0\nreference_pressure_hpa = 900.0\n{overrides}"
        )
    }

    #[cfg(feature = "serde")]
    #[test]
    fn canonical_toml_matches_generated_compatibility_constants() {
        let parsed = ObservatoryCatalog::from_toml(BUNDLED_CATALOG_TOML).unwrap();
        assert_eq!(parsed, ObservatoryCatalog::builtin());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn custom_observatory_and_optional_atmosphere_load() {
        let catalog = ObservatoryCatalog::from_toml(&one_record(
            "reference_temperature_k = 280.5\nreference_relative_humidity = 0.42\n",
        ))
        .unwrap();
        let site: &Observatory = catalog.get("Test Site").unwrap();
        assert_eq!(site.geodetic.lon.value(), 1.0);
        assert_eq!(site.geodetic.lat.value(), 2.0);
        assert_eq!(site.geodetic.height.value(), 3.0);
        assert_eq!(site.reference_pressure.value(), 900.0);
        assert_eq!(site.reference_temperature.unwrap().value(), 280.5);
        assert_eq!(site.reference_relative_humidity, Some(0.42));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn optional_atmosphere_can_be_omitted() {
        let catalog = ObservatoryCatalog::from_toml(&one_record("")).unwrap();
        let site = catalog.get("Test Site").unwrap();
        assert!(site.reference_temperature.is_none());
        assert!(site.reference_relative_humidity.is_none());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn rejects_invalid_scientific_values() {
        for (field, valid, invalid) in [
            ("longitude_deg", "1.0", "181.0"),
            ("latitude_deg", "2.0", "-91.0"),
            ("height_m", "3.0", "10001.0"),
            ("reference_pressure_hpa", "900.0", "0.0"),
        ] {
            let input = one_record("").replacen(
                &format!("{field} = {valid}"),
                &format!("{field} = {invalid}"),
                1,
            );
            let error = ObservatoryCatalog::from_toml(&input).unwrap_err();
            assert!(error.to_string().contains(field), "{error}");
        }
        for extra in [
            "reference_temperature_k = 0.0\n",
            "reference_relative_humidity = 1.1\n",
        ] {
            let error = ObservatoryCatalog::from_toml(&one_record(extra)).unwrap_err();
            assert!(error.to_string().contains("reference_"), "{error}");
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn rejects_duplicate_names() {
        let input = format!("{}{}", one_record(""), one_record(""));
        let error = ObservatoryCatalog::from_toml(&input).unwrap_err();
        assert!(error.to_string().contains("duplicate observatory name"));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn malformed_toml_is_useful() {
        let error = ObservatoryCatalog::from_toml("[[observatory]\n").unwrap_err();
        assert!(error.to_string().contains("TOML"));
        assert!(error.source().is_some());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn path_errors_include_the_path() {
        let path = Path::new("definitely-missing-observatory-catalog.toml");
        let error = ObservatoryCatalog::from_path(path).unwrap_err();
        assert!(error.to_string().contains(path.to_str().unwrap()));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn external_path_loads_unknown_observatory() {
        let path = std::env::temp_dir().join(format!(
            "siderust-observatories-{}-{}.toml",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        std::fs::write(&path, one_record("")).unwrap();
        let catalog = ObservatoryCatalog::from_path(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(catalog.get("Test Site").unwrap().name, "Test Site");
    }
}
