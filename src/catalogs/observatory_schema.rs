// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Private TOML schema and validation shared by build-time and runtime loading.

use serde::Deserialize;
use std::collections::HashSet;
use std::fmt;

/// Raw TOML observatory catalog representation.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct CatalogDto {
    /// Observatory records in source order.
    pub(super) observatory: Vec<ObservatoryDto>,
}

/// Raw TOML observatory record before conversion to typed Siderust quantities.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ObservatoryDto {
    /// Human-readable observatory name.
    pub(super) name: String,
    /// WGS84 geodetic longitude in degrees.
    pub(super) longitude_deg: f64,
    /// WGS84 geodetic latitude in degrees.
    pub(super) latitude_deg: f64,
    /// WGS84 ellipsoidal height in metres.
    pub(super) height_m: f64,
    /// Reference atmospheric pressure in hectopascals.
    pub(super) reference_pressure_hpa: f64,
    /// Optional reference temperature in kelvin.
    pub(super) reference_temperature_k: Option<f64>,
    /// Optional relative humidity as a fraction in `[0, 1]`.
    pub(super) reference_relative_humidity: Option<f64>,
}

/// Validation failure for one raw observatory field.
#[derive(Debug)]
pub(super) struct ValidationError {
    /// One-based record number.
    pub(super) record: usize,
    /// Observatory name, or a placeholder when the name itself is invalid.
    pub(super) name: String,
    /// Invalid TOML field.
    pub(super) field: &'static str,
    /// Human-readable validation requirement.
    pub(super) reason: String,
}

impl fmt::Display for ValidationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "invalid observatory record {} (`{}`), field `{}`: {}",
            self.record, self.name, self.field, self.reason
        )
    }
}

/// Parses the canonical `[[observatory]]` TOML representation.
pub(super) fn parse_catalog(input: &str) -> Result<CatalogDto, toml::de::Error> {
    toml::from_str(input)
}

/// Validates every observatory record and rejects duplicate exact names.
pub(super) fn validate_catalog(catalog: &CatalogDto) -> Result<(), ValidationError> {
    let mut names = HashSet::new();
    for (index, record) in catalog.observatory.iter().enumerate() {
        let record_number = index + 1;
        validate_record(record_number, record)?;
        if !names.insert(record.name.as_str()) {
            return Err(ValidationError {
                record: record_number,
                name: record.name.clone(),
                field: "name",
                reason: "duplicate observatory name".into(),
            });
        }
    }
    Ok(())
}

fn validate_record(record_number: usize, record: &ObservatoryDto) -> Result<(), ValidationError> {
    let display_name = if record.name.trim().is_empty() {
        "<empty>"
    } else {
        record.name.as_str()
    };
    if record.name.trim().is_empty() {
        return Err(ValidationError {
            record: record_number,
            name: display_name.to_owned(),
            field: "name",
            reason: "must not be empty".into(),
        });
    }

    validate_range(
        record_number,
        display_name,
        "longitude_deg",
        record.longitude_deg,
        -180.0,
        180.0,
    )?;
    validate_range(
        record_number,
        display_name,
        "latitude_deg",
        record.latitude_deg,
        -90.0,
        90.0,
    )?;
    validate_range(
        record_number,
        display_name,
        "height_m",
        record.height_m,
        -500.0,
        10_000.0,
    )?;
    validate_range(
        record_number,
        display_name,
        "reference_pressure_hpa",
        record.reference_pressure_hpa,
        f64::MIN_POSITIVE,
        1_100.0,
    )?;
    if let Some(value) = record.reference_temperature_k {
        validate_range(
            record_number,
            display_name,
            "reference_temperature_k",
            value,
            f64::MIN_POSITIVE,
            400.0,
        )?;
    }
    if let Some(value) = record.reference_relative_humidity {
        validate_range(
            record_number,
            display_name,
            "reference_relative_humidity",
            value,
            0.0,
            1.0,
        )?;
    }
    Ok(())
}

fn validate_range(
    record: usize,
    name: &str,
    field: &'static str,
    value: f64,
    min: f64,
    max: f64,
) -> Result<(), ValidationError> {
    if !value.is_finite() || value < min || value > max {
        return Err(ValidationError {
            record,
            name: name.to_owned(),
            field,
            reason: format!("must be finite and in [{min}, {max}], got {value}"),
        });
    }
    Ok(())
}
