// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Validates the canonical observatory TOML and generates compatibility constants.

use serde::Deserialize;
use std::collections::HashSet;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Catalog {
    observatory: Vec<Record>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Record {
    name: String,
    longitude_deg: f64,
    latitude_deg: f64,
    height_m: f64,
    reference_pressure_hpa: f64,
    reference_temperature_k: Option<f64>,
    reference_relative_humidity: Option<f64>,
}

fn validate(record: &Record) -> Result<(), String> {
    if record.name.trim().is_empty() {
        return Err("name must not be empty".into());
    }
    validate_range("longitude_deg", record.longitude_deg, -180.0, 180.0)?;
    validate_range("latitude_deg", record.latitude_deg, -90.0, 90.0)?;
    validate_range("height_m", record.height_m, -500.0, 10_000.0)?;
    validate_range(
        "reference_pressure_hpa",
        record.reference_pressure_hpa,
        f64::MIN_POSITIVE,
        1_100.0,
    )?;
    if let Some(value) = record.reference_temperature_k {
        validate_range("reference_temperature_k", value, f64::MIN_POSITIVE, 400.0)?;
    }
    if let Some(value) = record.reference_relative_humidity {
        validate_range("reference_relative_humidity", value, 0.0, 1.0)?;
    }
    Ok(())
}

fn validate_range(field: &str, value: f64, min: f64, max: f64) -> Result<(), String> {
    if !value.is_finite() || value < min || value > max {
        return Err(format!(
            "{field} must be finite and in [{min}, {max}], got {value}"
        ));
    }
    Ok(())
}

fn symbol_for(name: &str) -> Option<&'static str> {
    match name {
        "El Paranal Observatory" => Some("EL_PARANAL"),
        "Roque de los Muchachos Observatory" => Some("ROQUE_DE_LOS_MUCHACHOS"),
        "Mauna Kea Observatory" => Some("MAUNA_KEA"),
        "La Silla Observatory" => Some("LA_SILLA_OBSERVATORY"),
        _ => None,
    }
}

fn option_quantity(value: Option<f64>, constructor: &str) -> String {
    value.map_or_else(
        || "None".into(),
        |value| format!("Some({constructor}::new({value:?}))"),
    )
}

fn observatory_expression(record: &Record) -> String {
    let temperature = option_quantity(record.reference_temperature_k, "Kelvins");
    let humidity = record
        .reference_relative_humidity
        .map_or_else(|| "None".into(), |value| format!("Some({value:?})"));
    format!(
        "Observatory {{\n        name: Cow::Borrowed({name:?}),\n        geodetic: Geodetic::new_raw(Degrees::new({lon:?}), Degrees::new({lat:?}), Meters::new({height:?})),\n        reference_pressure: Hectopascals::new({pressure:?}),\n        reference_temperature: {temperature},\n        reference_relative_humidity: {humidity},\n    }}",
        name = record.name,
        lon = record.longitude_deg,
        lat = record.latitude_deg,
        height = record.height_m,
        pressure = record.reference_pressure_hpa,
    )
}

fn main() {
    const SOURCE: &str = "data/observatories.toml";
    println!("cargo:rerun-if-changed={SOURCE}");

    let raw = fs::read_to_string(SOURCE).unwrap_or_else(|error| {
        panic!("failed to read canonical observatory catalog {SOURCE}: {error}")
    });
    let catalog: Catalog = toml::from_str(&raw)
        .unwrap_or_else(|error| panic!("invalid canonical observatory catalog {SOURCE}: {error}"));

    let mut names = HashSet::new();
    let mut compatibility_symbols = HashSet::new();
    let mut entries = Vec::new();
    let mut generated = String::from("// @generated from data/observatories.toml; do not edit.\n");
    for (index, record) in catalog.observatory.iter().enumerate() {
        validate(record).unwrap_or_else(|error| {
            panic!(
                "invalid observatory record {} (`{}`): {error}",
                index + 1,
                record.name
            )
        });
        assert!(
            names.insert(record.name.as_str()),
            "duplicate observatory name `{}`",
            record.name
        );
        let expression = observatory_expression(record);
        if let Some(symbol) = symbol_for(&record.name) {
            compatibility_symbols.insert(symbol);
            writeln!(
                generated,
                "/// The bundled {}, generated from `data/observatories.toml`.\npub const {symbol}: Observatory = {expression};",
                record.name
            )
            .expect("writing generated Rust to a String cannot fail");
            entries.push(symbol.to_owned());
        } else {
            entries.push(expression);
        }
    }

    for required in [
        "EL_PARANAL",
        "ROQUE_DE_LOS_MUCHACHOS",
        "MAUNA_KEA",
        "LA_SILLA_OBSERVATORY",
    ] {
        assert!(
            compatibility_symbols.contains(required),
            "canonical catalog is missing the observatory for compatibility symbol {required}"
        );
    }
    writeln!(
        generated,
        "fn generated_builtin_observatories() -> Vec<Observatory> {{\n    vec![\n        {}\n    ]\n}}",
        entries.join(",\n        ")
    )
    .expect("writing generated Rust to a String cannot fail");
    let output = PathBuf::from(env::var_os("OUT_DIR").expect("Cargo sets OUT_DIR"))
        .join("observatories_generated.rs");
    fs::write(&output, generated)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", output.display()));
}
