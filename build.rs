// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Validates the canonical observatory TOML and generates compatibility constants.

#[path = "src/catalogs/observatory_schema.rs"]
mod observatory_schema;

use observatory_schema::{parse_catalog, validate_catalog, ObservatoryDto};
use std::collections::HashSet;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

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

fn observatory_expression(record: &ObservatoryDto) -> String {
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
    println!("cargo:rerun-if-changed=src/catalogs/observatory_schema.rs");

    let raw = fs::read_to_string(SOURCE).unwrap_or_else(|error| {
        panic!("failed to read canonical observatory catalog {SOURCE}: {error}")
    });
    let catalog = parse_catalog(&raw)
        .unwrap_or_else(|error| panic!("invalid canonical observatory catalog {SOURCE}: {error}"));
    validate_catalog(&catalog)
        .unwrap_or_else(|error| panic!("invalid canonical observatory catalog {SOURCE}: {error}"));

    let mut compatibility_symbols = HashSet::new();
    let mut entries = Vec::new();
    let mut generated = String::from("// @generated from data/observatories.toml; do not edit.\n");
    for record in &catalog.observatory {
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
