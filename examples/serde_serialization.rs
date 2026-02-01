// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Serialization and Deserialization Example
//!
//! This example demonstrates how to serialize and deserialize siderust types
//! using serde with JSON format. It covers:
//! - Julian dates and time types
//! - Cartesian coordinates (positions and directions)
//! - Spherical coordinates
//! - Round-trip serialization (serialize → deserialize → verify)
//! - Saving/loading astronomical data to/from files

use qtty::*;
use serde::{Deserialize, Serialize};
use serde_json;
use siderust::coordinates::{cartesian, frames, spherical};
use siderust::time::JulianDate;
use std::fs;

#[derive(Debug, Serialize, Deserialize)]
struct ObservationData {
    /// Julian date of the observation
    julian_date: JulianDate,
    /// Observer position in geocentric ICRS coordinates
    observer_position:
        cartesian::Position<siderust::coordinates::centers::Geocentric, frames::ICRS, Kilometer>,
    /// Target direction in ICRS frame
    target_direction: cartesian::Direction<frames::ICRS>,
    /// Target spherical coordinates (RA/Dec)
    target_spherical:
        spherical::Position<siderust::coordinates::centers::Geocentric, frames::ICRS, Kilometer>,
}

fn main() {
    println!("=== Siderust Serialization/Deserialization Example ===\n");

    // =========================================================================
    // 1. Serializing Julian Dates
    // =========================================================================
    println!("1. JULIAN DATE SERIALIZATION");
    println!("----------------------------");

    let j2000 = JulianDate::J2000;
    let json = serde_json::to_string(&j2000).expect("Failed to serialize JulianDate");
    println!("JulianDate J2000 serialized to JSON: {}", json);

    let recovered: JulianDate =
        serde_json::from_str(&json).expect("Failed to deserialize JulianDate");
    println!("Deserialized back: JD = {:.6}", recovered.value());
    println!(
        "Round-trip successful: {}\n",
        (j2000.value() - recovered.value()).abs() < 1e-12
    );

    // =========================================================================
    // 2. Serializing Cartesian Positions
    // =========================================================================
    println!("2. CARTESIAN POSITION SERIALIZATION");
    println!("------------------------------------");

    let position = cartesian::Position::<
        siderust::coordinates::centers::Geocentric,
        frames::ICRS,
        Kilometer,
    >::new(6_700.0, 0.0, 0.0);

    let json =
        serde_json::to_string_pretty(&position).expect("Failed to serialize cartesian position");
    println!("Cartesian position serialized to JSON:");
    println!("{}", json);

    let recovered: cartesian::Position<
        siderust::coordinates::centers::Geocentric,
        frames::ICRS,
        Kilometer,
    > = serde_json::from_str(&json).expect("Failed to deserialize cartesian position");
    println!("\nDeserialized position:");
    println!("  X = {:.3} km", recovered.x().value());
    println!("  Y = {:.3} km", recovered.y().value());
    println!("  Z = {:.3} km", recovered.z().value());
    println!();

    // =========================================================================
    // 3. Serializing Cartesian Directions
    // =========================================================================
    println!("3. CARTESIAN DIRECTION SERIALIZATION");
    println!("------------------------------------");

    let direction = cartesian::Direction::<frames::ICRS>::new(1.0, 0.0, 0.0);

    let json = serde_json::to_string_pretty(&direction).expect("Failed to serialize direction");
    println!("Direction vector serialized to JSON:");
    println!("{}", json);

    let recovered: cartesian::Direction<frames::ICRS> =
        serde_json::from_str(&json).expect("Failed to deserialize direction");
    println!("\nDeserialized direction (normalized):");
    println!("  x̂ = {:.6}", recovered.x());
    println!("  ŷ = {:.6}", recovered.y());
    println!("  ẑ = {:.6}", recovered.z());
    println!();

    // =========================================================================
    // 4. Serializing Spherical Coordinates
    // =========================================================================
    println!("4. SPHERICAL COORDINATE SERIALIZATION");
    println!("-------------------------------------");

    // Create a spherical position (RA=0°, Dec=+45°, r=1 AU)
    let spherical_pos = spherical::Position::<
        siderust::coordinates::centers::Heliocentric,
        frames::ICRS,
        AstronomicalUnit,
    >::new(
        0.0 * DEG,  // RA
        45.0 * DEG, // Dec
        1.0 * AU,   // distance
    );

    let json = serde_json::to_string_pretty(&spherical_pos).expect("Failed to serialize spherical");
    println!("ICRS spherical position serialized to JSON (uses 'ra'/'dec'):");
    println!("{}", json);

    let recovered: spherical::Position<
        siderust::coordinates::centers::Heliocentric,
        frames::ICRS,
        AstronomicalUnit,
    > = serde_json::from_str(&json).expect("Failed to deserialize spherical");
    println!("\nDeserialized spherical coordinates:");
    println!("  RA  = {:.3}°", recovered.ra().value());
    println!("  Dec = {:.3}°", recovered.dec().value());
    println!("  Distance = {:.3} AU", recovered.distance().value());
    println!();

    // =========================================================================
    // 4b. Frame-Specific Field Names
    // =========================================================================
    println!("4b. FRAME-SPECIFIC FIELD NAMES");
    println!("------------------------------");
    println!("Different coordinate frames use domain-appropriate JSON field names:\n");

    // Ecliptic coordinates use lon/lat
    let ecliptic_pos = spherical::Position::<
        siderust::coordinates::centers::Heliocentric,
        frames::Ecliptic,
        AstronomicalUnit,
    >::new(
        120.0 * DEG, // Ecliptic longitude
        30.0 * DEG,  // Ecliptic latitude
        2.5 * AU,
    );
    let json = serde_json::to_string_pretty(&ecliptic_pos).unwrap();
    println!("Ecliptic coordinates (uses 'lon'/'lat'):");
    println!("{}\n", json);

    // Horizontal coordinates use az/alt
    let horizontal_dir = spherical::Direction::<frames::Horizontal>::new(
        45.0 * DEG,  // Altitude
        180.0 * DEG, // Azimuth (South)
    );
    let json = serde_json::to_string_pretty(&horizontal_dir).unwrap();
    println!("Horizontal direction (uses 'az'/'alt'):");
    println!("{}\n", json);

    // Galactic coordinates use l/b
    let galactic_pos = spherical::Position::<
        siderust::coordinates::centers::Barycentric,
        frames::Galactic,
        Parsec,
    >::new_raw(
        0.0 * DEG,   // Galactic latitude (b) - in plane
        0.0 * DEG,   // Galactic longitude (l) - toward center
        8000.0 * PC, // Distance to galactic center
    );
    let json = serde_json::to_string_pretty(&galactic_pos).unwrap();
    println!("Galactic coordinates (uses 'l'/'b'):");
    println!("{}\n", json);

    // Demonstrate round-trip with ecliptic
    let recovered: spherical::Position<
        siderust::coordinates::centers::Heliocentric,
        frames::Ecliptic,
        AstronomicalUnit,
    > = serde_json::from_str(&serde_json::to_string(&ecliptic_pos).unwrap()).unwrap();
    println!("Ecliptic round-trip successful:");
    println!(
        "  lon = {:.1}°, lat = {:.1}°",
        recovered.lon().value(),
        recovered.lat().value()
    );
    println!();

    // =========================================================================
    // 5. Serializing Complex Observation Data
    // =========================================================================
    println!("5. COMPLEX OBSERVATION DATA");
    println!("---------------------------");

    let observation = ObservationData {
        julian_date: JulianDate::new(2_459_580.5), // 2022-01-01 00:00:00 TT
        observer_position: cartesian::Position::new(
            6_371.0, // Earth radius at equator
            0.0, 0.0,
        ),
        target_direction: cartesian::Direction::new(
            0.707, // Pointing roughly northeast
            0.707, 0.0,
        ),
        target_spherical: spherical::Position::<
            siderust::coordinates::centers::Geocentric,
            frames::ICRS,
            Kilometer,
        >::new(
            45.0 * DEG,     // RA = 45°
            30.0 * DEG,     // Dec = 30°
            384_400.0 * KM, // Moon distance
        ),
    };

    let json = serde_json::to_string_pretty(&observation).expect("Failed to serialize observation");
    println!("Complete observation data serialized:");
    println!("{}", json);

    let recovered: ObservationData =
        serde_json::from_str(&json).expect("Failed to deserialize observation");
    println!("\nDeserialized observation:");
    println!("  JD = {:.6}", recovered.julian_date.value());
    println!(
        "  Observer at ({:.1}, {:.1}, {:.1}) km",
        recovered.observer_position.x().value(),
        recovered.observer_position.y().value(),
        recovered.observer_position.z().value()
    );
    println!(
        "  Target RA  = {:.1}°",
        recovered.target_spherical.ra().value()
    );
    println!(
        "  Target Dec = {:.1}°",
        recovered.target_spherical.dec().value()
    );
    println!();

    // =========================================================================
    // 6. Saving to File and Loading Back
    // =========================================================================
    println!("6. FILE I/O");
    println!("-----------");

    let filename = "/tmp/siderust_observation.json";

    // Save to file
    let json = serde_json::to_string_pretty(&observation).expect("Failed to serialize for file");
    fs::write(filename, json).expect("Failed to write file");
    println!("Observation data saved to: {}", filename);

    // Load from file
    let file_content = fs::read_to_string(filename).expect("Failed to read file");
    let loaded: ObservationData =
        serde_json::from_str(&file_content).expect("Failed to deserialize from file");
    println!("Data loaded from file:");
    println!("  JD = {:.6}", loaded.julian_date.value());
    println!(
        "  Match: {}",
        (loaded.julian_date.value() - observation.julian_date.value()).abs() < 1e-12
    );
    println!();

    // =========================================================================
    // 7. Working with Collections
    // =========================================================================
    println!("7. SERIALIZING COLLECTIONS");
    println!("--------------------------");

    let observations = vec![
        JulianDate::new(2_459_580.5),
        JulianDate::new(2_459_581.5),
        JulianDate::new(2_459_582.5),
    ];

    let json = serde_json::to_string_pretty(&observations).expect("Failed to serialize vector");
    println!("Multiple Julian dates:");
    println!("{}", json);

    let recovered: Vec<JulianDate> =
        serde_json::from_str(&json).expect("Failed to deserialize vector");
    println!("\nDeserialized {} Julian dates:", recovered.len());
    for (i, jd) in recovered.iter().enumerate() {
        println!("  [{}] JD = {:.1}", i, jd.value());
    }
    println!();

    println!("=== Example Complete ===");
}
