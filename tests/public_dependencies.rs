//! Compile-time checks for crates exposed by Siderust's public API.

use siderust::{affn, chrono, optica, principia, qtty, tempoch};

#[test]
fn public_api_dependencies_are_reexported() {
    let _: qtty::Meters = qtty::Meters::new(1.0);
    let _: Option<affn::Rotation3> = None;
    let _: chrono::DateTime<chrono::Utc> = chrono::DateTime::UNIX_EPOCH;
    let _: Option<principia::IntegratorTolerances> = None;
    let _: Option<tempoch::JulianDate<tempoch::TT>> = None;

    // Merely naming an optica API type verifies that the crate is reachable;
    // constructing a spectrum would add irrelevant fixture data to this test.
    let _: Option<optica::grid::OutOfRange> = None;
}
