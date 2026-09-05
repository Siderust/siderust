// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use siderust::catalogs::observatories::{
    ObservatoryCatalog, CTAO_NORTH, CTAO_SOUTH, EL_PARANAL, ROQUE_DE_LOS_MUCHACHOS,
};

#[test]
fn ctao_generated_constants_match_builtin_catalog() {
    let catalog = ObservatoryCatalog::builtin();

    assert_eq!(catalog.get("CTAO North"), Some(&CTAO_NORTH));
    assert_eq!(catalog.get("CTAO South"), Some(&CTAO_SOUTH));
}

#[test]
fn ctao_records_match_authoritative_array_coordinates() {
    assert_eq!(CTAO_NORTH.geodetic.lon.value(), -17.892005);
    assert_eq!(CTAO_NORTH.geodetic.lat.value(), 28.762164);
    assert_eq!(CTAO_NORTH.geodetic.height.value(), 2240.2);
    assert_eq!(CTAO_NORTH.reference_pressure.value(), 775.0);

    assert_eq!(CTAO_SOUTH.geodetic.lon.value(), -70.31634444444444);
    assert_eq!(CTAO_SOUTH.geodetic.lat.value(), -24.683427777777776);
    assert_eq!(CTAO_SOUTH.geodetic.height.value(), 2184.6);
    assert_eq!(CTAO_SOUTH.reference_pressure.value(), 780.0);
}

#[test]
fn ctao_sites_remain_distinct_from_nearby_observatories() {
    assert_ne!(CTAO_NORTH.geodetic, ROQUE_DE_LOS_MUCHACHOS.geodetic);
    assert_ne!(CTAO_SOUTH.geodetic, EL_PARANAL.geodetic);
}
