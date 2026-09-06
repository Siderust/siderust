// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression coverage for the public observatory catalog composition API
//! when Siderust is built without default Cargo features.

use siderust::catalogs::{ObservatoryCatalog, ObservatoryCatalogError};

#[test]
fn public_catalog_composition_api_works_without_serde() {
    let mut catalog = ObservatoryCatalog::default();
    catalog.extend(ObservatoryCatalog::builtin()).unwrap();
    assert!(!catalog.is_empty());
    assert!(catalog.get("El Paranal Observatory").is_some());

    let error = catalog.extend(ObservatoryCatalog::builtin()).unwrap_err();
    assert!(matches!(
        error,
        ObservatoryCatalogError::DuplicateName { ref name }
            if name == "El Paranal Observatory"
    ));
}
