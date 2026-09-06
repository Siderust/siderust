// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use siderust::catalogs::{ObservatoryCatalog, ObservatoryCatalogError};

#[test]
fn public_catalog_composition_api_works_without_serde() {
    let mut catalog = ObservatoryCatalog::builtin();
    let baseline = catalog.clone();

    let error = catalog.extend(baseline).unwrap_err();
    assert!(matches!(
        error,
        ObservatoryCatalogError::DuplicateName { ref name }
            if name == "El Paranal Observatory"
    ));
}
