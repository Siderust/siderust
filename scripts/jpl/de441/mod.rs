// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE441 build-time pipeline, thin configuration wrapper.
//!
//! All logic lives in the shared JPL build pipeline.
//! Only activated when the `de441` Cargo feature is enabled.

use crate::jpl_pipeline as pipeline;
use std::path::Path;

const CONFIG: pipeline::DeConfig = pipeline::DeConfig {
    prefix: "de441",
    label: "DE441",
    bsp_url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
    bsp_filename: "de441_part-2.bsp",
    min_bsp_size: 1_500_000_000,
    download_timeout: 5400,
    size_hint: "~1.65 GB",
};

/// Run the DE441 build pipeline.
pub(crate) fn run(data_dir: &Path) -> anyhow::Result<()> {
    pipeline::run(&CONFIG, data_dir)
}
