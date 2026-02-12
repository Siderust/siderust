// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE441 build-time pipeline — thin configuration wrapper.
//!
//! All logic lives in [`pipeline`](../pipeline.rs).
//! Only activated when the `de441` Cargo feature is enabled.

#[path = "../daf.rs"]
mod daf;
#[path = "../spk.rs"]
mod spk;
#[path = "../pipeline.rs"]
mod pipeline;

use std::path::Path;

const CONFIG: pipeline::DeConfig = pipeline::DeConfig {
    prefix: "de441",
    label: "DE441",
    bsp_url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp",
    bsp_filename: "de441_part-2.bsp",
    min_bsp_size: 1_500_000_000,
    download_timeout: 5400,
    size_hint: "~1.65 GB",
    repo_subdir: "scripts/de441/dataset",
};

/// Run the DE441 build pipeline.
pub fn run(data_dir: &Path) -> anyhow::Result<()> {
    pipeline::run(&CONFIG, data_dir)
}
