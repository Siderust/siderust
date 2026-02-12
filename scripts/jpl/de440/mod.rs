// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! DE440 build-time pipeline — thin configuration wrapper.
//!
//! All logic lives in [`pipeline`](../pipeline.rs).
//! Only activated when the `de440` Cargo feature is enabled.

#[path = "../daf.rs"]
mod daf;
#[path = "../spk.rs"]
mod spk;
#[path = "../pipeline.rs"]
mod pipeline;

use std::path::Path;

const CONFIG: pipeline::DeConfig = pipeline::DeConfig {
    prefix: "de440",
    label: "DE440",
    bsp_url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    bsp_filename: "de440.bsp",
    min_bsp_size: 100_000_000,
    download_timeout: 900,
    size_hint: "~120 MB",
    repo_subdir: "scripts/de440/dataset",
};

/// Run the DE440 build pipeline.
pub fn run(data_dir: &Path) -> anyhow::Result<()> {
    pipeline::run(&CONFIG, data_dir)
}
