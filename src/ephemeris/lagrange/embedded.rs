// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded Sun-Earth Lagrange Chebyshev kernels loaded from the archive.
//!
//! The five SCK binary files (`l1.sck`..`l5.sck`) live in the
//! [Siderust Archive](https://github.com/Siderust/archive) under
//! `lagrange/vsop87/` and are embedded at compile time via `include_bytes!`.
//! They are lazily parsed on first access so the parsing cost is paid only
//! when the Lagrange API is first called.
//!
//! # Archive requirement
//!
//! Building with the `lagrange-centers` feature requires a checkout of the
//! archive repository. `build.rs` resolves its location, in order:
//!
//! 1. `SIDERUST_ARCHIVE_ROOT` environment variable, if it points at a
//!    directory containing `MANIFEST.toml`.
//! 2. `./archive/` next to this crate.
//! 3. `../archive/` next to this crate.
//!
//! If none of the above contains `lagrange/vsop87/l[1-5].sck`, the build
//! fails with an actionable `compile_error!` rather than producing an opaque
//! `include_bytes!` not-found error.

use std::sync::LazyLock;

use crate::formats::sck::parse_sck;

pub(crate) const NCOEFF: usize = 8;

// Each `.sck` file is ~475 KB (2283 records × 26 f64s × 8 bytes + 64-byte header).
// L1_BYTES..L5_BYTES are declared by build.rs in OUT_DIR/lagrange_paths.rs.
include!(concat!(env!("OUT_DIR"), "/lagrange_paths.rs"));

static L1: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(L1_BYTES)
        .expect("embedded L1 SCK must be valid")
        .records
});
static L2: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(L2_BYTES)
        .expect("embedded L2 SCK must be valid")
        .records
});
static L3: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(L3_BYTES)
        .expect("embedded L3 SCK must be valid")
        .records
});
static L4: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(L4_BYTES)
        .expect("embedded L4 SCK must be valid")
        .records
});
static L5: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(L5_BYTES)
        .expect("embedded L5 SCK must be valid")
        .records
});

/// SHA-256 of the original concatenated f64 payload from which the SCK files were derived.
pub(crate) const CHECKSUM: &str =
    "69c8eddcaeedc6f2e906a2ac5a4a2c2c7d97f491c08ce2b5ac6f4038e1966846";

#[must_use]
pub(crate) fn records_l1() -> &'static [f64] {
    (*L1).as_slice()
}

#[must_use]
pub(crate) fn records_l2() -> &'static [f64] {
    (*L2).as_slice()
}

#[must_use]
pub(crate) fn records_l3() -> &'static [f64] {
    (*L3).as_slice()
}

#[must_use]
pub(crate) fn records_l4() -> &'static [f64] {
    (*L4).as_slice()
}

#[must_use]
pub(crate) fn records_l5() -> &'static [f64] {
    (*L5).as_slice()
}
