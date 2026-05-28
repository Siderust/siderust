// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded Sun-Earth Lagrange Chebyshev kernels loaded from the archive.
//!
//! The five SCK binary files in `archive/lagrange/vsop87/` are embedded at
//! compile time via `include_bytes!`.  They are lazily parsed on first access
//! so that the parsing cost is paid only when the Lagrange API is first called.
//!
//! # Archive requirement
//!
//! The `archive/` git submodule must be checked out before compiling with the
//! `lagrange-centers` feature:
//!
//! ```text
//! git submodule update --init --recursive
//! ```

use std::sync::LazyLock;

use crate::formats::sck::parse_sck;

pub(crate) const NCOEFF: usize = 8;

// Each `.sck` file is ~475 KB (2283 records × 26 f64s × 8 bytes + 64-byte header).
static L1_BYTES: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/archive/lagrange/vsop87/l1.sck"
));
static L2_BYTES: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/archive/lagrange/vsop87/l2.sck"
));
static L3_BYTES: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/archive/lagrange/vsop87/l3.sck"
));
static L4_BYTES: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/archive/lagrange/vsop87/l4.sck"
));
static L5_BYTES: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/archive/lagrange/vsop87/l5.sck"
));

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
