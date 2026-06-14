// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded Sun-Earth Lagrange Chebyshev kernels.
//!
//! The five SCK binary files (`l1.sck`..`l5.sck`) are embedded in
//! [`siderust_archive`] under `src/lagrange/vsop87/` and activated by the
//! `lagrange` feature of that crate (pulled in automatically when this crate
//! is built with `lagrange-centers`).
//!
//! The bytes are lazily parsed on first access so the parsing cost is paid only
//! when the Lagrange API is first called.

use std::sync::LazyLock;

use crate::archive::lagrange::vsop87 as sck_data;
use crate::formats::sck::parse_sck;

pub(crate) const NCOEFF: usize = 8;

static L1: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(sck_data::L1_BYTES)
        .expect("embedded L1 SCK must be valid")
        .records
});
static L2: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(sck_data::L2_BYTES)
        .expect("embedded L2 SCK must be valid")
        .records
});
static L3: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(sck_data::L3_BYTES)
        .expect("embedded L3 SCK must be valid")
        .records
});
static L4: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(sck_data::L4_BYTES)
        .expect("embedded L4 SCK must be valid")
        .records
});
static L5: LazyLock<Vec<f64>> = LazyLock::new(|| {
    parse_sck(sck_data::L5_BYTES)
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
