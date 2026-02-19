// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Embedded IERS EOP table from `finals2000A.all`, built at compile time.
//!
//! This module wraps the auto-generated `iers_eop_data.rs` and provides
//! a binary-search lookup with linear interpolation for any MJD within
//! the table's coverage.

#[allow(clippy::approx_constant)]
#[rustfmt::skip]
mod iers_eop_data {
    include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/generated/iers_eop_data.rs"));
}

pub use iers_eop_data::{EopEntry, EOP_TABLE};

/// Look up EOP values by MJD (UTC), interpolating linearly between
/// daily entries when the MJD falls between two table rows.
///
/// Returns `None` if `mjd` is outside the table range.
pub fn lookup(table: &[EopEntry], mjd: f64) -> Option<EopEntry> {
    if table.is_empty() {
        return None;
    }

    let first = table.first().unwrap().mjd;
    let last = table.last().unwrap().mjd;

    if mjd < first || mjd > last {
        return None;
    }

    // Binary search for the interval containing `mjd`.
    // The table has daily spacing (1.0 MJD), so we can also compute
    // a direct index guess for speed — but binary search is robust for
    // any spacing.
    let idx = table.partition_point(|e| e.mjd < mjd);

    if idx == 0 {
        return Some(table[0]);
    }
    if idx >= table.len() {
        return Some(*table.last().unwrap());
    }

    let lo = &table[idx - 1];
    let hi = &table[idx];

    // Exact match?
    if (hi.mjd - mjd).abs() < 1e-10 {
        return Some(*hi);
    }
    if (lo.mjd - mjd).abs() < 1e-10 {
        return Some(*lo);
    }

    // Linear interpolation
    let dt = hi.mjd - lo.mjd;
    let t = (mjd - lo.mjd) / dt;

    Some(EopEntry {
        mjd,
        xp: lo.xp + t * (hi.xp - lo.xp),
        yp: lo.yp + t * (hi.yp - lo.yp),
        dut1: lo.dut1 + t * (hi.dut1 - lo.dut1),
        dx: lo.dx + t * (hi.dx - lo.dx),
        dy: lo.dy + t * (hi.dy - lo.dy),
    })
}
