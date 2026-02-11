// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Segment lookup and coordinate evaluation for DE440 data.
//!
//! Given a Julian Date, locates the correct Chebyshev record and evaluates
//! the polynomial to produce position [x, y, z] in km (ICRF) and optionally
//! velocity [vx, vy, vz] in km/s.

use super::chebyshev;

/// J2000 epoch as Julian Date (TT/TDB).
const J2000_JD: f64 = 2_451_545.0;

/// Seconds per day.
const SECONDS_PER_DAY: f64 = 86_400.0;

/// Convert a Julian Date (TDB) to TDB seconds past J2000.
#[inline]
fn jd_to_et(jd: f64) -> f64 {
    (jd - J2000_JD) * SECONDS_PER_DAY
}

/// Evaluate position [x, y, z] in km (ICRF) for a single body at epoch `jd`.
///
/// # Arguments
/// - `jd`: Julian Date (TDB scale). The caller should apply TT→TDB if needed.
/// - `init`: segment initial epoch (TDB seconds past J2000)
/// - `intlen`: interval length (seconds)
/// - `ncoeff`: Chebyshev coefficients per coordinate
/// - `rsize`: doubles per record
/// - `n_records`: number of records
/// - `record_fn`: function to get the i-th record as a &[f64] slice
///
/// # Panics
/// Panics if `jd` is outside the segment's time span.
pub fn position(
    jd: f64,
    init: f64,
    intlen: f64,
    ncoeff: usize,
    _rsize: usize,
    n_records: usize,
    record_fn: fn(usize) -> &'static [f64],
) -> [f64; 3] {
    let et = jd_to_et(jd);

    // Find the record index
    let idx_f = (et - init) / intlen;
    let idx = (idx_f as usize).min(n_records - 1);

    let record = record_fn(idx);

    // First two values are MID and RADIUS
    let mid = record[0];
    let radius = record[1];

    // Normalize time to tau ∈ [-1, 1]
    let tau = (et - mid) / radius;

    // Coefficients for X, Y, Z start at offset 2
    let x_coeffs = &record[2..2 + ncoeff];
    let y_coeffs = &record[2 + ncoeff..2 + 2 * ncoeff];
    let z_coeffs = &record[2 + 2 * ncoeff..2 + 3 * ncoeff];

    [
        chebyshev::evaluate(x_coeffs, tau),
        chebyshev::evaluate(y_coeffs, tau),
        chebyshev::evaluate(z_coeffs, tau),
    ]
}

/// Evaluate velocity [vx, vy, vz] in km/day (ICRF) for a single body.
///
/// The derivative of the Chebyshev polynomial gives df/dtau. We multiply
/// by dtau/dt = 1/RADIUS to get df/dt in km/s, then by SECONDS_PER_DAY
/// to get km/day.
pub fn velocity(
    jd: f64,
    init: f64,
    intlen: f64,
    ncoeff: usize,
    _rsize: usize,
    n_records: usize,
    record_fn: fn(usize) -> &'static [f64],
) -> [f64; 3] {
    let et = jd_to_et(jd);

    let idx_f = (et - init) / intlen;
    let idx = (idx_f as usize).min(n_records - 1);

    let record = record_fn(idx);
    let mid = record[0];
    let radius = record[1];
    let tau = (et - mid) / radius;

    let x_coeffs = &record[2..2 + ncoeff];
    let y_coeffs = &record[2 + ncoeff..2 + 2 * ncoeff];
    let z_coeffs = &record[2 + 2 * ncoeff..2 + 3 * ncoeff];

    // df/dt = (df/dtau) * (dtau/dt) = (df/dtau) / radius
    // DE440 radius is in seconds, so df/dt is in km/s.
    // Convert to km/day by multiplying by SECONDS_PER_DAY.
    let scale = SECONDS_PER_DAY / radius;

    [
        chebyshev::evaluate_derivative(x_coeffs, tau) * scale,
        chebyshev::evaluate_derivative(y_coeffs, tau) * scale,
        chebyshev::evaluate_derivative(z_coeffs, tau) * scale,
    ]
}

/// Evaluate both position and velocity in one pass.
///
/// Returns `([x, y, z] km, [vx, vy, vz] km/day)`.
pub fn position_velocity(
    jd: f64,
    init: f64,
    intlen: f64,
    ncoeff: usize,
    _rsize: usize,
    n_records: usize,
    record_fn: fn(usize) -> &'static [f64],
) -> ([f64; 3], [f64; 3]) {
    let et = jd_to_et(jd);

    let idx_f = (et - init) / intlen;
    let idx = (idx_f as usize).min(n_records - 1);

    let record = record_fn(idx);
    let mid = record[0];
    let radius = record[1];
    let tau = (et - mid) / radius;

    let x_coeffs = &record[2..2 + ncoeff];
    let y_coeffs = &record[2 + ncoeff..2 + 2 * ncoeff];
    let z_coeffs = &record[2 + 2 * ncoeff..2 + 3 * ncoeff];

    let scale = SECONDS_PER_DAY / radius;

    let (px, vx) = chebyshev::evaluate_both(x_coeffs, tau);
    let (py, vy) = chebyshev::evaluate_both(y_coeffs, tau);
    let (pz, vz) = chebyshev::evaluate_both(z_coeffs, tau);

    (
        [px, py, pz],
        [vx * scale, vy * scale, vz * scale],
    )
}
