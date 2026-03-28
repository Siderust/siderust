// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! **Qtty-aware FFI entry points**, accept `QttyQuantity` for typed angles and thresholds.
//!
//! These functions accept [`qtty_ffi::QttyQuantity`] values instead of raw `f64` degrees,
//! providing compile-time dimension checking on the FFI boundary when the caller's
//! language supports it.
//!
//! # Example (C++)
//!
//! ```cpp
//! // Using qtty-cpp's Quantity type
//! auto threshold = qtty::Degrees(10.0);
//! auto status = siderust_above_threshold_qty(
//!     subject, observer, window,
//!     threshold.to_ffi(),  // Convert to qtty_quantity_t
//!     opts, &out, &count
//! );
//! ```
//!
//! # Dimension Safety
//!
//! These functions validate that the input quantity has an angle dimension and
//! return [`SiderustStatus::InvalidDimension`] if the check fails.

use crate::altitude::{periods_to_c, window_from_c};
use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use qtty_ffi::{QttyQuantity, UnitId};
use siderust::AltitudePeriodsProvider;
use tempoch::ModifiedJulianDate;

// Re-export QttyQuantity so callers can use it via siderust_ffi::QttyQuantity.
pub use qtty_ffi::QttyQuantity as SiderustQttyQuantity;

/// Helper to convert a qtty angle quantity to degrees.
/// Returns `Err(InvalidDimension)` if the quantity is not an angle.
fn qty_to_degrees(qty: QttyQuantity) -> Result<Degrees, SiderustStatus> {
    // Check if the unit is an angle unit
    let deg_result = qty.convert_to(UnitId::Degree);
    match deg_result {
        Some(converted) => Ok(Degrees::new(converted.value)),
        None => Err(SiderustStatus::InvalidDimension),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude, qtty-aware variants
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of any subject at an instant, with threshold as a `QttyQuantity`.
///
/// The `mjd` parameter is still a raw `f64` Modified Julian Date. Use the
/// tempoch-ffi functions to work with typed time values.
#[no_mangle]
pub extern "C" fn siderust_altitude_at_qty(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    mjd: f64,
    out: *mut QttyQuantity,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        dispatch_subject!(subject, |p| {
            let alt_rad = p
                .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd));
            let alt_deg = alt_rad.to::<Degree>();
            unsafe {
                *out = QttyQuantity::new(alt_deg.value(), UnitId::Degree);
            }
            SiderustStatus::Ok
        })
    }}
}

/// Periods when a subject is above a threshold altitude (qtty-aware).
///
/// The `threshold` parameter must be an angle quantity (degrees, radians, arcminutes, etc.).
/// Returns `InvalidDimension` if the quantity is not an angle.
#[no_mangle]
pub extern "C" fn siderust_above_threshold_qty(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold: QttyQuantity,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let threshold_deg = match qty_to_degrees(threshold) {
            Ok(d) => d,
            Err(e) => return e,
        };
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            periods_to_c(
                siderust::above_threshold(
                    p,
                    &observer.to_rust(),
                    window,
                    threshold_deg,
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

/// Periods when a subject is below a threshold altitude (qtty-aware).
///
/// The `threshold` parameter must be an angle quantity (degrees, radians, arcminutes, etc.).
/// Returns `InvalidDimension` if the quantity is not an angle.
#[no_mangle]
pub extern "C" fn siderust_below_threshold_qty(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold: QttyQuantity,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let threshold_deg = match qty_to_degrees(threshold) {
            Ok(d) => d,
            Err(e) => return e,
        };
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            periods_to_c(
                siderust::below_threshold(
                    p,
                    &observer.to_rust(),
                    window,
                    threshold_deg,
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    fn paris() -> SiderustGeodetict {
        SiderustGeodetict {
            lon_deg: 2.3522,
            lat_deg: 48.8566,
            height_m: 35.0,
        }
    }

    fn sun_subject() -> SiderustSubject {
        SiderustSubject::body(SiderustBody::Sun)
    }

    #[test]
    fn altitude_at_qty_sun() {
        let mut out = QttyQuantity::new(0.0, UnitId::Degree);
        let st = siderust_altitude_at_qty(sun_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert_eq!(out.unit, UnitId::Degree);
        assert!(out.value.is_finite());
    }

    #[test]
    fn above_threshold_qty_accepts_radians() {
        // Pass threshold in radians - should be converted to degrees internally
        let threshold = QttyQuantity::new(0.1745, UnitId::Radian); // ~10 degrees
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60001.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold_qty(
            sun_subject(),
            paris(),
            window,
            threshold,
            SiderustSearchOpts::default(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if !out.is_null() && count > 0 {
            unsafe { crate::ffi_utils::free_boxed_slice(out, count) };
        }
    }

    #[test]
    fn above_threshold_qty_rejects_non_angle() {
        // Pass a length quantity instead of an angle - should fail
        let threshold = QttyQuantity::new(100.0, UnitId::Meter);
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60001.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold_qty(
            sun_subject(),
            paris(),
            window,
            threshold,
            SiderustSearchOpts::default(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::InvalidDimension);
    }
}
