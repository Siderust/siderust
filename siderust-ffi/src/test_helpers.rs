// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared test helpers used across FFI test modules.
//!
//! This module is only compiled under `#[cfg(test)]`.

use crate::types::*;

/// Observer at Paris (lat 48.85°N, lon 2.35°E, 35 m).
pub fn paris() -> SiderustGeodetict {
    SiderustGeodetict {
        lon_deg: 2.35,
        lat_deg: 48.85,
        height_m: 35.0,
    }
}

/// A one-day time window starting at MJD 60000.
pub fn one_day_window() -> TempochPeriodMjd {
    TempochPeriodMjd {
        start_mjd: 60000.0,
        end_mjd: 60001.0,
    }
}

/// Default search options (tight tolerance, no explicit scan step).
pub fn default_opts() -> SiderustSearchOpts {
    SiderustSearchOpts {
        time_tolerance_days: 1e-9,
        scan_step_days: 0.0,
        has_scan_step: false,
    }
}

/// Vega's ICRS coordinates as a `SiderustSphericalDir`.
pub fn icrs_vega() -> SiderustSphericalDir {
    SiderustSphericalDir {
        polar_deg: 38.78,
        azimuth_deg: 279.23,
        frame: SiderustFrame::ICRS,
    }
}
