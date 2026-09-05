// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI access to observatories generated from Siderust's canonical TOML catalog.

use crate::error::SiderustStatus;
use crate::types::SiderustGeodetict;
use siderust::catalogs::observatories;

/// Fill `out` with the Roque de los Muchachos observatory (La Palma, Spain).
#[no_mangle]
pub extern "C" fn siderust_observatory_roque_de_los_muchachos(
    out: *mut SiderustGeodetict,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: raw-pointer use follows this function's C ABI preconditions.
        unsafe {
            *out = SiderustGeodetict::from_rust(&observatories::ROQUE_DE_LOS_MUCHACHOS.geodetic);
        }
        SiderustStatus::Ok

    }}
}

/// Fill `out` with the El Paranal observatory (Chile).
#[no_mangle]
pub extern "C" fn siderust_observatory_el_paranal(out: *mut SiderustGeodetict) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: raw-pointer use follows this function's C ABI preconditions.
        unsafe {
            *out = SiderustGeodetict::from_rust(&observatories::EL_PARANAL.geodetic);
        }
        SiderustStatus::Ok

    }}
}

/// Fill `out` with the Mauna Kea observatory (Hawaiʻi, USA).
#[no_mangle]
pub extern "C" fn siderust_observatory_mauna_kea(out: *mut SiderustGeodetict) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: raw-pointer use follows this function's C ABI preconditions.
        unsafe {
            *out = SiderustGeodetict::from_rust(&observatories::MAUNA_KEA.geodetic);
        }
        SiderustStatus::Ok

    }}
}

/// Fill `out` with the La Silla observatory (Chile).
#[no_mangle]
pub extern "C" fn siderust_observatory_la_silla(out: *mut SiderustGeodetict) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: raw-pointer use follows this function's C ABI preconditions.
        unsafe {
            *out = SiderustGeodetict::from_rust(&observatories::LA_SILLA_OBSERVATORY.geodetic);
        }
        SiderustStatus::Ok

    }}
}

/// Create a custom geodetic position (WGS84/ECEF).
#[no_mangle]
pub extern "C" fn siderust_geodetic_new(
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
    out: *mut SiderustGeodetict,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: raw-pointer use follows this function's C ABI preconditions.
        unsafe {
            *out = SiderustGeodetict {
                lon_deg,
                lat_deg,
                height_m,
            };
        }
        SiderustStatus::Ok

    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    fn uninit_geodetic() -> SiderustGeodetict {
        SiderustGeodetict {
            lon_deg: 0.0,
            lat_deg: 0.0,
            height_m: 0.0,
        }
    }

    #[test]
    fn roque_de_los_muchachos_returns_ok() {
        let mut out = uninit_geodetic();
        let s = siderust_observatory_roque_de_los_muchachos(&mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.lon_deg, -17.8925);
        assert_eq!(out.lat_deg, 28.7543);
        assert_eq!(out.height_m, 2396.0);
    }

    #[test]
    fn null_ptr_roque() {
        assert_eq!(
            siderust_observatory_roque_de_los_muchachos(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn el_paranal_returns_ok() {
        let mut out = uninit_geodetic();
        assert_eq!(
            siderust_observatory_el_paranal(&mut out),
            SiderustStatus::Ok
        );
        assert_eq!(out.lon_deg, -70.4043);
        assert_eq!(out.lat_deg, -24.6272);
        assert_eq!(out.height_m, 2635.0);
    }

    #[test]
    fn null_ptr_el_paranal() {
        assert_eq!(
            siderust_observatory_el_paranal(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn mauna_kea_returns_ok() {
        let mut out = uninit_geodetic();
        assert_eq!(siderust_observatory_mauna_kea(&mut out), SiderustStatus::Ok);
        assert_eq!(out.lon_deg, -155.4681);
        assert_eq!(out.lat_deg, 19.8207);
        assert_eq!(out.height_m, 4207.0);
    }

    #[test]
    fn null_ptr_mauna_kea() {
        assert_eq!(
            siderust_observatory_mauna_kea(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn la_silla_returns_ok() {
        let mut out = uninit_geodetic();
        assert_eq!(siderust_observatory_la_silla(&mut out), SiderustStatus::Ok);
        assert_eq!(out.lon_deg, -70.7346);
        assert_eq!(out.lat_deg, -29.2584);
        assert_eq!(out.height_m, 2400.0);
    }

    #[test]
    fn null_ptr_la_silla() {
        assert_eq!(
            siderust_observatory_la_silla(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn geodetic_new_stores_values() {
        let mut out = uninit_geodetic();
        let s = siderust_geodetic_new(10.5, -20.3, 1500.0, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert!((out.lon_deg - 10.5).abs() < 1e-12);
        assert!((out.lat_deg - (-20.3)).abs() < 1e-12);
        assert!((out.height_m - 1500.0).abs() < 1e-12);
    }

    #[test]
    fn null_ptr_geodetic_new() {
        assert_eq!(
            siderust_geodetic_new(0.0, 0.0, 0.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }
}
