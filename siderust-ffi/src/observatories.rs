// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for observatory constants.

use crate::error::SiderustStatus;
use crate::types::SiderustGeodetict;

/// Fill `out` with the Roque de los Muchachos observatory (La Palma, Spain).
#[no_mangle]
pub extern "C" fn siderust_observatory_roque_de_los_muchachos(
    out: *mut SiderustGeodetict,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe {
            *out = SiderustGeodetict::from_rust(&siderust::observatories::ROQUE_DE_LOS_MUCHACHOS);
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
        unsafe {
            *out = SiderustGeodetict::from_rust(&siderust::observatories::EL_PARANAL);
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
        unsafe {
            *out = SiderustGeodetict::from_rust(&siderust::observatories::MAUNA_KEA);
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
        unsafe {
            *out = SiderustGeodetict::from_rust(&siderust::observatories::LA_SILLA_OBSERVATORY);
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
        // La Palma ~-17.9°, ~28.8°
        assert!((out.lon_deg - (-17.89)).abs() < 1.0);
        assert!((out.lat_deg - 28.76).abs() < 1.0);
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
        // Paranal ~-70.4°, ~-24.6°
        assert!((out.lon_deg - (-70.4)).abs() < 1.0);
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
        // Mauna Kea ~-155.5°, ~19.8°
        assert!((out.lon_deg - (-155.5)).abs() < 1.0);
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
        // La Silla ~-70.7°, ~-29.3°
        assert!(out.lat_deg < 0.0, "La Silla is in the southern hemisphere");
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
