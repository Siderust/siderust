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
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out = SiderustGeodetict::from_rust(&siderust::observatories::ROQUE_DE_LOS_MUCHACHOS);
    }
    SiderustStatus::Ok
}

/// Fill `out` with the El Paranal observatory (Chile).
#[no_mangle]
pub extern "C" fn siderust_observatory_el_paranal(out: *mut SiderustGeodetict) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out = SiderustGeodetict::from_rust(&siderust::observatories::EL_PARANAL);
    }
    SiderustStatus::Ok
}

/// Fill `out` with the Mauna Kea observatory (Hawaiʻi, USA).
#[no_mangle]
pub extern "C" fn siderust_observatory_mauna_kea(out: *mut SiderustGeodetict) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out = SiderustGeodetict::from_rust(&siderust::observatories::MAUNA_KEA);
    }
    SiderustStatus::Ok
}

/// Fill `out` with the La Silla observatory (Chile).
#[no_mangle]
pub extern "C" fn siderust_observatory_la_silla(out: *mut SiderustGeodetict) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out = SiderustGeodetict::from_rust(&siderust::observatories::LA_SILLA_OBSERVATORY);
    }
    SiderustStatus::Ok
}

/// Create a custom geodetic position (WGS84/ECEF).
#[no_mangle]
pub extern "C" fn siderust_geodetic_new(
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
    out: *mut SiderustGeodetict,
) -> SiderustStatus {
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
}
