// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression tests for the constops FFI surface.

#![cfg(feature = "constops")]

use std::ffi::CStr;

use siderust_ffi::constops_ffi::{
    constops_buffer_free, constops_last_error, constops_passes_generate, ConstopsBuffer,
};

const ISS_LINE_1: &str = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
const ISS_LINE_2: &str = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

fn empty_buffer() -> ConstopsBuffer {
    ConstopsBuffer {
        data: std::ptr::null_mut(),
        len: 0,
    }
}

fn buffer_to_vec(buf: &ConstopsBuffer) -> Vec<u8> {
    if buf.data.is_null() || buf.len == 0 {
        return Vec::new();
    }
    // Copy out into a freshly-owned Vec so we can free the FFI buffer.
    unsafe { std::slice::from_raw_parts(buf.data, buf.len) }.to_vec()
}

/// Anchor the horizon at the canonical Vallado TLE epoch
/// (2008-264.51782528 UT) so SGP4 reliably produces passes.
fn iss_horizon_epoch_seconds() -> (i64, i64) {
    // 2008-09-20 12:25:40.1043 UT, give or take fractional seconds.
    // Day-of-year 264 in 2008 → 2008-09-20.
    // 0.51782528 day = 12h 25m 40.1043s.
    // We use chrono via constops's own dep tree to keep this transparent.
    use chrono::{TimeZone, Utc};
    let start = Utc
        .with_ymd_and_hms(2008, 9, 20, 12, 25, 40)
        .single()
        .expect("valid epoch");
    let end = start + chrono::Duration::seconds(86_400);
    (start.timestamp(), end.timestamp())
}

#[test]
fn passes_generate_roundtrip() {
    let (start_secs, end_secs) = iss_horizon_epoch_seconds();
    // tempoch's serde representation for `Time<UTC>` is RFC3339 (per
    // constops::api::e2e tests). Build the request via serde_json::json! so we
    // do not have to depend on tempoch directly here.
    let start = chrono::DateTime::<chrono::Utc>::from_timestamp(start_secs, 0)
        .unwrap()
        .to_rfc3339_opts(chrono::SecondsFormat::Secs, true);
    let end = chrono::DateTime::<chrono::Utc>::from_timestamp(end_secs, 0)
        .unwrap()
        .to_rfc3339_opts(chrono::SecondsFormat::Secs, true);

    let req = serde_json::json!({
        "tenant": "acme",
        "spacecraft": [{
            "id": "ISS",
            "constellation": null,
            "norad": null,
            "tenant": "acme",
            "radios": [{
                "band": "S",
                "downlink_bps": 2_000_000.0,
                "supports_ttc": true,
                "supports_payload": true
            }],
            "allowed_modes": ["Nominal", "Payload"],
            "budget": {
                "battery_capacity_wh": 2000.0,
                "payload_power": 40.0,
                "bus_power": 60.0,
                "solar_power": 250.0,
                "storage_capacity": 68_719_476_736u64,
                "payload_data_rate": 1_000_000.0
            }
        }],
        "tle": [{
            "spacecraft_id": "ISS",
            "line1": ISS_LINE_1,
            "line2": ISS_LINE_2,
            "name": "ISS (ZARYA)"
        }],
        "ground_assets": [{
            "provider": "aws",
            "station": "bochum",
            "bands": ["S"],
            "min_elevation": 5.0,
            "max_elevation": 90.0,
            "latency": "LowLatency",
            "historical_success": 0.95,
            "location": {
                "longitude": 7.22,
                "latitude": 51.48,
                "height": 150.0
            }
        }],
        "horizon": { "start": start, "end": end }
    });
    let req_bytes = serde_json::to_vec(&req).expect("serialise request");

    let mut out = empty_buffer();
    let rc = unsafe {
        constops_passes_generate(req_bytes.as_ptr(), req_bytes.len(), &mut out as *mut _)
    };
    if rc != 0 {
        let msg_ptr = constops_last_error();
        let msg = if msg_ptr.is_null() {
            "<no message>".to_string()
        } else {
            unsafe { CStr::from_ptr(msg_ptr) }
                .to_string_lossy()
                .into_owned()
        };
        panic!("constops_passes_generate failed (rc={rc}): {msg}");
    }

    let body = buffer_to_vec(&out);
    unsafe { constops_buffer_free(out) };

    let resp: serde_json::Value = serde_json::from_slice(&body).expect("decode response");
    let passes = resp
        .get("passes")
        .and_then(|p| p.as_array())
        .expect("response has passes array");
    assert!(
        !passes.is_empty(),
        "expected at least one ISS pass over Bochum, got {}",
        serde_json::to_string_pretty(&resp).unwrap_or_default()
    );
}

#[test]
fn error_path_returns_nonzero_and_sets_message() {
    let bad = b"{ this is not valid json";
    let mut out = empty_buffer();
    let rc = unsafe { constops_passes_generate(bad.as_ptr(), bad.len(), &mut out as *mut _) };
    assert_ne!(rc, 0, "expected non-zero status for invalid JSON");
    assert!(
        out.data.is_null() && out.len == 0,
        "out buffer must remain empty on failure"
    );

    let msg_ptr = constops_last_error();
    assert!(
        !msg_ptr.is_null(),
        "constops_last_error must be set on failure"
    );
    let msg = unsafe { CStr::from_ptr(msg_ptr) }.to_string_lossy();
    assert!(!msg.is_empty(), "error message must be non-empty");
    assert!(
        msg.contains("JSON") || msg.contains("json"),
        "error message should mention JSON, got: {msg}"
    );
}
