// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C ABI bridge for the [`constops`] planning loop.
//!
//! The wire format is JSON (UTF-8): every entry point takes a JSON request
//! buffer and produces a heap-allocated JSON response buffer that the caller
//! must release via [`constops_buffer_free`]. The JSON shapes exactly mirror
//! the corresponding axum handlers in [`constops::api`] — the FFI calls the
//! refactored `*_core` functions directly, so no HTTP layer is involved.
//!
//! Every entry point returns a [`ConstopsErrorCode`] as `i32`. On non-zero
//! return, the caller can fetch a human-readable message for the calling
//! thread via [`constops_last_error`].
//!
//! Memory ownership:
//! * Input buffers are borrowed for the duration of the call only.
//! * Output [`ConstopsBuffer`]s are produced by `Box::leak` of a boxed slice;
//!   they MUST be returned via [`constops_buffer_free`] (which reconstructs
//!   the box and drops it).

#![deny(unsafe_op_in_unsafe_fn)]
#![allow(clippy::missing_safety_doc)]

use std::cell::RefCell;
use std::ffi::{c_char, CString};
use std::sync::Arc;

use constops::api::{
    generate_passes_core, ingest_orbit_state_core, score_links_core, solve_schedule_core,
    validate_schedule_core, GeneratePassesRequest, IngestOrbitStateRequest, ScoreLinksRequest,
    SolveScheduleRequest, ValidateScheduleRequest,
};
use constops::ingest::OrbitRegistry;
use constops::model::ObjectiveWeights;

// ───────────────────────────── Status enum ──────────────────────────────────

/// Status codes returned by every `constops_*` FFI entry point.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstopsErrorCode {
    /// Success.
    Ok = 0,
    /// JSON encoding/decoding failure or malformed UTF-8.
    Json = 1,
    /// Pass generation failed (TLE, horizon, propagation).
    Pass = 2,
    /// Link scoring failed.
    Score = 3,
    /// Scheduler failed.
    Schedule = 4,
    /// Schedule validation failed.
    Validate = 5,
    /// Orbit-state ingest failed.
    Ingest = 6,
    /// Internal error (null pointer, panic, runtime construction failure).
    Internal = 99,
}

// ───────────────────────────── Output buffer ────────────────────────────────

/// Heap-allocated UTF-8 byte buffer returned to C callers.
///
/// `data` is owned by the FFI allocator until the caller releases it via
/// [`constops_buffer_free`]. A zero-length response uses `data == NULL` and
/// `len == 0`.
#[repr(C)]
pub struct ConstopsBuffer {
    /// Pointer to the first byte (or NULL when `len == 0`).
    pub data: *mut u8,
    /// Length in bytes.
    pub len: usize,
}

impl ConstopsBuffer {
    fn empty() -> Self {
        Self { data: std::ptr::null_mut(), len: 0 }
    }

    fn from_vec(mut v: Vec<u8>) -> Self {
        if v.is_empty() {
            return Self::empty();
        }
        v.shrink_to_fit();
        let boxed: Box<[u8]> = v.into_boxed_slice();
        let len = boxed.len();
        let data = Box::into_raw(boxed) as *mut u8;
        Self { data, len }
    }
}

// ───────────────────────────── Thread-local error ───────────────────────────

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
}

fn set_last_error<E: std::fmt::Display>(err: E) {
    let msg = err.to_string();
    let cstring = CString::new(msg.replace('\0', "?"))
        .unwrap_or_else(|_| CString::new("constops: error message contained NUL").unwrap());
    LAST_ERROR.with(|slot| *slot.borrow_mut() = Some(cstring));
}

fn clear_last_error() {
    LAST_ERROR.with(|slot| *slot.borrow_mut() = None);
}

// ───────────────────────────── Internal helpers ─────────────────────────────

/// Decode the input slice into `T` and capture the error message on failure.
fn decode_request<'a, T: serde::de::DeserializeOwned>(
    bytes: &'a [u8],
) -> Result<T, ConstopsErrorCode> {
    serde_json::from_slice::<T>(bytes).map_err(|e| {
        set_last_error(format_args!("invalid JSON request: {e}"));
        ConstopsErrorCode::Json
    })
}

/// Serialise `value` and write it through `out_ptr`.
fn write_response<T: serde::Serialize>(
    value: &T,
    out_ptr: *mut ConstopsBuffer,
) -> ConstopsErrorCode {
    let bytes = match serde_json::to_vec(value) {
        Ok(b) => b,
        Err(e) => {
            set_last_error(format_args!("response serialisation failed: {e}"));
            return ConstopsErrorCode::Json;
        }
    };
    // SAFETY: caller-supplied non-null, properly aligned `*mut ConstopsBuffer`
    // (validated by the entry point before we get here).
    unsafe { out_ptr.write(ConstopsBuffer::from_vec(bytes)) };
    ConstopsErrorCode::Ok
}

/// Validate the input pointers and return the borrowed request slice.
///
/// # Safety
/// `req` and `out` semantics match the public contract documented on each
/// entry point.
unsafe fn prepare<'a>(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> Result<&'a [u8], ConstopsErrorCode> {
    if out.is_null() {
        set_last_error("output pointer is null");
        return Err(ConstopsErrorCode::Internal);
    }
    // SAFETY: out is non-null (just checked).
    unsafe { out.write(ConstopsBuffer::empty()) };
    if req.is_null() && req_len != 0 {
        set_last_error("request pointer is null but length is non-zero");
        return Err(ConstopsErrorCode::Internal);
    }
    if req.is_null() {
        return Ok(&[]);
    }
    // SAFETY: caller asserts (req, req_len) describes a readable byte range.
    Ok(unsafe { std::slice::from_raw_parts(req, req_len) })
}

/// Run `body` while catching panics; on panic, set Internal error.
fn guard<F: FnOnce() -> ConstopsErrorCode>(body: F) -> ConstopsErrorCode {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(body)) {
        Ok(code) => code,
        Err(_) => {
            set_last_error("panic caught at constops FFI boundary");
            ConstopsErrorCode::Internal
        }
    }
}

// ───────────────────────────── Entry points ─────────────────────────────────

/// Generate satellite passes for a `GeneratePassesRequest` JSON payload.
///
/// On success returns 0 ([`ConstopsErrorCode::Ok`]) and writes a JSON-encoded
/// `GeneratePassesResponse` into `*out`. On failure returns a non-zero
/// [`ConstopsErrorCode`] and leaves `*out` zero-initialised.
///
/// # Safety
/// * `req` must point to `req_len` readable bytes (or be NULL with `req_len == 0`).
/// * `out` must point to a writable `ConstopsBuffer`.
/// * The caller must release `*out` via [`constops_buffer_free`] on success.
#[no_mangle]
pub unsafe extern "C" fn constops_passes_generate(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> i32 {
    guard(|| {
        clear_last_error();
        // SAFETY: forwarded contract from this function's docs.
        let bytes = match unsafe { prepare(req, req_len, out) } {
            Ok(b) => b,
            Err(code) => return code,
        };
        let request: GeneratePassesRequest = match decode_request(bytes) {
            Ok(r) => r,
            Err(code) => return code,
        };
        match generate_passes_core(request) {
            Ok(resp) => write_response(&resp, out),
            Err(e) => {
                set_last_error(&e);
                ConstopsErrorCode::Pass
            }
        }
    }) as i32
}

/// Score link opportunities for a `ScoreLinksRequest` JSON payload.
///
/// # Safety
/// See [`constops_passes_generate`].
#[no_mangle]
pub unsafe extern "C" fn constops_links_score(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> i32 {
    guard(|| {
        clear_last_error();
        // SAFETY: forwarded contract from this function's docs.
        let bytes = match unsafe { prepare(req, req_len, out) } {
            Ok(b) => b,
            Err(code) => return code,
        };
        let request: ScoreLinksRequest = match decode_request(bytes) {
            Ok(r) => r,
            Err(code) => return code,
        };
        match score_links_core(request, ObjectiveWeights::default()) {
            Ok(resp) => write_response(&resp, out),
            Err(e) => {
                set_last_error(&e);
                ConstopsErrorCode::Score
            }
        }
    }) as i32
}

/// Solve a schedule for a `SolveScheduleRequest` JSON payload.
///
/// # Safety
/// See [`constops_passes_generate`].
#[no_mangle]
pub unsafe extern "C" fn constops_schedule_solve(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> i32 {
    guard(|| {
        clear_last_error();
        // SAFETY: forwarded contract from this function's docs.
        let bytes = match unsafe { prepare(req, req_len, out) } {
            Ok(b) => b,
            Err(code) => return code,
        };
        let request: SolveScheduleRequest = match decode_request(bytes) {
            Ok(r) => r,
            Err(code) => return code,
        };
        match solve_schedule_core(request, ObjectiveWeights::default()) {
            Ok(resp) => write_response(&resp, out),
            Err(e) => {
                set_last_error(&e);
                ConstopsErrorCode::Schedule
            }
        }
    }) as i32
}

/// Validate a schedule for a `ValidateScheduleRequest` JSON payload.
///
/// The response body is a JSON-encoded `ValidationReport`.
///
/// # Safety
/// See [`constops_passes_generate`].
#[no_mangle]
pub unsafe extern "C" fn constops_schedule_validate(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> i32 {
    guard(|| {
        clear_last_error();
        // SAFETY: forwarded contract from this function's docs.
        let bytes = match unsafe { prepare(req, req_len, out) } {
            Ok(b) => b,
            Err(code) => return code,
        };
        let request: ValidateScheduleRequest = match decode_request(bytes) {
            Ok(r) => r,
            Err(code) => return code,
        };
        match validate_schedule_core(request) {
            Ok(report) => write_response(&report, out),
            Err(e) => {
                set_last_error(&e);
                ConstopsErrorCode::Validate
            }
        }
    }) as i32
}

/// Ingest an orbit state for an `IngestOrbitStateRequest` JSON payload.
///
/// Each call uses a fresh ephemeral [`OrbitRegistry`]; the returned
/// `stored_id` is monotonic within that registry only. Cross-call lineage is
/// not preserved by the FFI (use the REST API for persistent lineage).
///
/// # Safety
/// See [`constops_passes_generate`].
#[no_mangle]
pub unsafe extern "C" fn constops_orbit_ingest(
    req: *const u8,
    req_len: usize,
    out: *mut ConstopsBuffer,
) -> i32 {
    guard(|| {
        clear_last_error();
        // SAFETY: forwarded contract from this function's docs.
        let bytes = match unsafe { prepare(req, req_len, out) } {
            Ok(b) => b,
            Err(code) => return code,
        };
        let request: IngestOrbitStateRequest = match decode_request(bytes) {
            Ok(r) => r,
            Err(code) => return code,
        };
        let runtime = match tokio::runtime::Builder::new_current_thread().enable_all().build() {
            Ok(rt) => rt,
            Err(e) => {
                set_last_error(format_args!("tokio runtime build failed: {e}"));
                return ConstopsErrorCode::Internal;
            }
        };
        let registry = Arc::new(OrbitRegistry::new());
        let result = runtime.block_on(ingest_orbit_state_core(&registry, request));
        match result {
            Ok(resp) => write_response(&resp, out),
            Err(e) => {
                set_last_error(&e);
                ConstopsErrorCode::Ingest
            }
        }
    }) as i32
}

// ───────────────────────────── Free / error ─────────────────────────────────

/// Release a [`ConstopsBuffer`] previously returned by a `constops_*` call.
///
/// Calling this on a zero-initialised buffer (`data == NULL`, `len == 0`) is a
/// no-op. Calling it twice on the same buffer is undefined behaviour.
///
/// # Safety
/// `buf.data` must have been produced by a successful `constops_*` call and
/// not previously freed.
#[no_mangle]
pub unsafe extern "C" fn constops_buffer_free(buf: ConstopsBuffer) {
    if !buf.data.is_null() && buf.len > 0 {
        // SAFETY: contract states `(data, len)` came from `Box::into_raw` of a
        // boxed slice produced by `ConstopsBuffer::from_vec`.
        let _ = unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(buf.data, buf.len)) };
    }
}

/// Return the last error message recorded on the calling thread.
///
/// The returned pointer is owned by the FFI and remains valid until the next
/// `constops_*` call on the same thread or process exit. Returns NULL when
/// no error is recorded.
#[no_mangle]
pub extern "C" fn constops_last_error() -> *const c_char {
    LAST_ERROR.with(|slot| match slot.borrow().as_ref() {
        Some(cstr) => cstr.as_ptr(),
        None => std::ptr::null(),
    })
}
