// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared FFI utilities: null-pointer guard macro, generic Vec→C converter and
//! generic free helper.  All other modules should use these rather than
//! duplicating the same patterns.

use crate::error::SiderustStatus;

// ═══════════════════════════════════════════════════════════════════════════
// Panic-catching FFI guard
// ═══════════════════════════════════════════════════════════════════════════

/// Wrap an FFI function body so that any Rust panic is caught and converted to
/// [`SiderustStatus::InternalPanic`] instead of unwinding across the `extern "C"`
/// boundary (which is undefined behaviour).
///
/// Usage:
///
/// ```ignore
/// #[no_mangle]
/// pub extern "C" fn siderust_some_function(args...) -> SiderustStatus {
///     ffi_guard! {
///         // … body that returns SiderustStatus …
///     }
/// }
/// ```
#[macro_export]
macro_rules! ffi_guard {
    ($body:block) => {{
        let result = ::std::panic::catch_unwind(::std::panic::AssertUnwindSafe(|| $body));
        match result {
            Ok(status) => status,
            Err(_) => $crate::error::SiderustStatus::InternalPanic,
        }
    }};
}

// ═══════════════════════════════════════════════════════════════════════════
// Null-pointer guard
// ═══════════════════════════════════════════════════════════════════════════

/// Early-return `NullPointer` if any of the listed pointer expressions is null.
///
/// ```ignore
/// check_out!(out);            // single pointer
/// check_out!(out, count);     // two pointers
/// ```
#[macro_export]
macro_rules! check_out {
    ($out:expr) => {
        if $out.is_null() {
            return $crate::error::SiderustStatus::NullPointer;
        }
    };
    ($out:expr, $count:expr) => {
        if $out.is_null() || $count.is_null() {
            return $crate::error::SiderustStatus::NullPointer;
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Body dispatch macro
// ═══════════════════════════════════════════════════════════════════════════

/// Dispatch an action to the concrete body type selected by a
/// [`SiderustBody`](crate::types::SiderustBody) discriminant.
///
/// The bound `$provider` will be an owned instance of the matching body type
/// (e.g. `Sun`, `Moon`, `solar_system::Mars`, …).
///
/// ```ignore
/// dispatch_body!(body, |b| {
///     b.altitude_at(&observer, mjd).value()
/// })
/// ```
#[macro_export]
macro_rules! dispatch_body {
    ($body:expr, |$provider:ident| $action:expr) => {
        match $body {
            $crate::types::SiderustBody::Sun => {
                let $provider = siderust::bodies::Sun;
                $action
            }
            $crate::types::SiderustBody::Moon => {
                let $provider = siderust::bodies::Moon;
                $action
            }
            $crate::types::SiderustBody::Mercury => {
                let $provider = siderust::bodies::solar_system::Mercury;
                $action
            }
            $crate::types::SiderustBody::Venus => {
                let $provider = siderust::bodies::solar_system::Venus;
                $action
            }
            $crate::types::SiderustBody::Mars => {
                let $provider = siderust::bodies::solar_system::Mars;
                $action
            }
            $crate::types::SiderustBody::Jupiter => {
                let $provider = siderust::bodies::solar_system::Jupiter;
                $action
            }
            $crate::types::SiderustBody::Saturn => {
                let $provider = siderust::bodies::solar_system::Saturn;
                $action
            }
            $crate::types::SiderustBody::Uranus => {
                let $provider = siderust::bodies::solar_system::Uranus;
                $action
            }
            $crate::types::SiderustBody::Neptune => {
                let $provider = siderust::bodies::solar_system::Neptune;
                $action
            }
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Subject dispatch macro
// ═══════════════════════════════════════════════════════════════════════════

/// Dispatch an action to the concrete provider identified by a
/// [`SiderustSubject`](crate::types::SiderustSubject).
///
/// For `Body` subjects, this expands into [`dispatch_body!`].
/// For `Star` and `Target`, `$provider` is a reference (`&Star` / `&ICRS`).
/// For `Icrs`, it is a reference to a local `ICRS` value.
///
/// **`$provider` is always a reference**, so pass it directly to free
/// functions (e.g. `above_threshold(p, …)`, **not** `&p`).  Method calls
/// (`p.altitude_at(…)`) auto-deref and work either way.
///
/// On null `star_handle` or `target_handle` the macro returns
/// `SiderustStatus::NullPointer`.
///
/// ```ignore
/// dispatch_subject!(subject, |p| {
///     siderust::above_threshold(p, &observer, window, threshold, opts)
/// })
/// ```
#[macro_export]
macro_rules! dispatch_subject {
    ($subject:expr, |$provider:ident| $action:expr) => {{
        let __subj = &$subject;
        match __subj.kind {
            $crate::types::SiderustSubjectKind::Body => {
                $crate::dispatch_body!(__subj.body, |__body_owned| {
                    let $provider = &__body_owned;
                    $action
                })
            }
            $crate::types::SiderustSubjectKind::Star => {
                if __subj.star_handle.is_null() {
                    return $crate::error::SiderustStatus::NullPointer;
                }
                let $provider = unsafe { &(*__subj.star_handle).inner };
                $action
            }
            $crate::types::SiderustSubjectKind::Icrs => {
                let __icrs_owned = siderust::coordinates::spherical::direction::ICRS::new(
                    qtty::Degrees::new(__subj.icrs_dir.azimuth_deg),
                    qtty::Degrees::new(__subj.icrs_dir.polar_deg),
                );
                let $provider = &__icrs_owned;
                $action
            }
            $crate::types::SiderustSubjectKind::Target => {
                if __subj.target_handle.is_null() {
                    return $crate::error::SiderustStatus::NullPointer;
                }
                let $provider = unsafe { &(*__subj.target_handle).dir };
                $action
            }
        }
    }};
}

// ═══════════════════════════════════════════════════════════════════════════
// FfiFrom trait
// ═══════════════════════════════════════════════════════════════════════════

/// Convert a Rust domain type into an FFI-compatible value.
///
/// Implement this trait for every FFI struct that mirrors a Rust type.  The
/// blanket call `U::ffi_from(r)` can then be passed directly to [`vec_to_c`].
pub trait FfiFrom<R>: Sized {
    fn ffi_from(r: &R) -> Self;
}

// ═══════════════════════════════════════════════════════════════════════════
// Vec → C array
// ═══════════════════════════════════════════════════════════════════════════

/// Convert a `Vec<R>` into a heap-allocated C array, writing the pointer and
/// element count into caller-supplied out-parameters.
///
/// Returns [`SiderustStatus::NullPointer`] if either output pointer is null.
///
/// # Safety
/// The returned array must be freed by the matching `siderust_*_free`
/// function, which reconstructs the boxed slice from `(ptr, count)` before
/// dropping it.
pub fn vec_to_c<R, U, F>(
    items: Vec<R>,
    conv: F,
    out: *mut *mut U,
    count: *mut usize,
) -> SiderustStatus
where
    U: Copy,
    F: Fn(&R) -> U,
{
    if out.is_null() || count.is_null() {
        return SiderustStatus::NullPointer;
    }
    let v: Vec<U> = items.iter().map(conv).collect();
    let len = v.len();
    unsafe {
        *out = Box::into_raw(v.into_boxed_slice()) as *mut _;
        *count = len;
    }
    SiderustStatus::Ok
}

// ═══════════════════════════════════════════════════════════════════════════
// Generic free
// ═══════════════════════════════════════════════════════════════════════════

/// Free a contiguous C array that was allocated by [`vec_to_c`].
///
/// # Safety
/// `ptr` and `count` must have been produced by the same call that returned
/// the array.  The pointer must not be used after this call.
pub unsafe fn free_boxed_slice<T>(ptr: *mut T, count: usize) {
    if !ptr.is_null() && count > 0 {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, count));
    }
}
