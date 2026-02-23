// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared FFI utilities: null-pointer guard macro, generic Vec→C converter and
//! generic free helper.  All other modules should use these rather than
//! duplicating the same patterns.

use crate::error::SiderustStatus;

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
