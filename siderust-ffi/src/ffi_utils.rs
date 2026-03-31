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
/// For `Star` and `GenericTarget`, `$provider` is a shared reference to the
/// underlying Rust value.
/// For `Icrs`, it is a reference to a local `ICRS` value.
///
/// **`$provider` is always a reference**, so pass it directly to free
/// functions (e.g. `above_threshold(p, …)`, **not** `&p`).  Method calls
/// (`p.altitude_at(…)`) auto-deref and work either way.
///
/// On null `star_handle` or `generic_target_handle` the macro returns
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
            $crate::types::SiderustSubjectKind::GenericTarget => {
                if __subj.generic_target_handle.is_null() {
                    return $crate::error::SiderustStatus::NullPointer;
                }
                // Dispatch on the inner position, not the whole CoordinateWithPM,
                // because AltitudePeriodsProvider is implemented for the position type.
                let $provider = unsafe { &(*__subj.generic_target_handle).inner.position };
                $action
            }
        }
    }};
}

// ═══════════════════════════════════════════════════════════════════════════
// FfiFrom trait (Rust → FFI)
// ═══════════════════════════════════════════════════════════════════════════

/// Convert a Rust domain type into an FFI-compatible value.
///
/// Implement this trait for every FFI struct that mirrors a Rust type.  The
/// blanket call `U::ffi_from(r)` can then be passed directly to [`vec_to_c`].
pub trait FfiFrom<R>: Sized {
    /// Convert a reference to a Rust domain value into `Self`.
    fn ffi_from(r: &R) -> Self;
}

// ═══════════════════════════════════════════════════════════════════════════
// TryFromFfi trait (FFI → Rust)
// ═══════════════════════════════════════════════════════════════════════════

/// Convert an FFI-compatible value into a Rust domain type.
///
/// This is the reverse of [`FfiFrom`]: it validates and converts a C-side
/// POD value into the corresponding rich Rust type, returning an appropriate
/// [`SiderustStatus`] on failure.
pub trait TryFromFfi<F>: Sized {
    /// Try to convert an FFI value into a Rust domain value.
    ///
    /// Returns `Ok(rust_value)` on success or `Err(status)` with a descriptive
    /// error code when the FFI value is invalid (e.g. unknown enum discriminant,
    /// out-of-range field).
    fn try_from_ffi(ffi: &F) -> Result<Self, SiderustStatus>;
}

// ═══════════════════════════════════════════════════════════════════════════
// Generic pointer / output helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Write a value through a non-null output pointer.
///
/// Returns [`SiderustStatus::NullPointer`] if `out` is null, otherwise writes
/// `value` and returns [`SiderustStatus::Ok`].
///
/// # Safety
///
/// `out` must be properly aligned and point to valid (possibly uninitialised)
/// memory for one `T`.  This is always the case for output parameters in the
/// siderust-ffi C API.
#[inline]
pub unsafe fn write_out<T>(out: *mut T, value: T) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { out.write(value) };
    SiderustStatus::Ok
}

/// Read through a non-null pointer, returning `NullPointer` on null.
///
/// # Safety
///
/// `ptr` must point to a valid, initialised `T` if it is non-null.
#[inline]
pub unsafe fn read_nonnull<'a, T>(ptr: *const T) -> Result<&'a T, SiderustStatus> {
    if ptr.is_null() {
        return Err(SiderustStatus::NullPointer);
    }
    Ok(unsafe { &*ptr })
}

// ═══════════════════════════════════════════════════════════════════════════
// Declarative FFI enum macro
// ═══════════════════════════════════════════════════════════════════════════

/// Define a `#[repr(i32)]` FFI enum from a single declarative table.
///
/// Generates the enum itself, `from_raw(i32) -> Option<Self>`, and
/// optionally `from_str(&str) -> Option<Self>` + `as_str() -> &'static str`
/// when string mappings are provided.
///
/// # Syntax
///
/// ```ignore
/// ffi_enum! {
///     /// Doc comment on the enum.
///     pub enum SiderustBody {
///         Sun     = 0 => "sun",
///         Moon    = 1 => "moon",
///         Mercury = 2 => "mercury",
///     }
/// }
/// ```
///
/// The `=> "canonical"` part is optional.  When omitted, `from_str`/`as_str`
/// are not generated.  String aliases can be added with `| "alias1" | "alias2"`:
///
/// ```ignore
/// ffi_enum! {
///     pub enum SiderustCenter {
///         Barycentric  = 1 => "barycentric" | "ssb" | "solar_system_barycenter",
///         Heliocentric = 2 => "heliocentric" | "solar" | "sun",
///     }
/// }
/// ```
#[macro_export]
macro_rules! ffi_enum {
    // ── With string mappings ────────────────────────────────────────────
    (
        $(#[$meta:meta])*
        $vis:vis enum $Name:ident {
            $(
                $(#[$vmeta:meta])*
                $Variant:ident = $disc:expr
                    => $canonical:literal $(| $alias:literal)*
            ),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[repr(i32)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        $vis enum $Name {
            $(
                $(#[$vmeta])*
                $Variant = $disc,
            )+
        }

        impl $Name {
            /// Decode a raw `i32` discriminant.
            ///
            /// Returns `None` if the value does not match any known variant.
            pub fn from_raw(raw: i32) -> Option<Self> {
                match raw {
                    $( $disc => Some(Self::$Variant), )+
                    _ => None,
                }
            }

            /// Parse a name string (case-insensitive).
            ///
            /// Returns `None` if the string does not match any known variant.
            pub fn parse_name(s: &str) -> Option<Self> {
                let lower = s.to_ascii_lowercase();
                match lower.as_str() {
                    $( $canonical $(| $alias)* => Some(Self::$Variant), )+
                    _ => None,
                }
            }

            /// Return the canonical string name of this variant.
            pub fn as_str(&self) -> &'static str {
                match self {
                    $( Self::$Variant => $canonical, )+
                }
            }
        }
    };

    // ── Without string mappings ─────────────────────────────────────────
    (
        $(#[$meta:meta])*
        $vis:vis enum $Name:ident {
            $(
                $(#[$vmeta:meta])*
                $Variant:ident = $disc:expr
            ),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[repr(i32)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        $vis enum $Name {
            $(
                $(#[$vmeta])*
                $Variant = $disc,
            )+
        }

        impl $Name {
            /// Decode a raw `i32` discriminant.
            ///
            /// Returns `None` if the value does not match any known variant.
            pub fn from_raw(raw: i32) -> Option<Self> {
                match raw {
                    $( $disc => Some(Self::$Variant), )+
                    _ => None,
                }
            }
        }
    };
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
        // SAFETY: caller guarantees `ptr` and `count` originate from the same
        // `vec_to_c` allocation and have not been freed before.
        let _ = unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(ptr, count)) };
    }
}
