// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Opaque transformation context handles for the C ABI.
//!
//! The current FFI context stores the runtime Earth-orientation model preset.
//! Additional runtime providers can be layered onto the same handle shape later
//! without forcing per-model symbol families in the public ABI.

use crate::error::SiderustStatus;
use crate::types::SiderustEarthOrientationModel;

/// Opaque FFI context used by model-sensitive transform entry points.
pub struct SiderustContext {
    pub(crate) model: SiderustEarthOrientationModel,
}

impl Default for SiderustContext {
    fn default() -> Self {
        Self {
            model: SiderustEarthOrientationModel::Iau2006A,
        }
    }
}

/// Create a context using the default siderust transform model.
#[no_mangle]
pub extern "C" fn siderust_context_create_default(
    out: *mut *mut SiderustContext,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let ctx = Box::new(SiderustContext::default());
        unsafe { *out = Box::into_raw(ctx) };
        SiderustStatus::Ok
    }}
}

/// Create a context with an explicit Earth-orientation model preset.
#[no_mangle]
pub extern "C" fn siderust_context_create_with_model(
    model: SiderustEarthOrientationModel,
    out: *mut *mut SiderustContext,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let ctx = Box::new(SiderustContext { model });
        unsafe { *out = Box::into_raw(ctx) };
        SiderustStatus::Ok
    }}
}

/// Free a context handle previously returned by `siderust_context_create_*`.
///
/// # Safety
/// `handle` must be either null or a live pointer produced by this crate.
#[no_mangle]
pub unsafe extern "C" fn siderust_context_free(handle: *mut SiderustContext) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Retrieve the model preset stored in a context handle.
#[no_mangle]
pub extern "C" fn siderust_context_get_model(
    handle: *const SiderustContext,
    out: *mut SiderustEarthOrientationModel,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).model };
        SiderustStatus::Ok
    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    #[test]
    fn create_default_and_read_back_model() {
        let mut handle: *mut SiderustContext = ptr::null_mut();
        assert_eq!(
            siderust_context_create_default(&mut handle),
            SiderustStatus::Ok
        );
        assert!(!handle.is_null());

        let mut model = SiderustEarthOrientationModel::Iau2000A;
        assert_eq!(
            siderust_context_get_model(handle, &mut model),
            SiderustStatus::Ok
        );
        assert_eq!(model, SiderustEarthOrientationModel::Iau2006A);

        unsafe { siderust_context_free(handle) };
    }

    #[test]
    fn create_with_model() {
        let mut handle: *mut SiderustContext = ptr::null_mut();
        assert_eq!(
            siderust_context_create_with_model(SiderustEarthOrientationModel::Iau2006, &mut handle),
            SiderustStatus::Ok
        );
        assert!(!handle.is_null());

        let mut model = SiderustEarthOrientationModel::Iau2000A;
        assert_eq!(
            siderust_context_get_model(handle, &mut model),
            SiderustStatus::Ok
        );
        assert_eq!(model, SiderustEarthOrientationModel::Iau2006);

        unsafe { siderust_context_free(handle) };
    }

    #[test]
    fn null_pointer_checks() {
        assert_eq!(
            siderust_context_create_default(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_context_create_with_model(
                SiderustEarthOrientationModel::Iau2006A,
                ptr::null_mut()
            ),
            SiderustStatus::NullPointer
        );

        let mut handle: *mut SiderustContext = ptr::null_mut();
        assert_eq!(
            siderust_context_create_default(&mut handle),
            SiderustStatus::Ok
        );

        let mut model = SiderustEarthOrientationModel::Iau2000A;
        assert_eq!(
            siderust_context_get_model(ptr::null(), &mut model),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_context_get_model(handle, ptr::null_mut()),
            SiderustStatus::NullPointer
        );

        unsafe { siderust_context_free(handle) };
        unsafe { siderust_context_free(ptr::null_mut()) };
    }
}
