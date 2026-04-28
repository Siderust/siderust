// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Twilight classification FFI.
//!
//! Mirrors [`siderust::TwilightPhase`] and [`siderust::twilight_classification`]
//! into a flat C ABI.

use crate::error::SiderustStatus;
use qtty::*;
use siderust::calculus::solar::classification::{twilight_classification, TwilightPhase};

ffi_enum! {
    /// Sky condition derived from the Sun's altitude.
    ///
    /// Boundaries are inclusive on the upper (brighter) side, following the
    /// IAU/USNO convention: the Sun exactly at `0°` is [`Civil`], exactly at
    /// `-6°` is [`Nautical`], etc.
    ///
    /// [`Civil`]: SiderustTwilightPhase::Civil
    /// [`Nautical`]: SiderustTwilightPhase::Nautical
    pub enum SiderustTwilightPhase {
        /// Sun above the horizon: altitude > 0°.
        Day = 0 => "day",
        /// Civil twilight: -6° < altitude ≤ 0°.
        Civil = 1 => "civil",
        /// Nautical twilight: -12° < altitude ≤ -6°.
        Nautical = 2 => "nautical",
        /// Astronomical twilight: -18° < altitude ≤ -12°.
        Astronomical = 3 => "astronomical",
        /// Full darkness: altitude ≤ -18°.
        Dark = 4 => "dark" | "night",
    }
}

impl SiderustTwilightPhase {
    /// Convert from the canonical Rust enum.
    #[inline]
    pub fn from_rust(p: TwilightPhase) -> Self {
        match p {
            TwilightPhase::Day => Self::Day,
            TwilightPhase::Civil => Self::Civil,
            TwilightPhase::Nautical => Self::Nautical,
            TwilightPhase::Astronomical => Self::Astronomical,
            TwilightPhase::Dark => Self::Dark,
        }
    }
}

/// Classify the sky condition from the Sun's altitude in degrees.
///
/// Writes the resulting phase to `out` and returns [`SiderustStatus::Ok`] on
/// success. Returns [`SiderustStatus::NullPointer`] if `out` is null.
///
/// The boundary convention follows IAU/USNO; see [`SiderustTwilightPhase`].
#[no_mangle]
pub extern "C" fn siderust_twilight_classification_deg(
    altitude_deg: f64,
    out: *mut SiderustTwilightPhase,
) -> SiderustStatus {
    ffi_guard! {{
        check_out!(out);
        let phase = twilight_classification(Degrees::new(altitude_deg));
        unsafe { *out = SiderustTwilightPhase::from_rust(phase) };
        SiderustStatus::Ok
    }}
}

/// Classify the sky condition from the Sun's altitude in radians.
///
/// Writes the resulting phase to `out` and returns [`SiderustStatus::Ok`] on
/// success. Returns [`SiderustStatus::NullPointer`] if `out` is null.
#[no_mangle]
pub extern "C" fn siderust_twilight_classification_rad(
    altitude_rad: f64,
    out: *mut SiderustTwilightPhase,
) -> SiderustStatus {
    ffi_guard! {{
        check_out!(out);
        let phase = twilight_classification(qtty::Radians::new(altitude_rad));
        unsafe { *out = SiderustTwilightPhase::from_rust(phase) };
        SiderustStatus::Ok
    }}
}

#[cfg(test)]
mod tests {
    use super::*;

    fn classify_deg(alt: f64) -> SiderustTwilightPhase {
        let mut out = SiderustTwilightPhase::Day;
        let s = siderust_twilight_classification_deg(alt, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        out
    }

    fn classify_rad(alt: f64) -> SiderustTwilightPhase {
        let mut out = SiderustTwilightPhase::Day;
        let s = siderust_twilight_classification_rad(alt, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        out
    }

    #[test]
    fn deg_phases() {
        assert_eq!(classify_deg(10.0), SiderustTwilightPhase::Day);
        assert_eq!(classify_deg(0.0), SiderustTwilightPhase::Civil);
        assert_eq!(classify_deg(-6.0), SiderustTwilightPhase::Nautical);
        assert_eq!(classify_deg(-12.0), SiderustTwilightPhase::Astronomical);
        assert_eq!(classify_deg(-18.0), SiderustTwilightPhase::Dark);
        assert_eq!(classify_deg(-30.0), SiderustTwilightPhase::Dark);
    }

    #[test]
    fn rad_phases() {
        assert_eq!(
            classify_rad(45.0_f64.to_radians()),
            SiderustTwilightPhase::Day
        );
        assert_eq!(
            classify_rad(-9.0_f64.to_radians()),
            SiderustTwilightPhase::Nautical
        );
        assert_eq!(
            classify_rad(-18.0_f64.to_radians()),
            SiderustTwilightPhase::Dark
        );
    }

    #[test]
    fn null_out_returns_null_pointer() {
        let s = siderust_twilight_classification_deg(0.0, std::ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
        let s = siderust_twilight_classification_rad(0.0, std::ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn parse_name_roundtrip() {
        assert_eq!(
            SiderustTwilightPhase::parse_name("nautical"),
            Some(SiderustTwilightPhase::Nautical)
        );
        assert_eq!(SiderustTwilightPhase::Day.as_str(), "day");
    }
}
