// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Orientation model primitives for rotating-body reference frames.

use crate::time::JulianDate;
use qtty::Degrees;

/// IAU pole and prime meridian rotation parameters for a body.
///
/// These parameters define a body's orientation at J2000.0 and its rates
/// of change. Time *T* is in Julian centuries of 36525 days from J2000.0 TDB,
/// and *d* is in days from J2000.0 TDB.
#[derive(Debug, Clone, Copy)]
pub struct IauRotationParams {
    /// Right ascension of the north pole at J2000.0, in degrees.
    pub alpha0_deg: Degrees,
    /// Rate of change of alpha0, in degrees per Julian century.
    pub alpha0_rate: Degrees,
    /// Declination of the north pole at J2000.0, in degrees.
    pub delta0_deg: Degrees,
    /// Rate of change of delta0, in degrees per Julian century.
    pub delta0_rate: Degrees,
    /// Prime meridian angle at J2000.0, in degrees.
    pub w0_deg: Degrees,
    /// Rate of change of W, in degrees per day.
    pub w_rate: Degrees,
}

impl IauRotationParams {
    /// Computes the pole right ascension alpha0 at a given time.
    ///
    /// # Arguments
    /// * `jd` - Julian Date (TDB/TT-compatible epoch for IAU model evaluation).
    #[inline]
    pub fn alpha0(&self, jd: JulianDate) -> Degrees {
        let t_centuries = jd.julian_centuries().value();
        self.alpha0_deg + self.alpha0_rate * t_centuries
    }

    /// Computes the pole declination delta0 at a given time.
    ///
    /// # Arguments
    /// * `jd` - Julian Date (TDB/TT-compatible epoch for IAU model evaluation).
    #[inline]
    pub fn delta0(&self, jd: JulianDate) -> Degrees {
        let t_centuries = jd.julian_centuries().value();
        self.delta0_deg + self.delta0_rate * t_centuries
    }

    /// Computes the prime meridian angle W at a given time.
    ///
    /// # Arguments
    /// * `jd` - Julian Date (TDB/TT-compatible epoch for IAU model evaluation).
    #[inline]
    pub fn w(&self, jd: JulianDate) -> Degrees {
        let d_days = (jd - JulianDate::J2000).value();
        self.w0_deg + self.w_rate * d_days
    }
}

/// Marker trait for bodies that expose IAU rotation parameters.
///
/// This is intentionally a regular trait with an associated const so it remains
/// zero-cost under monomorphization and stable on Rust.
pub trait HasIauRotation {
    /// IAU rotation parameters for the body.
    const ROTATION: IauRotationParams;
}
