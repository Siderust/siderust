// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Body Orientation Primitives
//!
//! Type-level primitives for describing the orientation of rotating
//! Solar-System bodies via the standard IAU pole-and-prime-meridian model.
//!
//! ## Scientific scope
//!
//! The IAU Working Group on Cartographic Coordinates and Rotational Elements
//! recommends representing each body's orientation by the right ascension
//! and declination of its north pole at J2000.0 together with their secular
//! rates, and by the rotation angle `W` of the prime meridian as a function
//! of time. These six numbers define the body-fixed frame relative to the
//! ICRS, which is what is needed to convert between body-fixed and inertial
//! coordinates for planets, satellites and minor bodies.
//!
//! ## Technical scope
//!
//! [`IauRotationParams`] stores the six parameters in [`Degrees`] using the
//! IAU convention that pole rates are per Julian century while `W` evolves
//! per day, both measured from the J2000.0 TDB epoch. Helper methods
//! [`IauRotationParams::alpha0`], [`IauRotationParams::delta0`] and
//! [`IauRotationParams::w`] evaluate the linear models at a given Julian
//! Date. The [`HasIauRotation`] trait lets concrete body types expose their
//! parameters as an associated `const`, keeping the abstraction zero-cost
//! after monomorphisation.
//!
//! ## References
//!
//! * Archinal, B. A. et al., *Report of the IAU WG on Cartographic
//!   Coordinates and Rotational Elements*, Celest. Mech. Dyn. Astron.
//! * IAU 2009 / 2015 / 2018 reports on rotational elements

use crate::qtty::Degrees;
use crate::time::JulianDate;

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
        let t_centuries = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;
        self.alpha0_deg + self.alpha0_rate * t_centuries
    }

    /// Computes the pole declination delta0 at a given time.
    ///
    /// # Arguments
    /// * `jd` - Julian Date (TDB/TT-compatible epoch for IAU model evaluation).
    #[inline]
    pub fn delta0(&self, jd: JulianDate) -> Degrees {
        let t_centuries = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;
        self.delta0_deg + self.delta0_rate * t_centuries
    }

    /// Computes the prime meridian angle W at a given time.
    ///
    /// # Arguments
    /// * `jd` - Julian Date (TDB/TT-compatible epoch for IAU model evaluation).
    #[inline]
    pub fn w(&self, jd: JulianDate) -> Degrees {
        let d_days = (jd.raw() - JulianDate::J2000.raw()).value();
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
