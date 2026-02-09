// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Julian Date (`Time<JD>`) specific extensions.

use qtty::*;
use std::ops::Add;

use super::instant::Time;
use super::scales::{JD, MJD};

impl Time<JD> {
    /// J2000.0 epoch: 2000-01-01T12:00:00 TT  (JD 2 451 545.0).
    pub const J2000: Self = Self::new(2_451_545.0);

    /// One Julian year expressed in days.
    pub const JULIAN_YEAR: Days = Days::new(365.25);

    /// One Julian century expressed in days.
    pub const JULIAN_CENTURY: Days = Days::new(36_525.0);

    /// One Julian millennium expressed in days.
    pub const JULIAN_MILLENNIUM: Days = Days::new(365_250.0);

    /// Julian millennia since J2000.0 (used by VSOP87).
    #[inline]
    pub fn julian_millennias(&self) -> Millennia {
        Millennia::new(
            ((*self - Self::J2000) / Self::JULIAN_MILLENNIUM)
                .simplify()
                .value(),
        )
    }

    /// Julian centuries since J2000.0 (used by nutation, precession, sidereal time).
    #[inline]
    pub fn julian_centuries(&self) -> Centuries {
        Centuries::new(
            ((*self - Self::J2000) / Self::JULIAN_CENTURY)
                .simplify()
                .value(),
        )
    }

    /// Julian years since J2000.0.
    #[inline]
    pub fn julian_years(&self) -> JulianYears {
        JulianYears::new(
            ((*self - Self::J2000) / Self::JULIAN_YEAR)
                .simplify()
                .value(),
        )
    }

    /// Converts JD(TT) → JD(TDB) by applying the Fairhead & Bretagnon
    /// periodic correction (≈ 1.7 ms amplitude).
    ///
    /// This is a function-level correction; the `TDB` and `TT` *markers*
    /// are numerically identical on the JD axis for the purpose of epoch
    /// conversions, but this method gives the exact value when you need
    /// sub-millisecond accuracy.
    pub fn tt_to_tdb(jd_tt: Self) -> Self {
        let e = (357.53 + 0.985_600_28 * (jd_tt - Self::J2000).value()).to_radians();
        let delta_t = Days::new((1.658e-3 * e.sin() + 1.4e-6 * (2.0 * e).sin()) / 86_400.0);
        jd_tt + delta_t
    }

    /// Convenience: MJD value corresponding to this JD.
    ///
    /// Kept as a convenience wrapper for `self.to::<MJD>()`.
    #[inline]
    pub fn to_mjd(&self) -> Time<MJD> {
        self.to::<MJD>()
    }
}

// ── From / Into conversions for qtty time quantities (on JD only) ────────

impl Add<Years> for Time<JD> {
    type Output = Self;
    fn add(self, years: Years) -> Self {
        // Treat `Years` here as Julian years for JD arithmetic stability.
        self + Days::new(years.value() * Self::JULIAN_YEAR.value())
    }
}

impl From<JulianYears> for Time<JD> {
    fn from(years: JulianYears) -> Self {
        Self::J2000 + years.to::<Day>()
    }
}

impl From<Time<JD>> for JulianYears {
    fn from(jd: Time<JD>) -> Self {
        jd.julian_years()
    }
}

impl From<Centuries> for Time<JD> {
    fn from(centuries: Centuries) -> Self {
        // `Centuries` are interpreted as Julian centuries relative to J2000.
        Self::J2000 + Days::new(centuries.value() * Self::JULIAN_CENTURY.value())
    }
}

impl From<Time<JD>> for Centuries {
    fn from(jd: Time<JD>) -> Self {
        jd.julian_centuries()
    }
}

impl From<Millennia> for Time<JD> {
    fn from(millennia: Millennia) -> Self {
        // `Millennia` are interpreted as Julian millennia relative to J2000.
        Self::J2000 + Days::new(millennia.value() * Self::JULIAN_MILLENNIUM.value())
    }
}

impl From<Time<JD>> for Millennia {
    fn from(jd: Time<JD>) -> Self {
        jd.julian_millennias()
    }
}
