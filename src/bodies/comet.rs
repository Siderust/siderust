// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Comet module
//!
//! ## Scientific scope
//!
//! This module models **comets** — small Solar System bodies on highly
//! eccentric orbits that develop a visible coma and tail when near the Sun.
//! Each comet is described by:
//! - a short human‑readable designation,
//! - an approximate dust/ion tail length,
//! - six Keplerian orbital elements with an epoch (see [`KeplerianOrbit`]), and
//! - a flag indicating whether elements are heliocentric or barycentric.
//!
//! | Constant | Designation | Reference frame | Period (yr) | Notes |
//! |----------|-------------|-----------------|-------------|-------|
//! | [`HALLEY`] | 1P/Halley | Heliocentric | 75.3 | Archetype of periodic comets |
//! | [`ENCKE`] | 2P/Encke | Heliocentric | 3.30 | Shortest‑period named comet |
//! | [`HALE_BOPP`] | C/1995 O1 | **Barycentric** | ~2 530 | Great comet of 1997 |
//!
//! ## Technical scope
//!
//! - [`Comet`] — runtime struct (lifetimed name, [`Kilometers`] tail, orbit, frame).
//! - [`OrbitFrame`] — enum discriminating heliocentric vs. barycentric elements.
//! - [`CometBuilder`] — builder for runtime construction.
//! - [`Comet::period_years`] — returns a typed [`Years`] value via Kepler's
//!   third law (`P = a^(3/2)` with `a` in AU, `P` in Julian years).
//!
//! ## References
//!
//! - JPL Small‑Body Database. <https://ssd.jpl.nasa.gov/tools/sbdb_query.html>
//! - Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann‑Bell.

use crate::astro::orbit::KeplerianOrbit;
use crate::qtty::{AstronomicalUnits, Degrees, Kilometers, Years};

/// Indicates whether orbital elements are given **with respect to the Solar‑System barycentre**
/// or the heliocentre (Sun‑centred).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OrbitFrame {
    /// Elements are relative to the Sun's centre of mass.
    Heliocentric,
    /// Elements are relative to the Solar System barycentre.
    Barycentric,
}

/// Represents a **Comet** characterised by its name, approximate tail length and
/// orbital elements in a specific reference frame.
#[derive(Clone, Debug)]
pub struct Comet<'a> {
    /// Short human-readable designation (e.g. `"1P/Halley"`).
    pub name: &'a str,
    /// Approximate maximum tail length in kilometres.
    pub tail_length: Kilometers,
    /// Keplerian orbital elements.
    pub orbit: KeplerianOrbit,
    /// Reference frame for the orbital elements.
    pub reference: OrbitFrame,
}

impl<'a> Comet<'a> {
    /// Compile‑time constructor for constant comets.
    pub const fn new_const(
        name: &'a str,
        tail_length: Kilometers,
        orbit: KeplerianOrbit,
        reference: OrbitFrame,
    ) -> Self {
        Self {
            name,
            tail_length,
            orbit,
            reference,
        }
    }

    /// Fluent builder entry‑point.
    pub fn builder() -> CometBuilder<'a> {
        CometBuilder::default()
    }
}

// -------------------------------------------------------------------------------------------------
//  Builder pattern
// -------------------------------------------------------------------------------------------------

/// Builder for runtime construction of [`Comet`] values.
#[derive(Default, Debug, Clone)]
pub struct CometBuilder<'a> {
    name: Option<&'a str>,
    tail_length: Option<Kilometers>,
    orbit: Option<KeplerianOrbit>,
    reference: Option<OrbitFrame>,
}

impl<'a> CometBuilder<'a> {
    /// Set the comet name/designation.
    pub fn name(mut self, name: &'a str) -> Self {
        self.name = Some(name);
        self
    }
    /// Set the tail length.
    pub fn tail_length(mut self, len: Kilometers) -> Self {
        self.tail_length = Some(len);
        self
    }
    /// Set the Keplerian orbital elements.
    pub fn orbit(mut self, orbit: KeplerianOrbit) -> Self {
        self.orbit = Some(orbit);
        self
    }
    /// Set the reference frame for orbital elements.
    pub fn reference(mut self, frame: OrbitFrame) -> Self {
        self.reference = Some(frame);
        self
    }

    /// Build the [`Comet`]; panics if required fields are missing.
    pub fn build(self) -> Comet<'a> {
        Comet {
            name: self.name.expect("missing name"),
            tail_length: self.tail_length.expect("missing tail_length"),
            orbit: self.orbit.expect("missing orbit"),
            reference: self.reference.unwrap_or(OrbitFrame::Heliocentric),
        }
    }
}

// -------------------------------------------------------------------------------------------------
//  Derived helpers
// -------------------------------------------------------------------------------------------------

impl<'a> Comet<'a> {
    /// Approximate orbital period as a typed [`Years`] value using Kepler’s third law.
    ///
    /// Uses `P = a^(3/2)` with `a` in AU and `P` in Julian years (μ_sun ≈ 1 in these units).
    ///
    /// # Returns
    ///
    /// A [`Years`] quantity equal to the orbital period.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use siderust::bodies::comet::HALLEY;
    /// let p = HALLEY.period_years();
    /// assert!((p.value() - 75.0).abs() < 2.0, "Halley period ~75 yr, got {}", p.value());
    /// ```
    pub fn period_years(&self) -> Years {
        // μ_sun ≈ 1 (when a in AU, P in years): P = sqrt(a^3)
        let a = self.orbit.shape().semi_major_axis();
        let a_val = a.value();
        Years::new(a_val * a_val.sqrt())
    }
}

// -------------------------------------------------------------------------------------------------
//  Famous comets (orbital elements from JPL Small‑Body DB, epoch J2000)
// -------------------------------------------------------------------------------------------------

/// **1P/Halley** – archetype periodic comet (heliocentric elements).
pub const HALLEY: Comet = Comet::new_const(
    "1P/Halley",
    Kilometers::new(1.0e7),
    KeplerianOrbit::new(
        AstronomicalUnits::new(17.834144),
        0.967_142_9,
        Degrees::new(162.262_69),
        Degrees::new(58.42008),
        Degrees::new(111.33249),
        Degrees::new(38.384),
        crate::J2000,
    ),
    OrbitFrame::Heliocentric,
);

/// **2P/Encke** – shortest‑period named comet (heliocentric).
pub const ENCKE: Comet = Comet::new_const(
    "2P/Encke",
    Kilometers::new(5.0e6),
    KeplerianOrbit::new(
        AstronomicalUnits::new(2.215_080),
        0.850_219,
        Degrees::new(11.780),
        Degrees::new(334.568),
        Degrees::new(186.545),
        Degrees::new(163.105),
        crate::J2000,
    ),
    OrbitFrame::Heliocentric,
);

/// **C/1995 O1 (Hale‑Bopp)** – great comet of 1997 (barycentric elements for long‑period accuracy).
pub const HALE_BOPP: Comet = Comet::new_const(
    "C/1995 O1 (Hale‑Bopp)",
    Kilometers::new(1.0e8),
    KeplerianOrbit::new(
        AstronomicalUnits::new(286.538),
        0.995_086,
        Degrees::new(89.425),
        Degrees::new(282.472),
        Degrees::new(130.589),
        Degrees::new(17.258),
        crate::J2000,
    ),
    OrbitFrame::Barycentric,
);

/// Preset catalogue of known comets (Halley, Encke, Hale-Bopp).
pub const COMET_PRESETS: &[&Comet] = &[&HALLEY, &ENCKE, &HALE_BOPP];
