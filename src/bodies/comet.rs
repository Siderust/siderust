// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Comet module
//!
//! Represents comets with their name, tail length, and orbital parameters.
//! - name: String identifier for the comet.
//! - tail_length: Tail length in kilometers (Kilometers).
//! - orbit: Orbital parameters (see [Orbit]).
//!
//! | Constant | Designation | Reference frame | Period (yr) | Notes |
//! |----------|-------------|-----------------|-------------|-------|
//! | [`HALLEY`] | 1P/Halley | Heliocentric | 75.3 | Archetype of periodic comets |
//! | [`ENCKE`] | 2P/Encke | Heliocentric | 3.30 | Shortest‑period named comet |
//! | [`HALE_BOPP`] | C/1995 O1 | **Barycentric** | ~2 530 | Great comet of 1997 |
//!
//! ## Example
//! ```rust
//! use siderust::bodies::comet::{HALLEY, OrbitFrame, CometBuilder};
//! use qtty::*;
//!
//! println!("Halley's tail is at least {:?} long", HALLEY.tail_length);
//!
//! // Build a custom comet with barycentric elements
//! let custom = CometBuilder::default()
//!     .name("C/2025 A1 Customius")
//!     .tail_length(Kilometers::new(1.2e7))
//!     .reference(OrbitFrame::Barycentric)
//!     .orbit( HALLEY.orbit )
//!     .build();
//! ```
//!
//! ---

use crate::astro::orbit::Orbit;
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, Degrees, Kilometers};

/// Indicates whether orbital elements are given **with respect to the Solar‑System barycentre**
/// or the heliocentre (Sun‑centred).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OrbitFrame {
    Heliocentric,
    Barycentric,
}

/// Represents a **Comet** characterised by its name, approximate tail length and
/// orbital elements in a specific reference frame.
#[derive(Clone, Debug)]
pub struct Comet<'a> {
    pub name: &'a str,
    pub tail_length: Kilometers,
    pub orbit: Orbit,
    pub reference: OrbitFrame,
}

impl<'a> Comet<'a> {
    /// Compile‑time constructor for constant comets.
    pub const fn new_const(
        name: &'a str,
        tail_length: Kilometers,
        orbit: Orbit,
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

#[derive(Default, Debug, Clone)]
pub struct CometBuilder<'a> {
    name: Option<&'a str>,
    tail_length: Option<Kilometers>,
    orbit: Option<Orbit>,
    reference: Option<OrbitFrame>,
}

impl<'a> CometBuilder<'a> {
    pub fn name(mut self, name: &'a str) -> Self {
        self.name = Some(name);
        self
    }
    pub fn tail_length(mut self, len: Kilometers) -> Self {
        self.tail_length = Some(len);
        self
    }
    pub fn orbit(mut self, orbit: Orbit) -> Self {
        self.orbit = Some(orbit);
        self
    }
    pub fn reference(mut self, frame: OrbitFrame) -> Self {
        self.reference = Some(frame);
        self
    }

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
    /// Approximate orbital period in Julian years (`365.25 d`) using Kepler’s third law.
    pub fn period_years(&self) -> f64 {
        // μ_sun ≈ 1 (when a in AU, P in years): P = sqrt(a^3)
        let a = self.orbit.semi_major_axis.value();
        a.powf(1.5)
    }
}

// -------------------------------------------------------------------------------------------------
//  Famous comets (orbital elements from JPL Small‑Body DB, epoch J2000)
// -------------------------------------------------------------------------------------------------

/// **1P/Halley** – archetype periodic comet (heliocentric elements).
pub const HALLEY: Comet = Comet::new_const(
    "1P/Halley",
    Kilometers::new(1.0e7),
    Orbit::new(
        AstronomicalUnits::new(17.834144),
        0.967_142_9,
        Degrees::new(162.262_69),
        Degrees::new(58.42008),
        Degrees::new(111.33249),
        Degrees::new(38.384),
        JulianDate::J2000,
    ),
    OrbitFrame::Heliocentric,
);

/// **2P/Encke** – shortest‑period named comet (heliocentric).
pub const ENCKE: Comet = Comet::new_const(
    "2P/Encke",
    Kilometers::new(5.0e6),
    Orbit::new(
        AstronomicalUnits::new(2.215_080),
        0.850_219,
        Degrees::new(11.780),
        Degrees::new(334.568),
        Degrees::new(186.545),
        Degrees::new(163.105),
        JulianDate::J2000,
    ),
    OrbitFrame::Heliocentric,
);

/// **C/1995 O1 (Hale‑Bopp)** – great comet of 1997 (barycentric elements for long‑period accuracy).
pub const HALE_BOPP: Comet = Comet::new_const(
    "C/1995 O1 (Hale‑Bopp)",
    Kilometers::new(1.0e8),
    Orbit::new(
        AstronomicalUnits::new(286.538),
        0.995_086,
        Degrees::new(89.425),
        Degrees::new(282.472),
        Degrees::new(130.589),
        Degrees::new(17.258),
        JulianDate::J2000,
    ),
    OrbitFrame::Barycentric,
);

pub const COMET_PRESETS: &[&Comet] = &[&HALLEY, &ENCKE, &HALE_BOPP];
