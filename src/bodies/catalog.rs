// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Star Catalog Module
//!
//! This module defines a **compile‑time** catalog of some of the brightest and most
//! frequently referenced stars in the night sky.  Each entry is a `const`
//! [`Star`], fully initialised with its
//! * distance,
//! * mass,
//! * radius,
//! * bolometric luminosity, and
//! * **equatorial coordinates referenced to the J2000.0 epoch**.
//!
//! All quantities use strongly typed physical‑unit wrappers from the crate’s
//! [`units`] module, providing both type‑safety and convenience.  Because every
//! field is available at compile‑time, the constants can be inlined into
//! `const fn` contexts or embedded in `no_std` targets where dynamic allocation
//! is undesirable.
//!
//! ## Coordinate system
//! Positions are expressed in the Earth‑centred *mean equator and equinox of
//! J2000.0*.  Right‑ascension (RA) and declination (Dec) are stored in
//! [`Degrees`]; distance is expressed in [`LightYear`].  A star’s position is
//! wrapped in a [`Target`] tagged with [`JulianDate::J2000`] so that downstream
//! algorithms can perform precession or proper‑motion corrections when needed.
//!
//! ## Example
//! ```rust
//! use siderust::bodies::catalog::VEGA;
//! println!(
//!     "{} is {:.1} ly away and shines with {:.0} L☉",
//!     VEGA.name,
//!     VEGA.distance,
//!     VEGA.luminosity
//! );
//! ```
//!
//! ## Stars included
//! | Constant | Common / Bayer name | Constellation | V‑mag |
//! |----------|---------------------|---------------|-------|
//! | [`VEGA`] | Vega / α Lyrae | Lyra | +0.03 |
//! | [`POLARIS`] | Polaris / α UMi | Ursa Minor | +1.98 |
//! | [`SIRIUS`] | Sirius / α CMa | Canis Major | −1.46 |
//! | [`CANOPUS`] | Canopus / α Car | Carina | −0.72 |
//! | [`ARCTURUS`] | Arcturus / α Boo | Boötes | −0.04 |
//! | [`RIGEL`] | Rigel / β Ori | Orion | +0.12 |
//! | [`BETELGEUSE`] | Betelgeuse / α Ori | Orion | +0.50 (var) |
//! | [`PROCYON`] | Procyon / α CMi | Canis Minor | +0.34 |
//! | [`ALDEBARAN`] | Aldebaran / α Tau | Taurus | +0.86 |
//! | [`ALTAIR`] | Altair / α Aql | Aquila | +0.76 |
//!
//! ---

use super::Star;
use crate::coordinates::spherical::position::EquatorialMeanJ2000;
use crate::targets::CoordinateWithPM;
use crate::time::JulianDate;
use qtty::length::nominal::SolarRadiuses;
use qtty::*;

/// **Vega** (α Lyrae) — an A0 V star that defines the zero‑point of the Johnson
/// V photometric system.
///
/// * **Distance:** 25 ly
/// * **Mass:** 2.135 M☉
/// * **Radius:** 2.59 R☉
/// * **Luminosity:** 40 L☉
/// * **Coordinates (J2000):** RA 279.2347°, Dec +38.7837°
pub const VEGA: Star<'static> = Star::new_const(
    "Vega",
    LightYears::new(25.0),
    SolarMasses::new(2.135),
    SolarRadiuses::new(2.59),
    SolarLuminosities::new(40.12),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(38.7837),  // Dec (polar)
            Degrees::new(279.2347), // RA (azimuth)
            LightYears::new(25.0),
        ),
        JulianDate::J2000,
    ),
);

/// **Polaris** (α Ursae Minoris) — the current North Celestial Pole star and a
/// classical Cepheid variable.
///
/// * **Distance:** 433 ly
/// * **Mass:** 6.5 M☉
/// * **Radius:** 46 R☉
/// * **Luminosity:** 2 500 L☉
/// * **Coordinates (J2000):** RA 37.95456°, Dec +89.26411°
pub const POLARIS: Star<'static> = Star::new_const(
    "Polaris",
    LightYears::new(433.0),
    SolarMasses::new(6.5),
    SolarRadiuses::new(46.0),
    SolarLuminosities::new(2500.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(89.26410897), // Dec (polar)
            Degrees::new(37.95456067), // RA (azimuth)
            LightYears::new(433.0),
        ),
        JulianDate::J2000,
    ),
);

/// **Sirius A** (α Canis Majoris) — the brightest star in the night sky,
/// accompanied by the white‑dwarf Sirius B.
///
/// * **Distance:** 8.6 ly
/// * **Mass:** 2.063 M☉
/// * **Radius:** 1.713 R☉
/// * **Luminosity:** 24.7 L☉
/// * **Coordinates (J2000):** RA 101.28716°, Dec −16.71612°
pub const SIRIUS: Star<'static> = Star::new_const(
    "Sirius",
    LightYears::new(8.6),
    SolarMasses::new(2.063),
    SolarRadiuses::new(1.713),
    SolarLuminosities::new(24.7),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(-16.716115867), // Dec (polar)
            Degrees::new(101.28715533),  // RA (azimuth)
            LightYears::new(8.6),
        ),
        JulianDate::J2000,
    ),
);

/// **Canopus** (α Carinae) — a yellow‑white F‑type super‑giant and the second
/// brightest star in the night sky.
///
/// * **Distance:** 310 ly
/// * **Mass:** 8 M☉
/// * **Radius:** 71 R☉
/// * **Luminosity:** 13 600 L☉
/// * **Coordinates (J2000):** RA 95.98788°, Dec −52.69566°
pub const CANOPUS: Star<'static> = Star::new_const(
    "Canopus",
    LightYears::new(310.0),
    SolarMasses::new(8.0),
    SolarRadiuses::new(71.0),
    SolarLuminosities::new(13_600.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(-52.69566111), // Dec (polar)
            Degrees::new(95.98787778),  // RA (azimuth)
            LightYears::new(310.0),
        ),
        JulianDate::J2000,
    ),
);

/// **Arcturus** (α Boötis) — a K‑type red giant and the brightest star in the
/// northern celestial hemisphere.
///
/// * **Distance:** 36.7 ly
/// * **Mass:** 1.1 M☉
/// * **Radius:** 26 R☉
/// * **Luminosity:** 170 L☉
/// * **Coordinates (J2000):** RA 213.9153°, Dec +19.1825°
pub const ARCTURUS: Star<'static> = Star::new_const(
    "Arcturus",
    LightYears::new(36.7),
    SolarMasses::new(1.1),
    SolarRadiuses::new(26.0),
    SolarLuminosities::new(170.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(19.1825),  // Dec (polar)
            Degrees::new(213.9153), // RA (azimuth)
            LightYears::new(36.7),
        ),
        JulianDate::J2000,
    ),
);

/// **Rigel** (β Orionis) — a blue super‑giant that outshines most of the Milky
/// Way despite its great distance.
///
/// * **Distance:** 860 ly
/// * **Mass:** 17 M☉
/// * **Radius:** 78.9 R☉
/// * **Luminosity:** 120 000 L☉
/// * **Coordinates (J2000):** RA 78.63447°, Dec −08.20164°
pub const RIGEL: Star<'static> = Star::new_const(
    "Rigel",
    LightYears::new(860.0),
    SolarMasses::new(17.0),
    SolarRadiuses::new(78.9),
    SolarLuminosities::new(120_000.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(-8.20163889), // Dec (polar)
            Degrees::new(78.634467),   // RA (azimuth)
            LightYears::new(860.0),
        ),
        JulianDate::J2000,
    ),
);

/// **Betelgeuse** (α Orionis) — a red super‑giant nearing the end of its life;
/// expected to explode as a type‑II supernova within the next million years.
///
/// * **Distance:** 548 ly
/// * **Mass:** 11.6 M☉
/// * **Radius:** 724 R☉
/// * **Luminosity:** 14 000 L☉
/// * **Coordinates (J2000):** RA 88.79294°, Dec +07.40706°
pub const BETELGEUSE: Star<'static> = Star::new_const(
    "Betelgeuse",
    LightYears::new(548.0),
    SolarMasses::new(11.6),
    SolarRadiuses::new(724.0),
    SolarLuminosities::new(14_000.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(7.407064),  // Dec (polar)
            Degrees::new(88.792939), // RA (azimuth)
            LightYears::new(548.0),
        ),
        JulianDate::J2000,
    ),
);

/// **Procyon A** (α Canis Minoris) — a binary system consisting of a F5 IV–V
/// sub‑giant and a faint white‑dwarf companion.
///
/// * **Distance:** 11.5 ly
/// * **Mass:** 1.499 M☉
/// * **Radius:** 2.048 R☉
/// * **Luminosity:** 6.93 L☉
/// * **Coordinates (J2000):** RA 114.82549°, Dec +05.22499°
pub const PROCYON: Star<'static> = Star::new_const(
    "Procyon",
    LightYears::new(11.5),
    SolarMasses::new(1.499),
    SolarRadiuses::new(2.048),
    SolarLuminosities::new(6.93),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(5.224993),   // Dec (polar)
            Degrees::new(114.825493), // RA (azimuth)
            LightYears::new(11.5),
        ),
        JulianDate::J2000,
    ),
);

/// **Aldebaran** (α Tauri) — an orange K‑type giant lying in front of the Hyades
/// cluster.
///
/// * **Distance:** 65.1 ly
/// * **Mass:** 1.16 M☉
/// * **Radius:** 45.1 R☉
/// * **Luminosity:** 439 L☉
/// * **Coordinates (J2000):** RA 68.98016°, Dec +16.50930°
pub const ALDEBARAN: Star<'static> = Star::new_const(
    "Aldebaran",
    LightYears::new(65.1),
    SolarMasses::new(1.16),
    SolarRadiuses::new(45.1),
    SolarLuminosities::new(439.0),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(16.509302), // Dec (polar)
            Degrees::new(68.980163), // RA (azimuth)
            LightYears::new(65.1),
        ),
        JulianDate::J2000,
    ),
);

/// **Altair** (α Aquilae) — a fast‑rotating A‑type main‑sequence star that bulges
/// noticeably at its equator.
///
/// * **Distance:** 16.7 ly
/// * **Mass:** 1.86 M☉
/// * **Radius:** 1.79 R☉
/// * **Luminosity:** 10.6 L☉
/// * **Coordinates (J2000):** RA 297.69583°, Dec +08.86832°
pub const ALTAIR: Star<'static> = Star::new_const(
    "Altair",
    LightYears::new(16.7),
    SolarMasses::new(1.86),
    SolarRadiuses::new(1.79),
    SolarLuminosities::new(10.6),
    CoordinateWithPM::<EquatorialMeanJ2000<LightYear>>::new_static(
        EquatorialMeanJ2000::<LightYear>::new_raw(
            Degrees::new(8.868321),   // Dec (polar)
            Degrees::new(297.695827), // RA (azimuth)
            LightYears::new(16.7),
        ),
        JulianDate::J2000,
    ),
);
