// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Programmatic constructor for [`TLE`] records.
//!
//! [`TleBuilder`] is the synthesis counterpart of the parser in
//! [`crate::formats::tle::parse_tle`]. It is used by `siderust-sgp4` tests and by any
//! caller that needs to fabricate plausible TLEs (regression suites,
//! property tests, scenario generators).
//!
//! All public setters take **typed** quantities; there is no `f64`-only
//! entry point for fields that have a typed representation.

use chrono::{Datelike, Timelike};
use qtty::angular::Degrees;
use qtty::angular::Turn;
use qtty::angular_rate::AngularRate;
use qtty::time::Day;
use tempoch::{Time, UTC};

use super::parse::compute_tle_checksum;
use super::TleError;
use crate::formats::tle::{Classification, InternationalDesignator, SatelliteNumber, TLE};

/// Typed builder for [`TLE`].
///
/// Required fields are `norad_id`, `international_designator`, `epoch`,
/// `inclination`, `raan`, `eccentricity`, `argument_of_perigee`,
/// `mean_anomaly`, and `mean_motion`. Optional fields default to sensible
/// values (`Classification::Unclassified`, zero drag terms, element-set
/// and revolution counters of zero, and no name).
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{Classification, InternationalDesignator, SatelliteNumber, TleBuilder};
/// use qtty::angular::Degrees;
/// use qtty::angular::Turn;
/// use qtty::angular_rate::AngularRate;
/// use qtty::time::Day;
/// use tempoch::{Time, UTC};
/// use chrono::{TimeZone, Utc};
///
/// let dt = Utc.with_ymd_and_hms(2008, 9, 20, 12, 25, 40).unwrap();
/// let tle = TleBuilder::new()
///     .norad_id(SatelliteNumber(25_544))
///     .international_designator(InternationalDesignator("98067A".into()))
///     .classification(Classification::Unclassified)
///     .epoch(Time::<UTC>::try_from_chrono(dt).unwrap())
///     .inclination(Degrees::new(51.6416))
///     .raan(Degrees::new(247.4627))
///     .eccentricity(0.0006703)
///     .argument_of_perigee(Degrees::new(130.5360))
///     .mean_anomaly(Degrees::new(325.0288))
///     .mean_motion(AngularRate::<Turn, Day>::new(15.72125391))
///     .name("ISS (ZARYA)")
///     .build()
///     .unwrap();
/// assert_eq!(tle.norad_id.0, 25_544);
/// ```
#[derive(Default, Clone, Debug)]
pub struct TleBuilder {
    name: Option<String>,
    norad_id: Option<SatelliteNumber>,
    classification: Option<Classification>,
    international_designator: Option<InternationalDesignator>,
    epoch: Option<Time<UTC>>,
    mean_motion_dot: f64,
    mean_motion_ddot: f64,
    bstar: f64,
    element_set_number: u16,
    revolution_number_at_epoch: u32,
    inclination: Option<Degrees>,
    raan: Option<Degrees>,
    eccentricity: Option<f64>,
    argument_of_perigee: Option<Degrees>,
    mean_anomaly: Option<Degrees>,
    mean_motion: Option<AngularRate<Turn, Day>>,
}

impl TleBuilder {
    /// Construct an empty builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the optional satellite name (line 0 of a 3LE).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().name("ISS (ZARYA)");
    /// ```
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Set the NORAD catalog id.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{SatelliteNumber, TleBuilder};
    /// let _b = TleBuilder::new().norad_id(SatelliteNumber(25_544));
    /// ```
    pub fn norad_id(mut self, n: SatelliteNumber) -> Self {
        self.norad_id = Some(n);
        self
    }

    /// Set the classification character.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{Classification, TleBuilder};
    /// let _b = TleBuilder::new().classification(Classification::Unclassified);
    /// ```
    pub fn classification(mut self, c: Classification) -> Self {
        self.classification = Some(c);
        self
    }

    /// Set the COSPAR (international) designator.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{InternationalDesignator, TleBuilder};
    /// let _b = TleBuilder::new()
    ///     .international_designator(InternationalDesignator("98067A".into()));
    /// ```
    pub fn international_designator(mut self, d: InternationalDesignator) -> Self {
        self.international_designator = Some(d);
        self
    }

    /// Set the epoch (UTC).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use tempoch::{Time, UTC};
    /// use chrono::Utc;
    /// let _b = TleBuilder::new().epoch(Time::<UTC>::from_chrono(Utc::now()));
    /// ```
    pub fn epoch(mut self, t: Time<UTC>) -> Self {
        self.epoch = Some(t);
        self
    }

    /// Set the first time derivative of mean motion (revolutions per day²).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().mean_motion_dot(-2.182e-5);
    /// ```
    pub fn mean_motion_dot(mut self, v: f64) -> Self {
        self.mean_motion_dot = v;
        self
    }

    /// Set the second time derivative of mean motion (revolutions per day³).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().mean_motion_ddot(0.0);
    /// ```
    pub fn mean_motion_ddot(mut self, v: f64) -> Self {
        self.mean_motion_ddot = v;
        self
    }

    /// Set the SGP4 BSTAR drag term (1 / Earth radii).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().bstar(-0.11606e-4);
    /// ```
    pub fn bstar(mut self, v: f64) -> Self {
        self.bstar = v;
        self
    }

    /// Set the element-set number.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().element_set_number(292);
    /// ```
    pub fn element_set_number(mut self, v: u16) -> Self {
        self.element_set_number = v;
        self
    }

    /// Set the revolution number at epoch.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().revolution_number_at_epoch(56_353);
    /// ```
    pub fn revolution_number_at_epoch(mut self, v: u32) -> Self {
        self.revolution_number_at_epoch = v;
        self
    }

    /// Set the inclination.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use qtty::angular::Degrees;
    /// let _b = TleBuilder::new().inclination(Degrees::new(51.6416));
    /// ```
    pub fn inclination(mut self, v: Degrees) -> Self {
        self.inclination = Some(v);
        self
    }

    /// Set the right ascension of the ascending node.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use qtty::angular::Degrees;
    /// let _b = TleBuilder::new().raan(Degrees::new(247.4627));
    /// ```
    pub fn raan(mut self, v: Degrees) -> Self {
        self.raan = Some(v);
        self
    }

    /// Set the eccentricity (dimensionless, 0 ≤ e < 1).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// let _b = TleBuilder::new().eccentricity(0.0006703);
    /// ```
    pub fn eccentricity(mut self, v: f64) -> Self {
        self.eccentricity = Some(v);
        self
    }

    /// Set the argument of perigee.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use qtty::angular::Degrees;
    /// let _b = TleBuilder::new().argument_of_perigee(Degrees::new(130.5360));
    /// ```
    pub fn argument_of_perigee(mut self, v: Degrees) -> Self {
        self.argument_of_perigee = Some(v);
        self
    }

    /// Set the mean anomaly.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use qtty::angular::Degrees;
    /// let _b = TleBuilder::new().mean_anomaly(Degrees::new(325.0288));
    /// ```
    pub fn mean_anomaly(mut self, v: Degrees) -> Self {
        self.mean_anomaly = Some(v);
        self
    }

    /// Set the typed mean motion (revolutions per day).
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::TleBuilder;
    /// use qtty::angular::Turn;
    /// use qtty::angular_rate::AngularRate;
    /// use qtty::time::Day;
    /// let _b = TleBuilder::new()
    ///     .mean_motion(AngularRate::<Turn, Day>::new(15.72125391));
    /// ```
    pub fn mean_motion(mut self, v: AngularRate<Turn, Day>) -> Self {
        self.mean_motion = Some(v);
        self
    }

    /// Finalize the builder into a [`TLE`].
    ///
    /// Returns [`TleError::BuilderMissingField`] if any mandatory field
    /// was left unset.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::formats::tle::{TleBuilder, TleError};
    /// assert!(matches!(TleBuilder::new().build(), Err(TleError::BuilderMissingField(_))));
    /// ```
    pub fn build(self) -> Result<TLE, TleError> {
        Ok(TLE {
            name: self.name,
            norad_id: self
                .norad_id
                .ok_or(TleError::BuilderMissingField("norad_id"))?,
            classification: self.classification.unwrap_or(Classification::Unclassified),
            international_designator: self
                .international_designator
                .ok_or(TleError::BuilderMissingField("international_designator"))?,
            epoch: self.epoch.ok_or(TleError::BuilderMissingField("epoch"))?,
            mean_motion_dot: self.mean_motion_dot,
            mean_motion_ddot: self.mean_motion_ddot,
            bstar: self.bstar,
            element_set_number: self.element_set_number,
            revolution_number_at_epoch: self.revolution_number_at_epoch,
            inclination: self
                .inclination
                .ok_or(TleError::BuilderMissingField("inclination"))?,
            raan: self.raan.ok_or(TleError::BuilderMissingField("raan"))?,
            eccentricity: self
                .eccentricity
                .ok_or(TleError::BuilderMissingField("eccentricity"))?,
            argument_of_perigee: self
                .argument_of_perigee
                .ok_or(TleError::BuilderMissingField("argument_of_perigee"))?,
            mean_anomaly: self
                .mean_anomaly
                .ok_or(TleError::BuilderMissingField("mean_anomaly"))?,
            mean_motion: self
                .mean_motion
                .ok_or(TleError::BuilderMissingField("mean_motion"))?,
        })
    }
}

/// Render a [`TLE`] back to its canonical 2-line ASCII representation.
///
/// The returned tuple is `(line1, line2)`; both are 69 ASCII characters
/// long with a valid checksum digit. Composing them with the satellite
/// `name` recovers the 3LE form.
///
/// # Examples
///
/// ```
/// use siderust::formats::tle::{format_tle, parse_tle};
/// let l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
/// let l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
/// let tle = parse_tle(l1, l2).unwrap();
/// let (rl1, rl2) = format_tle(&tle).unwrap();
/// assert_eq!(rl1, l1);
/// assert_eq!(rl2, l2);
/// ```
pub fn format_tle(tle: &TLE) -> Result<(String, String), TleError> {
    let cat = tle.norad_id.format_alpha5()?;
    let cls = tle.classification.as_char();
    let intl = format!("{:<8}", tle.international_designator.0);
    let dt = tle
        .epoch
        .try_to_chrono()
        .map_err(|e| TleError::EpochConversion(format!("{e:?}")))?;
    let year = dt.year();
    let yy = ((year % 100) + 100) % 100;
    let doy = dt.ordinal();
    let secs_of_day =
        dt.num_seconds_from_midnight() as f64 + dt.timestamp_subsec_nanos() as f64 * 1e-9;
    let frac_day = secs_of_day / 86_400.0;
    let day_full = doy as f64 + frac_day;
    let n_dot = format_signed_decimal(tle.mean_motion_dot);
    let n_ddot = format_assumed_decimal_exponent(tle.mean_motion_ddot);
    let bstar_s = format_assumed_decimal_exponent(tle.bstar);
    let elset = format!("{:>4}", tle.element_set_number);
    let mut l1 =
        format!("1 {cat}{cls} {intl} {yy:02}{day_full:012.8} {n_dot} {n_ddot} {bstar_s} 0 {elset}");
    debug_assert_eq!(
        l1.len(),
        68,
        "rendered TLE line1 prefix not 68 chars: {l1:?}"
    );
    let cs1 = compute_tle_checksum(&l1);
    l1.push_str(&cs1.to_string());

    let inc = format!("{:8.4}", tle.inclination.value());
    let raan = format!("{:8.4}", tle.raan.value());
    let ecc_scaled = (tle.eccentricity * 1.0e7).round() as u64;
    let ecc = format!("{ecc_scaled:07}");
    let argp = format!("{:8.4}", tle.argument_of_perigee.value());
    let m_ano = format!("{:8.4}", tle.mean_anomaly.value());
    let n_str = format!("{:11.8}", tle.mean_motion.value());
    let rev = format!("{:>5}", tle.revolution_number_at_epoch);
    let mut l2 = format!("2 {cat} {inc} {raan} {ecc} {argp} {m_ano} {n_str}{rev}");
    debug_assert_eq!(
        l2.len(),
        68,
        "rendered TLE line2 prefix not 68 chars: {l2:?}"
    );
    let cs2 = compute_tle_checksum(&l2);
    l2.push_str(&cs2.to_string());
    Ok((l1, l2))
}

fn format_signed_decimal(v: f64) -> String {
    let sign = if v.is_sign_negative() { '-' } else { ' ' };
    let mag = v.abs();
    let fract = (mag * 1.0e8).round() as u64;
    format!("{sign}.{fract:08}")
}

fn format_assumed_decimal_exponent(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return " 00000-0".to_string();
    }
    let sign = if v.is_sign_negative() { '-' } else { ' ' };
    let mag = v.abs();
    let mut exp = mag.log10().floor() as i32 + 1;
    let mut mantissa = mag / 10f64.powi(exp);
    let scaled = (mantissa * 1.0e5).round() as u64;
    if scaled >= 1_000_000 {
        mantissa /= 10.0;
        exp += 1;
    }
    let scaled = (mantissa * 1.0e5).round() as u64;
    let exp_sign = if exp < 0 { '-' } else { '+' };
    let exp_abs = exp.unsigned_abs();
    if exp_sign == '+' {
        format!("{sign}{scaled:05}+{exp_abs}")
    } else {
        format!("{sign}{scaled:05}-{exp_abs}")
    }
}
