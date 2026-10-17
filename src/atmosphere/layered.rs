// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Finite plane-parallel columns over Optica coefficients.
//!
//! Heights are absolute meters in one caller-declared vertical datum (MSL for
//! terrestrial site profiles); no geoid conversion, curvature or refraction is
//! implied. Ground and observer coincide at the first edge. A column does not
//! determine its vertical distribution: profile construction requires fractions.

use super::{mie_optical_depth, rayleigh_optical_depth_bodhaine99, AtmosphereProfile};
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames::Horizontal};
use crate::qtty::{unit::Meter, Albedos, Kilometers, Meters, Nanometers, OpticalDepths, Quantity};
use alloc::vec::Vec;
use optica::medium::{InverseLength, OpticalCoefficients, TryMedium};

/// Errors in an atmospheric optical column or request.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum LayeredAtmosphereError {
    /// Invalid finite, positive spectral interval or unsupported wavelength.
    #[error("wavelength must be finite, positive and inside the declared interval")]
    Wavelength,
    /// Nonfinite, reversed, out-of-domain or noncontiguous height edges.
    #[error("heights must be finite, ordered and inside the contiguous column")]
    Height,
    /// Incorrect number of layers or fractions.
    #[error("layer and edge dimensions do not agree")]
    Shape,
    /// Invalid or inconsistent Optica coefficients, including overflow.
    #[error("optical coefficients must be finite, nonnegative and internally consistent")]
    Coefficients,
    /// Invalid single-scattering albedo.
    #[error("single-scattering albedo must be finite and in [0, 1]")]
    Albedo,
    /// Fractions must be nonnegative and sum to one within rounding allowance.
    #[error("layer fractions must be finite, nonnegative and sum to one")]
    Fractions,
    /// Invalid site profile or inconsistent observer height.
    #[error("profile parameters or observer height are invalid")]
    Profile,
    /// The position has a different bound center, or nonfinite local coordinates.
    #[error("position must use the bound center and finite local coordinates")]
    Position,
}

type Result<T> = core::result::Result<T, LayeredAtmosphereError>;

/// Atmospheric constituent grouping, reusing Optica's coefficient type.
///
/// Values are inverse meters. Rayleigh has no absorption and molecular
/// absorption has no scattering. Private fields preserve validated invariants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphericOptics {
    rayleigh: OpticalCoefficients<Meter>,
    aerosol: OpticalCoefficients<Meter>,
    molecular: OpticalCoefficients<Meter>,
}

impl AtmosphericOptics {
    /// Validate complete constituent states, including public Optica fields.
    ///
    /// Rejects nonfinite/negative values, inconsistent extinction/SSA, arithmetic
    /// overflow, Rayleigh absorption and molecular scattering. Vacuum's stored
    /// SSA=0 is only Optica's computational sentinel; its physical ratio is undefined.
    pub fn try_new(
        rayleigh: OpticalCoefficients<Meter>,
        aerosol: OpticalCoefficients<Meter>,
        molecular: OpticalCoefficients<Meter>,
    ) -> Result<Self> {
        for c in [&rayleigh, &aerosol, &molecular] {
            validate_coefficients(c)?;
        }
        if rayleigh.sigma_a.value() != 0.0 || molecular.sigma_s.value() != 0.0 {
            return Err(LayeredAtmosphereError::Coefficients);
        }
        let optics = Self {
            rayleigh,
            aerosol,
            molecular,
        };
        optics.total()?;
        Ok(optics)
    }

    /// Construct from Rayleigh scattering, aerosol extinction/SSA and molecular
    /// absorption. All three coefficients are typed inverse meters.
    pub fn from_extinction(
        rayleigh: Quantity<InverseLength<Meter>>,
        aerosol_extinction: Quantity<InverseLength<Meter>>,
        aerosol_ssa: Albedos,
        molecular_absorption: Quantity<InverseLength<Meter>>,
    ) -> Result<Self> {
        let omega = aerosol_ssa.value();
        if !omega.is_finite() || !(0.0..=1.0).contains(&omega) {
            return Err(LayeredAtmosphereError::Albedo);
        }
        // Validate extinction before multiplication: a negative tiny value must
        // not disappear into signed zero when multiplied by an endpoint SSA.
        nonnegative(aerosol_extinction.value())?;
        Self::try_new(
            coefficients(0.0, rayleigh.value())?,
            coefficients(
                (1.0 - omega) * aerosol_extinction.value(),
                omega * aerosol_extinction.value(),
            )?,
            coefficients(molecular_absorption.value(), 0.0)?,
        )
    }

    /// Rayleigh coefficients (absorption exactly zero).
    pub fn rayleigh(&self) -> OpticalCoefficients<Meter> {
        self.rayleigh
    }
    /// Aerosol scattering and absorption, not extinction alone.
    pub fn aerosol(&self) -> OpticalCoefficients<Meter> {
        self.aerosol
    }
    /// Molecular absorption coefficients (scattering exactly zero).
    pub fn molecular(&self) -> OpticalCoefficients<Meter> {
        self.molecular
    }

    /// Sum constituents, checking arithmetic rather than trusting typed inputs.
    pub fn total(&self) -> Result<OpticalCoefficients<Meter>> {
        coefficients(
            self.aerosol.sigma_a.value() + self.molecular.sigma_a.value(),
            self.rayleigh.sigma_s.value() + self.aerosol.sigma_s.value(),
        )
    }
}

/// Constituent vertical optical depths, including separate scattering/absorption.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphericOpticalDepth {
    /// Rayleigh scattering depth.
    pub rayleigh: OpticalDepths,
    /// Aerosol scattering depth.
    pub aerosol_scattering: OpticalDepths,
    /// Aerosol absorption depth.
    pub aerosol_absorption: OpticalDepths,
    /// Molecular absorption depth.
    pub molecular_absorption: OpticalDepths,
}

impl AtmosphericOpticalDepth {
    /// Total scattering optical depth.
    pub fn scattering(&self) -> OpticalDepths {
        self.rayleigh + self.aerosol_scattering
    }
    /// Total absorption optical depth.
    pub fn absorption(&self) -> OpticalDepths {
        self.aerosol_absorption + self.molecular_absorption
    }
    /// Total extinction optical depth.
    pub fn extinction(&self) -> OpticalDepths {
        self.scattering() + self.absorption()
    }
}

#[derive(Debug, Clone)]
enum LayerModel {
    Constant(Vec<AtmosphericOptics>),
    Profile {
        profile: AtmosphereProfile,
        rayleigh_weights: Vec<f64>,
        aerosol_weights: Vec<f64>,
        ssa: Albedos,
    },
}

/// Validated, contiguous, piecewise-constant atmospheric layers.
///
/// No allocation occurs during queries. Wavelength-independent literal layers
/// or exact Bodhaine/Patat profile columns with explicit vertical fractions are
/// supported. The configured domain is a caller assertion, not site calibration.
/// Only the profile adapter restricts its optical approximation to 360–570 nm.
///
/// This is a fallible **vertical** API, not an infallible Optica `Medium`.
/// [`Self::at_origin`] binds an explicit center and absolute origin height for
/// Optica `TryMedium` queries in the local Horizontal frame. No geoid conversion
/// or mapping from ellipsoidal to MSL heights is guessed.
///
/// ```
/// use siderust::atmosphere::{AtmosphericOptics, LayeredAtmosphere};
/// use siderust::optica::medium::InverseLength;
/// use siderust::qtty::{unit::Meter, Albedos, Meters, Nanometers, Quantity};
/// let c = |v| Quantity::<InverseLength<Meter>>::new(v);
/// let optics = AtmosphericOptics::from_extinction(c(0.001), c(0.002),
///     Albedos::new(0.75), c(0.0))?;
/// let column = LayeredAtmosphere::try_new(&[Meters::new(1000.), Meters::new(2000.)],
///     [Nanometers::new(360.), Nanometers::new(570.)], vec![optics])?;
/// let local = column.coefficients(Nanometers::new(500.), Meters::new(1500.))?;
/// assert!((local.total()?.sigma_t.value() - 0.003).abs() < 1e-15);
/// let tau = column.vertical_optical_depth(Nanometers::new(500.),
///     Meters::new(1250.), Meters::new(1750.))?;
/// assert!((tau.extinction().value() - 1.5).abs() < 1e-14);
/// # Ok::<(), siderust::atmosphere::LayeredAtmosphereError>(())
/// ```
#[derive(Debug, Clone)]
pub struct LayeredAtmosphere {
    edges: Vec<Meters>,
    wavelengths: [Nanometers; 2],
    model: LayerModel,
}

impl LayeredAtmosphere {
    /// Construct wavelength-independent layers over an explicitly declared
    /// wavelength interval. A single shared edge vector excludes gaps/overlaps.
    /// All values are owned snapshots; later mutation of caller inputs is harmless.
    pub fn try_new(
        edges: &[Meters],
        wavelengths: [Nanometers; 2],
        optics: Vec<AtmosphericOptics>,
    ) -> Result<Self> {
        validate_domain(edges, wavelengths)?;
        if optics.len() + 1 != edges.len() {
            return Err(LayeredAtmosphereError::Shape);
        }
        Ok(Self {
            edges: edges.to_vec(),
            wavelengths,
            model: LayerModel::Constant(optics),
        })
    }

    /// Distribute an existing Rayleigh+aerosol profile using independent explicit
    /// fractions, and an explicit wavelength-independent aerosol SSA.
    ///
    /// Fractions for each constituent must sum to one within
    /// `8*n*f64::EPSILON`; only rounding within that allowance is normalized.
    /// Local pressure is passed to Bodhaine with **zero additional altitude**.
    /// The first edge must equal the profile observer altitude in meters (MSL).
    /// The original profile is retained for inspection; no site defaults change.
    /// Molecular absorption is zero here; literal layers can supply it explicitly.
    pub fn from_profile(
        profile: AtmosphereProfile,
        edges: &[Meters],
        wavelengths: [Nanometers; 2],
        rayleigh_fractions: &[f64],
        aerosol_fractions: &[f64],
        aerosol_ssa: Albedos,
    ) -> Result<Self> {
        validate_domain(edges, wavelengths)?;
        let p = profile;
        if !p.surface_pressure.value().is_finite()
            || p.surface_pressure.value() < 0.0
            || !p.rayleigh_scale_height.value().is_finite()
            || p.rayleigh_scale_height.value() <= 0.0
            || !p.observer_altitude.value().is_finite()
            || p.observer_altitude.to::<Meter>() != edges[0]
            || !p.mie_params.tau0.value().is_finite()
            || p.mie_params.tau0.value() < 0.0
            || !p.mie_params.alpha.is_finite()
            || !p.mie_params.lambda_ref.value().is_finite()
            || p.mie_params.lambda_ref.value() <= 0.0
            || wavelengths[0].value() < 360.0
            || wavelengths[1].value() > 570.0
        {
            return Err(LayeredAtmosphereError::Profile);
        }
        // The same constructor checks SSA even for a vacuum profile.
        AtmosphericOptics::from_extinction(
            Quantity::new(0.0),
            Quantity::new(0.0),
            aerosol_ssa,
            Quantity::new(0.0),
        )?;
        let rayleigh_weights = weights(edges, rayleigh_fractions)?;
        let aerosol_weights = weights(edges, aerosol_fractions)?;
        let column = Self {
            edges: edges.to_vec(),
            wavelengths,
            model: LayerModel::Profile {
                profile,
                rayleigh_weights,
                aerosol_weights,
                ssa: aerosol_ssa,
            },
        };
        for wavelength in wavelengths {
            for i in 0..edges.len() - 1 {
                column.layer(wavelength, i)?;
            }
        }
        Ok(column)
    }

    /// Bind a local Horizontal coordinate origin to an explicitly supplied
    /// absolute height in this column's datum. The origin can be outside the
    /// column; each subsequent position query must land inside it.
    pub fn at_origin<C: ReferenceCenter>(
        &self,
        center: C::Params,
        height: Meters,
    ) -> Result<AtmosphericColumnMedium<'_, C>> {
        if !height.value().is_finite() {
            return Err(LayeredAtmosphereError::Height);
        }
        Ok(AtmosphericColumnMedium {
            column: self,
            center,
            height,
        })
    }

    /// Absolute contiguous layer boundaries, including ground/observer and top.
    pub fn edges(&self) -> &[Meters] {
        &self.edges
    }
    /// Inclusive supported wavelength interval.
    pub fn wavelength_domain(&self) -> [Nanometers; 2] {
        self.wavelengths
    }
    /// Original column model, when constructed from an atmospheric profile.
    pub fn profile(&self) -> Option<&AtmosphereProfile> {
        match &self.model {
            LayerModel::Profile { profile, .. } => Some(profile),
            _ => None,
        }
    }

    /// Evaluate at absolute height. Internal edges belong to the layer above;
    /// the final top is included in the last layer. No extrapolation/clamping.
    pub fn coefficients(
        &self,
        wavelength: Nanometers,
        height: Meters,
    ) -> Result<AtmosphericOptics> {
        self.check_wavelength(wavelength)?;
        self.check_height(height)?;
        let i = self.edges[1..self.edges.len() - 1].partition_point(|edge| *edge <= height);
        self.layer(wavelength, i)
    }

    /// Exact partial-layer vertical integral in deterministic bottom-to-top order.
    /// Equal bounds return zero, including at the top; reversed/outside bounds fail.
    pub fn vertical_optical_depth(
        &self,
        wavelength: Nanometers,
        lower: Meters,
        upper: Meters,
    ) -> Result<AtmosphericOpticalDepth> {
        self.check_wavelength(wavelength)?;
        self.check_height(lower)?;
        self.check_height(upper)?;
        if lower > upper {
            return Err(LayeredAtmosphereError::Height);
        }
        let mut sums = [0.0; 4];
        for (i, bounds) in self.edges.windows(2).enumerate() {
            let thickness =
                upper.value().min(bounds[1].value()) - lower.value().max(bounds[0].value());
            if thickness <= 0.0 {
                continue;
            }
            let c = self.layer(wavelength, i)?;
            for (sum, beta) in sums.iter_mut().zip([
                c.rayleigh.sigma_s.value(),
                c.aerosol.sigma_s.value(),
                c.aerosol.sigma_a.value(),
                c.molecular.sigma_a.value(),
            ]) {
                *sum += beta * thickness;
                nonnegative(*sum)?;
            }
        }
        let out = AtmosphericOpticalDepth {
            rayleigh: OpticalDepths::new(sums[0]),
            aerosol_scattering: OpticalDepths::new(sums[1]),
            aerosol_absorption: OpticalDepths::new(sums[2]),
            molecular_absorption: OpticalDepths::new(sums[3]),
        };
        nonnegative(out.extinction().value())?;
        Ok(out)
    }

    fn layer(&self, wavelength: Nanometers, i: usize) -> Result<AtmosphericOptics> {
        match &self.model {
            LayerModel::Constant(values) => Ok(values[i]),
            LayerModel::Profile {
                profile: p,
                rayleigh_weights,
                aerosol_weights,
                ssa,
            } => {
                let r = rayleigh_optical_depth_bodhaine99(
                    wavelength,
                    p.surface_pressure,
                    Kilometers::new(0.0),
                    p.rayleigh_scale_height,
                )
                .value();
                let a = mie_optical_depth(&p.mie_params, wavelength).value();
                AtmosphericOptics::from_extinction(
                    Quantity::new(r * rayleigh_weights[i]),
                    Quantity::new(a * aerosol_weights[i]),
                    *ssa,
                    Quantity::new(0.0),
                )
            }
        }
    }

    fn check_wavelength(&self, w: Nanometers) -> Result<()> {
        if !w.value().is_finite() || w < self.wavelengths[0] || w > self.wavelengths[1] {
            Err(LayeredAtmosphereError::Wavelength)
        } else {
            Ok(())
        }
    }
    fn check_height(&self, h: Meters) -> Result<()> {
        if !h.value().is_finite() || h < self.edges[0] || h > self.edges[self.edges.len() - 1] {
            Err(LayeredAtmosphereError::Height)
        } else {
            Ok(())
        }
    }
}

fn nonnegative(value: f64) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        Err(LayeredAtmosphereError::Coefficients)
    } else {
        Ok(())
    }
}
fn coefficients(absorption: f64, scattering: f64) -> Result<OpticalCoefficients<Meter>> {
    let c = OpticalCoefficients::try_new(absorption, scattering)
        .map_err(|_| LayeredAtmosphereError::Coefficients)?;
    nonnegative(c.sigma_t.value())?;
    Ok(c)
}
fn validate_coefficients(c: &OpticalCoefficients<Meter>) -> Result<()> {
    let expected = coefficients(c.sigma_a.value(), c.sigma_s.value())?;
    if c.sigma_t != expected.sigma_t || c.ssa != expected.ssa {
        return Err(LayeredAtmosphereError::Coefficients);
    }
    Ok(())
}
fn validate_domain(edges: &[Meters], wavelengths: [Nanometers; 2]) -> Result<()> {
    if edges.len() < 2 {
        return Err(LayeredAtmosphereError::Shape);
    }
    if edges.iter().any(|h| !h.value().is_finite())
        || edges
            .windows(2)
            .any(|w| w[0] >= w[1] || !(w[1] - w[0]).value().is_finite())
    {
        return Err(LayeredAtmosphereError::Height);
    }
    if !wavelengths[0].value().is_finite()
        || !wavelengths[1].value().is_finite()
        || wavelengths[0].value() <= 0.0
        || wavelengths[0] >= wavelengths[1]
    {
        return Err(LayeredAtmosphereError::Wavelength);
    }
    Ok(())
}
fn weights(edges: &[Meters], fractions: &[f64]) -> Result<Vec<f64>> {
    if fractions.len() + 1 != edges.len() {
        return Err(LayeredAtmosphereError::Shape);
    }
    let sum: f64 = fractions.iter().sum();
    if fractions.iter().any(|v| !v.is_finite() || *v < 0.0)
        || !sum.is_finite()
        || (sum - 1.0).abs() > 8.0 * fractions.len() as f64 * f64::EPSILON
    {
        return Err(LayeredAtmosphereError::Fractions);
    }
    fractions
        .iter()
        .zip(edges.windows(2))
        .map(|(f, b)| {
            let w = (f / sum) / (b[1] - b[0]).value();
            nonnegative(w)?;
            Ok(w)
        })
        .collect()
}

/// Borrowed fallible Optica medium with an explicit local-to-absolute height map.
///
/// The local Horizontal frame's z axis points upward. Its position center must
/// match the bound center parameters. Latitude/longitude and a geoid model are
/// not used to infer MSL from an ellipsoidal height. Querying ignores horizontal
/// location only after requiring finite coordinates (horizontal homogeneity).
#[derive(Debug, Clone)]
pub struct AtmosphericColumnMedium<'a, C: ReferenceCenter> {
    column: &'a LayeredAtmosphere,
    center: C::Params,
    height: Meters,
}

impl<C: ReferenceCenter> TryMedium<C, Horizontal, Meter> for AtmosphericColumnMedium<'_, C> {
    type Error = LayeredAtmosphereError;

    fn try_coefficients(
        &self,
        p: Position<C, Horizontal, Meter>,
        wavelength: Nanometers,
    ) -> Result<OpticalCoefficients<Meter>> {
        if p.center_params() != &self.center
            || [p.x(), p.y(), p.z()].iter().any(|x| !x.value().is_finite())
        {
            return Err(LayeredAtmosphereError::Position);
        }
        self.column
            .coefficients(wavelength, self.height + p.z())?
            .total()
    }
}
