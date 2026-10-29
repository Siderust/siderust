// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Explicit geometric lunar observations with runtime SPK and binary PCK data.
use super::{SpiceContext, SpiceContextError};
use crate::bodies::solar_system::Moon;
use crate::coordinates::surface::{LatitudeConvention, ReferenceSurface, SurfaceSite};
use crate::coordinates::transform::HorizontalObservation;
use crate::formats::spice::{FrameClass, SpiceError};
use crate::qtty::{Degrees, Kilometer, Meters, Quantity};
use tempoch::{Time, TDB};

/// A target centre supported by the lunar geometric observer.
#[derive(Debug, Clone, Copy)]
pub enum LunarTarget {
    /// Sun centre, NAIF 10.
    Sun,
    /// Earth centre, NAIF 399.
    Earth,
}

/// A site in a named runtime lunar frame, bound to one immutable kernel set.
///
/// The caller selects mutually compatible SPK/BPC/FK kernels; filenames are
/// not evidence of compatibility. The frame chain must be Moon-centred and
/// terminate in a binary PCK with J2000 base. No analytical fallback occurs.
/// Site coordinates are interpreted in `frame`, never in analytical IAU_MOON.
/// Reuses the SPK evaluator underlying Siderust's runtime ephemerides directly,
/// avoiding TT round trips and ecliptic rotations for J2000 SPICE states.
pub struct LunarObserver<'a> {
    context: &'a SpiceContext,
    frame: i32,
    site: SurfaceSite<Moon>,
}

/// Epoch- and frame-tagged geometric result (`abcorr = NONE`).
pub struct LunarObservation<'a> {
    /// Local north/east/up observation, with finite-distance parallax.
    pub horizontal: HorizontalObservation<Moon, Kilometer>,
    /// Evaluation epoch in TDB (SPICE ET).
    pub epoch: Time<TDB>,
    /// NAIF frame ID in which site coordinates were interpreted.
    pub frame_id: i32,
    /// Target centre.
    pub target: LunarTarget,
    /// Loaded kernels used by the immutable observation context.
    pub kernels: Vec<&'a str>,
    /// Moon-centred target position in the selected fixed frame (km).
    pub target_fixed: [Quantity<Kilometer>; 3],
    /// Site-relative vector in the selected fixed frame (km).
    pub relative_fixed: [Quantity<Kilometer>; 3],
}

impl LunarObservation<'_> {
    /// Whether the target's centre is strictly above the geometric horizon.
    pub fn above_horizon(&self) -> bool {
        self.horizontal.alt().value() > 0.0
    }
}

impl<'a> LunarObserver<'a> {
    /// Bind a validated surface site to an explicit FK frame and kernel set.
    pub fn new(
        context: &'a SpiceContext,
        frame: &str,
        site: SurfaceSite<Moon>,
    ) -> Result<Self, SpiceContextError> {
        let id = context.frame_registry().frame_by_name(frame)?.id;
        let mut current = id;
        for _ in 0..64 {
            let spec = context.frame_registry().frame_by_id(current)?;
            if spec.center_id != 301 {
                return Err(invalid("lunar observer frame must be Moon-centred"));
            }
            if let Some((_, relative)) = &spec.tk_spec {
                current = context.frame_registry().frame_by_name(relative)?.id;
            } else if spec.class == FrameClass::Pck {
                if !context
                    .kernel_set()
                    .binary_pcks()
                    .any(|(_, b)| b.contains_class(spec.class_id))
                {
                    return Err(invalid(
                        "lunar observer requires binary PCK data for selected frame",
                    ));
                }
                if context.kernel_set().binary_pcks().any(|(_, b)| {
                    b.contains_class(spec.class_id) && !b.has_j2000_base(spec.class_id)
                }) {
                    return Err(invalid(
                        "lunar observer supports J2000-based binary PCK only",
                    ));
                }
                return Ok(Self {
                    context,
                    frame: id,
                    site,
                });
            } else {
                return Err(invalid("lunar frame chain must end in binary PCK"));
            }
        }
        Err(invalid("cyclic or excessively deep lunar frame graph"))
    }

    /// Define a radial site directly by distance from the Moon centre.
    /// Longitude is east-positive, latitude planetocentric; local up is radial.
    /// No ellipsoidal height or terrain normal is inferred.
    pub fn radial(
        context: &'a SpiceContext,
        frame: &str,
        lon: Degrees,
        lat: Degrees,
        radius: Meters,
    ) -> Result<Self, SpiceContextError> {
        let site = SurfaceSite::with_surface(
            lon,
            lat,
            Meters::new(0.0),
            ReferenceSurface::sphere(radius)?,
            LatitudeConvention::Planetocentric,
        )?;
        Self::new(context, frame, site)
    }

    /// Site geometry interpreted in this observer's runtime lunar frame.
    pub fn site(&self) -> &SurfaceSite<Moon> {
        &self.site
    }

    /// Compute geometric centre directions at a typed TDB epoch.
    /// Convert UTC/TT explicitly with `tempoch` before calling. Both state and
    /// orientation use the same ET seconds; no light time or aberration applies.
    pub fn observe(
        &self,
        target: LunarTarget,
        epoch: Time<TDB>,
    ) -> Result<LunarObservation<'a>, SpiceContextError> {
        let et = epoch.to_j2000s().raw().value();
        let matrix = self.context.rotation_naif(1, self.frame, et)?;
        let state = self.context.state_naif(
            match target {
                LunarTarget::Sun => 10,
                LunarTarget::Earth => 399,
            },
            301,
            et,
        )?;
        let target_fixed = apply(matrix, [state[0], state[1], state[2]]);
        let observer = self
            .site
            .body_fixed_cartesian()?
            .map(|q| q.to::<Kilometer>().value());
        let relative = core::array::from_fn(|i| target_fixed[i] - observer[i]);
        let horizontal = HorizontalObservation::from_fixed_relative(
            &self.site,
            relative.map(Quantity::<Kilometer>::new),
        )?;
        Ok(LunarObservation {
            horizontal,
            epoch,
            frame_id: self.frame,
            target,
            kernels: self.context.kernel_set().provenance(),
            target_fixed: target_fixed.map(Quantity::new),
            relative_fixed: relative.map(Quantity::new),
        })
    }
}

fn invalid(message: &str) -> SpiceContextError {
    SpiceError::Corrupted {
        message: message.into(),
    }
    .into()
}

fn apply(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    m.map(|r| r[0] * v[0] + r[1] * v[1] + r[2] * v[2])
}
