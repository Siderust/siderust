// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level SPICE kernel context.

use std::path::Path;

use affn::cartesian::{Position, Velocity};
use qtty::unit::Kilometer;
use qtty::{KmPerSecond, Quantity};
use tempoch::{Time, TDB};

use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::ICRS;
use crate::formats::spice::{
    fk::{FrameClass, TkSpec},
    CkKernel, FrameKernel, IkKernel, LeapSecondKernel, PckKernel, SclkKernel, SpiceError,
    SpkKernel, TextKernel,
};

use super::{FrameRegistry, KernelSet, SpiceContextError};

/// High-level SPICE kernel context.
///
/// `SpiceContext` owns a [`KernelSet`] and a [`FrameRegistry`]. It is the main
/// entry point for loading kernels and querying state, time offsets, and a
/// subset of frame rotations.
///
/// # Examples
///
/// ```rust,no_run
/// use siderust::spice::SpiceContext;
///
/// let mut ctx = SpiceContext::new();
/// ctx.load("de440.bsp")?;
/// let state = ctx.state_naif(399, 0, 0.0)?;
/// println!("Earth x = {} km", state[0]);
/// # Ok::<_, siderust::spice::SpiceContextError>(())
/// ```
pub struct SpiceContext {
    kernel_set: KernelSet,
    frame_registry: FrameRegistry,
}

impl SpiceContext {
    /// Create an empty context with no kernels loaded.
    pub fn new() -> Self {
        Self {
            kernel_set: KernelSet::new(),
            frame_registry: FrameRegistry::new(),
        }
    }

    /// Load a kernel from a filesystem path.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<(), SpiceContextError> {
        let path = path.as_ref();
        let alias = path.display().to_string();
        let bytes = std::fs::read(path)?;
        self.load_bytes(alias, bytes)
    }

    /// Load a kernel from in-memory bytes.
    pub fn load_bytes(
        &mut self,
        alias: impl Into<String>,
        bytes: Vec<u8>,
    ) -> Result<(), SpiceContextError> {
        let alias = alias.into();
        let locator = if bytes.len() >= 8 {
            std::str::from_utf8(&bytes[0..8])
                .unwrap_or("")
                .trim()
                .to_ascii_uppercase()
        } else {
            String::new()
        };

        if locator.starts_with("DAF/SPK") {
            let kernel = SpkKernel::from_bytes(bytes)?;
            self.kernel_set.add_spk(alias, kernel);
        } else if locator.starts_with("DAF/CK") {
            let kernel = CkKernel::from_bytes(bytes)?;
            self.kernel_set.add_ck(alias, kernel);
        } else if locator.starts_with("DAF/PCK") || locator.starts_with("DAF/BPC") {
            return Err(SpiceContextError::UnsupportedKernelQuery {
                message: format!(
                    "Binary PCK kernel '{alias}' loaded but BPC rotation evaluation is not implemented in V1; load a PCK text kernel instead"
                ),
            });
        } else if locator.starts_with("DAF/DSK") || locator.starts_with("DAF/EK") {
            return Err(SpiceContextError::UnsupportedKernelQuery {
                message: format!(
                    "DSK/EK kernel '{alias}': surface/table queries are not implemented in V1"
                ),
            });
        } else {
            let text = std::str::from_utf8(&bytes).map_err(|_| {
                SpiceContextError::Kernel(SpiceError::FormatParse(format!(
                    "'{alias}': not a recognised DAF binary kernel and not valid UTF-8 text"
                )))
            })?;
            self.load_text_kernel(alias, text)?;
        }
        Ok(())
    }

    fn load_text_kernel(&mut self, alias: String, text: &str) -> Result<(), SpiceContextError> {
        let tk = TextKernel::parse(text)?;
        let has_lsk = tk.get("DELTET/DELTA_T_A").is_some() || tk.get("DELTET/DELTA_AT").is_some();
        let has_fk = tk.data.keys().any(|key| key.starts_with("FRAME_"));
        let has_pck = tk.data.keys().any(|key| {
            key.starts_with("BODY") && (key.contains("_POLE_RA") || key.contains("_RADII"))
        });
        let has_sclk = tk
            .data
            .keys()
            .any(|key| key.starts_with("SCLK01_COEFFICIENTS_"));
        let has_ik = tk
            .data
            .keys()
            .any(|key| key.starts_with("INS") && key.contains("_FOV"));

        if has_lsk {
            self.kernel_set
                .set_lsk(alias.clone(), LeapSecondKernel::from_text(text)?);
        }
        if has_fk {
            let fk = FrameKernel::from_text(text)?;
            self.frame_registry.add_fk(&fk);
            self.kernel_set.add_fk(alias.clone(), fk);
        }
        if has_pck {
            self.kernel_set
                .add_pck_text(alias.clone(), PckKernel::from_text(text)?);
        }
        if has_sclk {
            for sclk in SclkKernel::from_text(text)? {
                self.kernel_set.add_sclk(alias.clone(), sclk);
            }
        }
        if has_ik {
            self.kernel_set.add_ik(alias, IkKernel::from_text(text)?);
        }
        Ok(())
    }

    /// Compute state `[x, y, z, vx, vy, vz]` in km and km/s.
    pub fn state_naif(
        &self,
        target: i32,
        center: i32,
        epoch_tdb_s: f64,
    ) -> Result<[f64; 6], SpiceContextError> {
        Ok(self.kernel_set.state_naif(target, center, epoch_tdb_s)?)
    }

    /// Rotation matrix (3×3, row-major) from frame `from_id` to frame `to_id`.
    pub fn rotation_naif(
        &self,
        from_id: i32,
        to_id: i32,
        epoch_tdb_s: f64,
    ) -> Result<[[f64; 3]; 3], SpiceContextError> {
        if from_id == to_id {
            return Ok(identity_3x3());
        }
        let from_to_j2000 = self.frame_to_j2000(from_id, epoch_tdb_s)?;
        let to_to_j2000 = self.frame_to_j2000(to_id, epoch_tdb_s)?;
        Ok(mat3_mul(mat3_transpose(to_to_j2000), from_to_j2000))
    }

    /// Convert TDB seconds past J2000 to the TAI−UTC offset at that epoch.
    pub fn tai_minus_utc(&self, epoch_tdb_s: f64) -> Result<f64, SpiceContextError> {
        let lsk = self
            .kernel_set
            .lsk()
            .ok_or_else(|| SpiceContextError::KernelNotLoaded {
                kernel_type: "LSK".to_string(),
            })?;
        Ok(lsk.tai_minus_utc_at(epoch_tdb_s))
    }

    /// Convert TDB seconds past J2000 to TDB − UTC at that epoch.
    pub fn tdb_minus_utc(&self, epoch_tdb_s: f64) -> Result<f64, SpiceContextError> {
        let lsk = self
            .kernel_set
            .lsk()
            .ok_or_else(|| SpiceContextError::KernelNotLoaded {
                kernel_type: "LSK".to_string(),
            })?;
        Ok(lsk.tdb_minus_utc_at(epoch_tdb_s))
    }

    /// Typed state query using `affn` position and velocity wrappers.
    #[allow(clippy::type_complexity)]
    pub fn state_typed<C>(
        &self,
        naif_body: i32,
        naif_center: i32,
        epoch: Time<TDB>,
    ) -> Result<(Position<C, ICRS, Kilometer>, Velocity<ICRS, KmPerSecond>), SpiceContextError>
    where
        C: ReferenceCenter<Params = ()>,
    {
        let tdb_s = epoch.to_j2000s().raw().value();
        let state = self.state_naif(naif_body, naif_center, tdb_s)?;
        let position = Position::<C, ICRS, Kilometer>::new(
            Quantity::<Kilometer>::new(state[0]),
            Quantity::<Kilometer>::new(state[1]),
            Quantity::<Kilometer>::new(state[2]),
        );
        let velocity = Velocity::<ICRS, KmPerSecond>::new(
            Quantity::<KmPerSecond>::new(state[3]),
            Quantity::<KmPerSecond>::new(state[4]),
            Quantity::<KmPerSecond>::new(state[5]),
        );
        Ok((position, velocity))
    }

    /// Borrow the underlying kernel set.
    pub fn kernel_set(&self) -> &KernelSet {
        &self.kernel_set
    }

    /// Borrow the frame registry.
    pub fn frame_registry(&self) -> &FrameRegistry {
        &self.frame_registry
    }

    fn frame_to_j2000(
        &self,
        frame_id: i32,
        epoch_tdb_s: f64,
    ) -> Result<[[f64; 3]; 3], SpiceContextError> {
        if frame_id == 1 {
            return Ok(identity_3x3());
        }

        let spec = self
            .frame_registry
            .frame_by_id(frame_id)
            .map_err(spice_error_to_context)?;

        if let Some((TkSpec::Matrix(matrix), relative_name)) = &spec.tk_spec {
            let relative = self
                .frame_registry
                .frame_by_name(relative_name)
                .map_err(spice_error_to_context)?;
            let relative_to_j2000 = self.frame_to_j2000(relative.id, epoch_tdb_s)?;
            return Ok(mat3_mul(relative_to_j2000, *matrix));
        }

        match spec.class {
            FrameClass::Inertial => Err(SpiceContextError::UnsupportedKernelQuery {
                message: format!("rotation from inertial frame {frame_id} to J2000 is not implemented"),
            }),
            FrameClass::Pck => self.iau_body_rotation(spec.class_id, epoch_tdb_s),
            FrameClass::Ck => Err(SpiceContextError::UnsupportedKernelQuery {
                message: format!("CK frame {frame_id} rotation requires CK/SCLK composition outside rotation_naif V1"),
            }),
            FrameClass::TwoVector | FrameClass::Dynamic => {
                Err(SpiceContextError::UnsupportedKernelQuery {
                    message: format!("frame {frame_id} is not implemented in V1"),
                })
            }
        }
    }

    fn iau_body_rotation(
        &self,
        body_id: i32,
        epoch_tdb_s: f64,
    ) -> Result<[[f64; 3]; 3], SpiceContextError> {
        for pck in self.kernel_set.pck_text_kernels().rev() {
            if let Some(orientation) = pck.body(body_id) {
                return Ok(orientation.rotation_to_j2000(epoch_tdb_s));
            }
        }
        Err(SpiceContextError::KernelNotLoaded {
            kernel_type: format!("PCK data for body {body_id}"),
        })
    }
}

impl Default for SpiceContext {
    fn default() -> Self {
        Self::new()
    }
}

fn spice_error_to_context(error: SpiceError) -> SpiceContextError {
    match error {
        SpiceError::UnknownFrame { description } => SpiceContextError::UnknownFrame { description },
        SpiceError::TimeConversion { message } => SpiceContextError::TimeConversion { message },
        SpiceError::UnsupportedKernelQuery { message } => {
            SpiceContextError::UnsupportedKernelQuery { message }
        }
        other => SpiceContextError::Kernel(other),
    }
}

fn identity_3x3() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn mat3_mul(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += left[i][k] * right[k][j];
            }
        }
    }
    out
}

fn mat3_transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]],
    ]
}
