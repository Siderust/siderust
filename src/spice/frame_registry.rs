// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame registry built from FK kernels.

use crate::formats::spice::{FrameClass, FrameKernel, FrameSpec, SpiceError};

/// Registry of known frames, populated from FK kernels.
///
/// Includes a small built-in table of well-known inertial frames.
#[derive(Debug)]
pub struct FrameRegistry {
    frames: Vec<FrameSpec>,
}

impl Default for FrameRegistry {
    fn default() -> Self {
        Self {
            frames: vec![FrameSpec {
                id: 1,
                name: "J2000".to_string(),
                class: FrameClass::Inertial,
                class_id: 1,
                center_id: 0,
                tk_spec: None,
            }],
        }
    }
}

impl FrameRegistry {
    /// Create an empty registry seeded with built-in frames.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add all frames from an FK kernel.
    pub fn add_fk(&mut self, fk: &FrameKernel) {
        self.frames.extend(fk.frames.iter().cloned());
    }

    /// Look up a frame by integer ID.
    pub fn frame_by_id(&self, id: i32) -> Result<&FrameSpec, SpiceError> {
        self.frames
            .iter()
            .rev()
            .find(|frame| frame.id == id)
            .ok_or_else(|| SpiceError::UnknownFrame {
                description: format!("frame ID {id}"),
            })
    }

    /// Look up a frame by name (case-insensitive).
    pub fn frame_by_name(&self, name: &str) -> Result<&FrameSpec, SpiceError> {
        let upper = name.to_ascii_uppercase();
        self.frames
            .iter()
            .rev()
            .find(|frame| frame.name == upper)
            .ok_or_else(|| SpiceError::UnknownFrame {
                description: name.to_string(),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::spice::{FrameClass, FrameKernel};

    #[test]
    fn builtin_j2000_lookup_by_id_and_name() {
        let registry = FrameRegistry::new();
        assert_eq!(registry.frame_by_id(1).unwrap().name, "J2000");
        assert_eq!(registry.frame_by_name("j2000").unwrap().id, 1);
        assert_eq!(
            registry.frame_by_name("J2000").unwrap().class,
            FrameClass::Inertial
        );
    }

    #[test]
    fn unknown_frame_errors_are_structured() {
        let registry = FrameRegistry::new();
        assert!(matches!(
            registry.frame_by_id(9_999).unwrap_err(),
            SpiceError::UnknownFrame { description } if description.contains("9999")
        ));
        assert!(matches!(
            registry.frame_by_name("NOPE").unwrap_err(),
            SpiceError::UnknownFrame { description } if description == "NOPE"
        ));
    }

    #[test]
    fn add_fk_extends_lookup() {
        let src = "\\begindata\nFRAME_1234_NAME = 'TEST'\nFRAME_1234_CLASS = 1\nFRAME_1234_CLASS_ID = 1234\nFRAME_1234_CENTER = 0\n";
        let fk = FrameKernel::from_text(src).unwrap();
        let mut registry = FrameRegistry::new();
        registry.add_fk(&fk);
        assert_eq!(registry.frame_by_id(1234).unwrap().name, "TEST");
        assert_eq!(registry.frame_by_name("test").unwrap().id, 1234);
    }
}
