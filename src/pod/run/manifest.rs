//! Run-level provenance.
//!
//! [`RunManifest`] captures everything needed to reproduce or audit a POD
//! run: tool version, configuration hash, hashed inputs, hashed outputs,
//! and the wall-clock interval. JSON serialisation is deterministic
//! (sorted keys, no incidental ordering) so two identical runs produce
//! byte-identical manifests.

use super::dataset::DatasetRef;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Reproducible record of a single POD run.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RunManifest {
    /// Free-form identifier (e.g. UUID) for the run.
    pub run_id: String,
    /// Tool version string (e.g. `"siderust-pod 0.0.0"`).
    pub tool_version: String,
    /// Hex SHA-256 of the canonical run configuration document.
    pub config_sha256: String,
    /// Inputs consumed by the run, sorted by `(kind, path)` for determinism.
    pub inputs: Vec<DatasetRef>,
    /// Artifacts produced by the run, sorted by `(kind, path)` for determinism.
    pub outputs: Vec<DatasetRef>,
    /// RFC-3339 UTC timestamp of run start.
    pub started_at: String,
    /// RFC-3339 UTC timestamp of run end.
    pub finished_at: String,
}

impl RunManifest {
    /// Sort inputs and outputs into the canonical order. Call before serialising.
    pub fn canonicalize(&mut self) {
        let key = |d: &DatasetRef| (d.kind.clone(), d.path.clone());
        self.inputs.sort_by_key(key);
        self.outputs.sort_by_key(&key);
    }

    /// Serialize to a deterministic pretty-printed JSON string.
    #[cfg(feature = "serde")]
    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        let mut clone = self.clone();
        clone.canonicalize();
        serde_json::to_string_pretty(&clone)
    }
}

/// Embedded JSON-Schema (draft-07) for a v1 [`RunManifest`].
///
/// The schema lives next to this source file under `schema/run_manifest.v1.json`
/// and is included verbatim at build time so the crate ships a single source
/// of truth without filesystem access at runtime.
#[cfg(feature = "serde")]
pub const RUN_MANIFEST_SCHEMA_V1: &str = include_str!("../schema/run_manifest.v1.json");

/// Embedded JSON-Schema (draft-07) for a v1 QC report (`qc.v1`).
#[cfg(feature = "serde")]
pub const QC_SCHEMA_V1: &str = include_str!("../schema/qc.v1.json");

#[cfg(all(test, feature = "serde"))]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn sample() -> RunManifest {
        RunManifest {
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            tool_version: "siderust-pod test".into(),
            config_sha256: "0".repeat(64),
            inputs: vec![
                DatasetRef {
                    path: PathBuf::from("/in/b.sp3"),
                    kind: "sp3".into(),
                    bytes: 1,
                    sha256: "1".repeat(64),
                },
                DatasetRef {
                    path: PathBuf::from("/in/a.sp3"),
                    kind: "sp3".into(),
                    bytes: 1,
                    sha256: "2".repeat(64),
                },
            ],
            outputs: vec![],
            started_at: "2026-05-11T22:00:00Z".into(),
            finished_at: "2026-05-11T22:00:01Z".into(),
        }
    }

    #[test]
    fn json_round_trip_is_deterministic() {
        let s1 = sample().to_json_pretty().unwrap();
        let s2 = sample().to_json_pretty().unwrap();
        assert_eq!(s1, s2);
        // Canonicalisation sorts inputs by (kind, path).
        let parsed: RunManifest = serde_json::from_str(&s1).unwrap();
        assert_eq!(parsed.inputs[0].path.to_str().unwrap(), "/in/a.sp3");
    }

    /// The embedded schema must parse as a JSON object and declare the keys
    /// the manifest serialises into; this prevents silent drift between
    /// `RunManifest` and `schema/run_manifest.v1.json`.
    #[test]
    fn embedded_schema_lists_all_required_fields() {
        let schema: serde_json::Value =
            serde_json::from_str(RUN_MANIFEST_SCHEMA_V1).expect("schema is JSON");
        let req = schema["required"].as_array().expect("required is array");
        let req: Vec<&str> = req.iter().map(|v| v.as_str().unwrap()).collect();
        for k in [
            "run_id",
            "tool_version",
            "config_sha256",
            "inputs",
            "outputs",
            "started_at",
            "finished_at",
        ] {
            assert!(req.contains(&k), "schema missing required field {k}");
        }
        let serialised = sample().to_json_pretty().unwrap();
        let v: serde_json::Value = serde_json::from_str(&serialised).unwrap();
        for k in &req {
            assert!(
                v.get(*k).is_some(),
                "manifest missing schema-required field {k}"
            );
        }
    }

    #[test]
    fn qc_schema_is_well_formed_v1() {
        let schema: serde_json::Value = serde_json::from_str(QC_SCHEMA_V1).unwrap();
        assert_eq!(
            schema["properties"]["schema"]["const"].as_str(),
            Some("qc.v1")
        );
        let req = schema["required"].as_array().unwrap();
        let req: Vec<&str> = req.iter().map(|v| v.as_str().unwrap()).collect();
        assert!(req.contains(&"schema"));
        assert!(req.contains(&"rtn"));
    }
}
