//! # Run manifest JSON writer
//!
//! ## Scientific scope
//!
//! Every POD run must be reproducible and auditable. The run manifest records
//! the complete provenance of a run — input files with content hashes, output
//! product paths, configuration hash, software version, and wall-clock
//! interval — so that any downstream consumer can verify exactly what was
//! computed and how.
//!
//! The scientific validity of outputs depends on correct provenance recording;
//! a manifest with wrong hashes cannot be distinguished from a forged result
//! at the file-system level.
//!
//! ## Technical scope
//!
//! [`ManifestWriter`] serialises a [`RunManifest`] to a deterministic,
//! pretty-printed JSON document and writes it to any [`std::io::Write`] sink.
//! "Deterministic" means that two calls with identical `RunManifest` values
//! produce byte-identical output regardless of the order in which `inputs` or
//! `outputs` were pushed onto the struct.
//!
//! The canonical order is obtained by calling
//! [`RunManifest::canonicalize`] before serialisation.
//!
//! ## References
//!
//! - Bray, T. (2017). The JavaScript Object Notation (JSON) Data Interchange
//!   Format. RFC 8259.
//! - Fisher, M. (2020). Reproducibility in Data Science. O'Reilly Media.

use super::error::PodProductsError;
use crate::pod::run::manifest::RunManifest;
use std::io::Write;

/// Writer that serialises a [`RunManifest`] to deterministic JSON.
///
/// # Examples
///
/// ```
/// use siderust::pod::run::manifest::RunManifest;
/// use siderust::pod::product::manifest::ManifestWriter;
///
/// let manifest = RunManifest {
///     run_id: "test-run".into(),
///     tool_version: "siderust 0.8.0".into(),
///     config_sha256: "0".repeat(64),
///     inputs: vec![],
///     outputs: vec![],
///     started_at: "2026-01-01T00:00:00Z".into(),
///     finished_at: "2026-01-01T00:01:00Z".into(),
/// };
/// let mut buf = Vec::new();
/// ManifestWriter::write(&mut buf, &manifest).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.contains("run_id"));
/// assert!(text.contains("test-run"));
/// ```
pub struct ManifestWriter;

impl ManifestWriter {
    /// Serialise `manifest` to pretty-printed JSON and write it to `w`.
    ///
    /// Inputs and outputs are sorted into canonical order before serialisation
    /// so that the output is deterministic regardless of insertion order.
    ///
    /// The output is terminated with a single newline.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Json`] if `serde_json` serialisation fails
    /// (in practice this can only happen for non-finite floats in the struct,
    /// which the `RunManifest` type does not contain).
    /// Returns [`PodProductsError::Io`] if writing to `w` fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::run::manifest::RunManifest;
    /// use siderust::pod::product::manifest::ManifestWriter;
    ///
    /// let m = RunManifest {
    ///     run_id: "00000000-0000-0000-0000-000000000000".into(),
    ///     tool_version: "siderust 0.8.0".into(),
    ///     config_sha256: "a".repeat(64),
    ///     inputs: vec![],
    ///     outputs: vec![],
    ///     started_at: "2026-01-01T00:00:00Z".into(),
    ///     finished_at: "2026-01-01T00:00:01Z".into(),
    /// };
    /// let mut buf = Vec::new();
    /// ManifestWriter::write(&mut buf, &m).unwrap();
    /// assert!(String::from_utf8(buf).unwrap().ends_with('\n'));
    /// ```
    pub fn write<W: Write>(w: &mut W, manifest: &RunManifest) -> Result<(), PodProductsError> {
        let json = manifest.to_json_pretty()?;
        w.write_all(json.as_bytes())?;
        w.write_all(b"\n")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_manifest() -> RunManifest {
        RunManifest {
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            tool_version: "siderust 0.8.0".into(),
            config_sha256: "0".repeat(64),
            inputs: vec![],
            outputs: vec![],
            started_at: "2026-01-01T00:00:00Z".into(),
            finished_at: "2026-01-01T00:00:01Z".into(),
        }
    }

    #[test]
    fn write_produces_valid_json() {
        let mut buf = Vec::new();
        ManifestWriter::write(&mut buf, &sample_manifest()).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let v: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(
            v["run_id"].as_str().unwrap(),
            "00000000-0000-0000-0000-000000000000"
        );
    }

    #[test]
    fn write_is_deterministic() {
        let mut a = Vec::new();
        let mut b = Vec::new();
        ManifestWriter::write(&mut a, &sample_manifest()).unwrap();
        ManifestWriter::write(&mut b, &sample_manifest()).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn write_ends_with_newline() {
        let mut buf = Vec::new();
        ManifestWriter::write(&mut buf, &sample_manifest()).unwrap();
        assert!(buf.ends_with(b"\n"));
    }

    #[test]
    fn write_contains_all_required_fields() {
        let mut buf = Vec::new();
        ManifestWriter::write(&mut buf, &sample_manifest()).unwrap();
        let text = String::from_utf8(buf).unwrap();
        for field in [
            "run_id",
            "tool_version",
            "config_sha256",
            "inputs",
            "outputs",
            "started_at",
            "finished_at",
        ] {
            assert!(text.contains(field), "missing field: {field}");
        }
    }
}
