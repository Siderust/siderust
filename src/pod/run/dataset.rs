//! Dataset reference with content hashing for run manifests.

use std::path::{Path, PathBuf};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use sha2::{Digest, Sha256};

/// Reference to an input or output dataset used by a POD run.
///
/// A `DatasetRef` records the canonical path, an opaque kind tag (e.g.
/// `"sp3"`, `"rinex_obs"`, `"qc_json"`), the byte length, and the SHA-256
/// content hash. It is the building block of [`crate::pod::run::manifest::RunManifest`].
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DatasetRef {
    /// Path to the dataset file.
    pub path: PathBuf,
    /// Opaque content type tag (e.g., `"sp3"`, `"rinex_obs"`).
    pub kind: String,
    /// File size in bytes.
    pub bytes: u64,
    /// Lower-case hex SHA-256 digest of the file contents.
    pub sha256: String,
}

impl DatasetRef {
    /// Hash a file on disk and build a [`DatasetRef`].
    pub fn from_file(path: impl AsRef<Path>, kind: impl Into<String>) -> std::io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let bytes_vec = std::fs::read(&path)?;
        let bytes = bytes_vec.len() as u64;
        let mut hasher = Sha256::new();
        hasher.update(&bytes_vec);
        let sha256 = hex::encode(hasher.finalize());
        Ok(Self {
            path,
            kind: kind.into(),
            bytes,
            sha256,
        })
    }

    /// Hash a byte buffer directly (useful for in-memory artifacts).
    pub fn from_bytes(path: impl Into<PathBuf>, kind: impl Into<String>, data: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(data);
        Self {
            path: path.into(),
            kind: kind.into(),
            bytes: data.len() as u64,
            sha256: hex::encode(hasher.finalize()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_bytes_hashes_deterministically() {
        let a = DatasetRef::from_bytes("/tmp/a", "test", b"hello world");
        let b = DatasetRef::from_bytes("/tmp/a", "test", b"hello world");
        assert_eq!(a, b);
        // Known SHA-256 of "hello world".
        assert_eq!(
            a.sha256,
            "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        );
        assert_eq!(a.bytes, 11);
    }
}
