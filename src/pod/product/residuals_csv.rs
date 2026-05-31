//! # Residual CSV writer
//!
//! ## Scientific scope
//!
//! Residual time series are a standard diagnostic product in POD,
//! allowing analysts to inspect measurement fit quality by epoch,
//! sensor and observable type. This module serializes residuals
//! computed by the estimation/QC layers into a stable flat table for
//! downstream plotting or audit workflows.
//!
//! ## Technical scope
//!
//! A single streaming writer is exposed:
//!
//! - [`ResidualRecord`] is the column model
//!   `(epoch_jd_tt, obs_type, satellite, residual_m, sigma_m, rejected)`.
//! - [`ResidualCsvWriter`] flushes the sink after every record so a
//!   process crash cannot silently drop already-written residuals.
//!
//! Formatting is intentionally simple and stable; statistical
//! aggregation remains the responsibility of [`crate::pod::qc`].
//!
//! ## References
//!
//! - Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical
//!   Orbit Determination*. Elsevier Academic Press.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and
//!   Applications* (4th ed.). Microcosm Press.

use super::error::PodProductsError;
use serde::{Deserialize, Serialize};
use std::io::Write;

/// One residual record in the standard POD column format.
///
/// The column order in the CSV output is:
/// `epoch_jd_tt, obs_type, satellite, residual_m, sigma_m, rejected`.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::residuals_csv::ResidualRecord;
///
/// let r = ResidualRecord {
///     epoch_jd_tt: 2_451_545.0,
///     obs_type: "C1C".into(),
///     satellite: "G01".into(),
///     residual_m: -0.15,
///     sigma_m: 0.30,
///     rejected: false,
/// };
/// assert_eq!(r.satellite, "G01");
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidualRecord {
    /// Epoch expressed as a Julian Date in TT time scale.
    pub epoch_jd_tt: f64,
    /// RINEX observation type code, e.g. `"C1C"`, `"L1C"`, `"P2"`.
    pub obs_type: String,
    /// Satellite identifier, e.g. `"G01"`, `"R05"`, `"E11"`.
    pub satellite: String,
    /// Post-fit residual in metres (positive = observed − computed).
    pub residual_m: f64,
    /// Measurement sigma used in the estimator (metres).
    pub sigma_m: f64,
    /// Whether this measurement was rejected by the outlier filter.
    pub rejected: bool,
}

/// Streaming CSV writer for POD residuals.
///
/// Writes one [`ResidualRecord`] at a time and flushes the underlying
/// sink after every record so that a crash preserves all records
/// written so far. The header is written on construction.
///
/// # Examples
///
/// ```
/// use siderust::pod::product::residuals_csv::{ResidualCsvWriter, ResidualRecord};
///
/// let mut buf = Vec::<u8>::new();
/// let mut w = ResidualCsvWriter::new(&mut buf).unwrap();
/// w.write_record(&ResidualRecord {
///     epoch_jd_tt: 2_451_545.0,
///     obs_type: "C1C".into(),
///     satellite: "G01".into(),
///     residual_m: -0.15,
///     sigma_m: 0.30,
///     rejected: false,
/// }).unwrap();
/// w.flush().unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.lines().next().unwrap().starts_with("epoch_jd_tt,"));
/// assert_eq!(text.lines().count(), 2);
/// ```
pub struct ResidualCsvWriter<W: Write> {
    writer: W,
}

impl<W: Write> ResidualCsvWriter<W> {
    /// Create a new writer, immediately writing the CSV header line.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Io`] if writing or flushing the
    /// header to `writer` fails.
    pub fn new(mut writer: W) -> Result<Self, PodProductsError> {
        writeln!(
            writer,
            "epoch_jd_tt,obs_type,satellite,residual_m,sigma_m,rejected"
        )?;
        writer.flush()?;
        Ok(Self { writer })
    }

    /// Write one residual record and flush the sink.
    ///
    /// The flush after each record guarantees that a process crash
    /// cannot silently drop already-written residuals.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Io`] if the write or flush fails.
    pub fn write_record(&mut self, r: &ResidualRecord) -> Result<(), PodProductsError> {
        writeln!(
            self.writer,
            "{:.9},{},{},{:.6},{:.6},{}",
            r.epoch_jd_tt, r.obs_type, r.satellite, r.residual_m, r.sigma_m, r.rejected,
        )?;
        self.writer.flush()?;
        Ok(())
    }

    /// Manually flush the underlying writer.
    ///
    /// Normally not needed because [`write_record`][Self::write_record]
    /// already flushes, but useful when the caller wants an explicit
    /// sync point.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Io`] if flushing fails.
    pub fn flush(&mut self) -> Result<(), PodProductsError> {
        self.writer.flush().map_err(PodProductsError::Io)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streaming_writer_header() {
        let mut buf = Vec::<u8>::new();
        ResidualCsvWriter::new(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(
            text.trim(),
            "epoch_jd_tt,obs_type,satellite,residual_m,sigma_m,rejected"
        );
    }

    #[test]
    fn streaming_writer_record_count() {
        let mut buf = Vec::<u8>::new();
        let mut w = ResidualCsvWriter::new(&mut buf).unwrap();
        for i in 0..3 {
            w.write_record(&ResidualRecord {
                epoch_jd_tt: 2_451_545.0 + i as f64,
                obs_type: "C1C".into(),
                satellite: "G01".into(),
                residual_m: 0.0,
                sigma_m: 1.0,
                rejected: false,
            })
            .unwrap();
        }
        let text = String::from_utf8(buf).unwrap();
        // 1 header + 3 data rows
        assert_eq!(text.lines().count(), 4);
    }

    #[test]
    fn streaming_writer_rejected_field() {
        let mut buf = Vec::<u8>::new();
        let mut w = ResidualCsvWriter::new(&mut buf).unwrap();
        w.write_record(&ResidualRecord {
            epoch_jd_tt: 2_451_545.0,
            obs_type: "P2".into(),
            satellite: "G03".into(),
            residual_m: 10.0,
            sigma_m: 1.0,
            rejected: true,
        })
        .unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(
            text.contains(",true"),
            "rejected=true must appear in output"
        );
    }
}
