//! # Residual Parquet writer (feature `parquet`)
//!
//! ## Scientific scope
//!
//! Parquet is a columnar binary format that supports efficient analytical
//! queries and compact storage for large residual tables. This module
//! provides the same schema as [`siderust::residuals_csv`] but in Parquet
//! format for use in big-data pipelines and analytical tools.
//!
//! The schema and column semantics are identical to the CSV surface:
//! `(epoch_jd_tt, obs_type, satellite, residual_m, sigma_m, rejected)`.
//!
//! ## Technical scope
//!
//! [`ResidualParquetWriter`] buffers records in memory and serialises them
//! as a single Parquet row group on [`finish`][ResidualParquetWriter::finish].
//! The writer requires a `Write + Seek` sink (e.g. `std::fs::File` or
//! `std::io::Cursor<Vec<u8>>`) because the Parquet format requires the
//! footer to be seeked back to write after the data blocks.
//!
//! Enable with `features = ["parquet"]` in your dependency declaration.
//!
//! ## References
//!
//! - Apache Software Foundation. (2024). Apache Parquet Format Specification.
//!   <https://parquet.apache.org/docs/file-format/>

use super::error::PodProductsError;
use super::residuals_csv::ResidualRecord;
use parquet::{
    data_type::{BoolType, ByteArray, ByteArrayType, DoubleType},
    file::{properties::WriterProperties, writer::SerializedFileWriter},
    schema::parser::parse_message_type,
};
use std::io::{Seek, Write};
use std::sync::Arc;

/// Parquet writer for POD residuals.
///
/// Records are buffered in memory and written as a single row group when
/// [`finish`][Self::finish] is called. The Parquet file will contain the
/// columns `epoch_jd_tt` (DOUBLE), `obs_type` (BYTE_ARRAY/UTF8),
/// `satellite` (BYTE_ARRAY/UTF8), `residual_m` (DOUBLE),
/// `sigma_m` (DOUBLE), and `rejected` (BOOLEAN).
///
/// # Examples
///
/// ```
/// use siderust::pod::product::residuals_parquet::ResidualParquetWriter;
/// use siderust::pod::product::residuals_csv::ResidualRecord;
/// use std::io::Cursor;
///
/// let buf = Cursor::new(Vec::<u8>::new());
/// let mut w = ResidualParquetWriter::new(buf);
/// w.push(ResidualRecord {
///     epoch_jd_tt: 2_451_545.0,
///     obs_type: "C1C".into(),
///     satellite: "G01".into(),
///     residual_m: -0.15,
///     sigma_m: 0.30,
///     rejected: false,
/// });
/// let cursor = w.finish().unwrap();
/// assert!(!cursor.into_inner().is_empty());
/// ```
pub struct ResidualParquetWriter<W: Write + Seek + Send> {
    records: Vec<ResidualRecord>,
    output: W,
}

impl<W: Write + Seek + Send> ResidualParquetWriter<W> {
    /// Create a new writer targeting `output`.
    ///
    /// No data is written until [`finish`][Self::finish] is called.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::residuals_parquet::ResidualParquetWriter;
    /// use std::io::Cursor;
    ///
    /// let _w = ResidualParquetWriter::new(Cursor::new(Vec::<u8>::new()));
    /// ```
    pub fn new(output: W) -> Self {
        Self {
            records: Vec::new(),
            output,
        }
    }

    /// Buffer one record for later writing.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::residuals_parquet::ResidualParquetWriter;
    /// use siderust::pod::product::residuals_csv::ResidualRecord;
    /// use std::io::Cursor;
    ///
    /// let mut w = ResidualParquetWriter::new(Cursor::new(Vec::<u8>::new()));
    /// w.push(ResidualRecord {
    ///     epoch_jd_tt: 2_451_545.0,
    ///     obs_type: "L1C".into(),
    ///     satellite: "E01".into(),
    ///     residual_m: 0.01,
    ///     sigma_m: 0.005,
    ///     rejected: false,
    /// });
    /// ```
    pub fn push(&mut self, record: ResidualRecord) {
        self.records.push(record);
    }

    /// Serialise all buffered records as a Parquet file and return the sink.
    ///
    /// Consuming `self` ensures that no additional records can be pushed
    /// after the file has been finalised.
    ///
    /// # Errors
    ///
    /// Returns [`PodProductsError::Parquet`] if schema parsing or any
    /// Parquet serialisation step fails, or [`PodProductsError::Io`] for
    /// underlying sink errors.
    ///
    /// # Examples
    ///
    /// ```
    /// use siderust::pod::product::residuals_parquet::ResidualParquetWriter;
    /// use siderust::pod::product::residuals_csv::ResidualRecord;
    /// use std::io::Cursor;
    ///
    /// let mut w = ResidualParquetWriter::new(Cursor::new(Vec::<u8>::new()));
    /// w.push(ResidualRecord {
    ///     epoch_jd_tt: 2_451_545.0,
    ///     obs_type: "C1C".into(),
    ///     satellite: "G01".into(),
    ///     residual_m: -0.15,
    ///     sigma_m: 0.30,
    ///     rejected: false,
    /// });
    /// let cursor = w.finish().unwrap();
    /// // The Parquet magic bytes are "PAR1".
    /// let bytes = cursor.into_inner();
    /// assert_eq!(&bytes[..4], b"PAR1");
    /// ```
    pub fn finish(self) -> Result<W, PodProductsError> {
        const SCHEMA: &str = "message residuals {
            REQUIRED DOUBLE epoch_jd_tt;
            REQUIRED BYTE_ARRAY obs_type (UTF8);
            REQUIRED BYTE_ARRAY satellite (UTF8);
            REQUIRED DOUBLE residual_m;
            REQUIRED DOUBLE sigma_m;
            REQUIRED BOOLEAN rejected;
        }";

        let schema =
            parse_message_type(SCHEMA).map_err(|e| PodProductsError::Parquet(e.to_string()))?;
        let props = WriterProperties::builder().build();
        let mut file_writer =
            SerializedFileWriter::new(self.output, Arc::new(schema), Arc::new(props))
                .map_err(|e| PodProductsError::Parquet(e.to_string()))?;

        {
            let mut rg = file_writer
                .next_row_group()
                .map_err(|e| PodProductsError::Parquet(e.to_string()))?;

            // epoch_jd_tt — DOUBLE
            {
                let vals: Vec<f64> = self.records.iter().map(|r| r.epoch_jd_tt).collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 0".into()))?;
                cw.typed::<DoubleType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            // obs_type — BYTE_ARRAY/UTF8
            {
                let vals: Vec<ByteArray> = self
                    .records
                    .iter()
                    .map(|r| ByteArray::from(r.obs_type.as_str()))
                    .collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 1".into()))?;
                cw.typed::<ByteArrayType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            // satellite — BYTE_ARRAY/UTF8
            {
                let vals: Vec<ByteArray> = self
                    .records
                    .iter()
                    .map(|r| ByteArray::from(r.satellite.as_str()))
                    .collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 2".into()))?;
                cw.typed::<ByteArrayType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            // residual_m — DOUBLE
            {
                let vals: Vec<f64> = self.records.iter().map(|r| r.residual_m).collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 3".into()))?;
                cw.typed::<DoubleType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            // sigma_m — DOUBLE
            {
                let vals: Vec<f64> = self.records.iter().map(|r| r.sigma_m).collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 4".into()))?;
                cw.typed::<DoubleType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            // rejected — BOOLEAN
            {
                let vals: Vec<bool> = self.records.iter().map(|r| r.rejected).collect();
                let mut cw = rg
                    .next_column()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?
                    .ok_or_else(|| PodProductsError::Parquet("missing column 5".into()))?;
                cw.typed::<BoolType>()
                    .write_batch(&vals, None, None)
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
                cw.close()
                    .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
            }

            rg.close()
                .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
        }

        let output = file_writer
            .into_inner()
            .map_err(|e| PodProductsError::Parquet(e.to_string()))?;
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_record(i: u32) -> ResidualRecord {
        ResidualRecord {
            epoch_jd_tt: 2_451_545.0 + f64::from(i),
            obs_type: "C1C".into(),
            satellite: "G01".into(),
            residual_m: f64::from(i) * 0.1,
            sigma_m: 0.3,
            rejected: i % 3 == 0,
        }
    }

    #[test]
    fn parquet_magic_bytes() {
        let mut w = ResidualParquetWriter::new(Cursor::new(Vec::new()));
        for i in 0..5 {
            w.push(sample_record(i));
        }
        let cursor = w.finish().unwrap();
        let bytes = cursor.into_inner();
        assert_eq!(&bytes[..4], b"PAR1", "Parquet magic bytes not found");
        assert_eq!(
            bytes.last_chunk::<4>().unwrap(),
            b"PAR1",
            "Parquet footer not found"
        );
    }

    #[test]
    fn empty_writer_produces_valid_parquet() {
        let w: ResidualParquetWriter<Cursor<Vec<u8>>> =
            ResidualParquetWriter::new(Cursor::new(Vec::new()));
        let cursor = w.finish().unwrap();
        let bytes = cursor.into_inner();
        assert_eq!(&bytes[..4], b"PAR1");
    }
}
