// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::coordinates::spherical::direction::ICRS;
use crate::qtty::{Degree, JulianYears, Radian, Radians};

use super::{CatalogRecord, LargeStarCatalog};

/// Candidate selection for each row of the left catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossMatchMode {
    /// Return every pair within the search radius.
    AllCandidates,
    /// Return the closest candidate per left row. Exact ties favor the lowest
    /// right row index. Several left rows may select the same right row.
    Nearest,
}

/// Policy for rows missing one or both proper-motion components when the
/// comparison epoch differs from the row's reference epoch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingProperMotion {
    /// Exclude the row from candidate evaluation.
    Skip,
    /// Use the recorded position, ignoring any available partial motion.
    AssumeStationary,
    /// Return an error identifying the incomplete row.
    Error,
}

/// Configuration for [`LargeStarCatalog::cross_match`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrossMatchOptions {
    /// Common reference epoch in Julian years (e.g. J2016.0).
    ///
    /// Supply all epochs on a consistent time scale. Gaia DR3 uses TCB.
    pub epoch: JulianYears,
    /// Inclusive angular search radius in [0, π] radians.
    pub radius: Radians,
    /// Candidate selection per left row.
    pub mode: CrossMatchMode,
    /// Treatment of incomplete proper motion away from the reference epoch.
    pub missing_proper_motion: MissingProperMotion,
}

/// Identity and motion information for one endpoint of a candidate pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CrossMatchSource {
    /// Zero-based row index in the corresponding input catalog.
    pub index: usize,
    /// Source identifier as stored in that catalog. IDs can repeat.
    pub source_id: u64,
    /// Whether incomplete proper motion required a stationary assumption.
    pub assumed_stationary: bool,
}

/// One geometric candidate pair evaluated at [`CrossMatchOptions::epoch`].
///
/// The two counts describe all pairs inside the search radius, including
/// candidates omitted by [`CrossMatchMode::Nearest`]. A count greater than one
/// indicates geometric ambiguity on that side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrossMatchCandidate {
    /// Row from the catalog on which the search was called.
    pub left: CrossMatchSource,
    /// Row from the catalog supplied as the other input.
    pub right: CrossMatchSource,
    /// Great-circle separation of the two directions at the comparison epoch.
    pub separation: Radians,
    /// Number of right rows within the radius of this left row.
    pub left_candidate_count: usize,
    /// Number of left rows within the radius of this right row.
    pub right_candidate_count: usize,
}

/// Input catalog containing the row involved in a cross-match error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossMatchSide {
    /// Catalog on which the search was called.
    Left,
    /// Catalog supplied as the other input.
    Right,
}

/// Invalid input or unavailable astrometry during cross-matching.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CrossMatchError {
    /// The comparison epoch must be finite.
    #[error("comparison epoch must be finite")]
    InvalidEpoch,
    /// The angular radius must be finite and lie in the interval [0, π].
    #[error("search radius must be finite and in [0, π] radians")]
    InvalidRadius,
    /// A value used in astrometric evaluation is out of range or non-finite.
    #[error("invalid {field} in {side:?} catalog at row {index}")]
    InvalidRecord {
        /// Input catalog containing the invalid row.
        side: CrossMatchSide,
        /// Zero-based row index in that catalog.
        index: usize,
        /// Catalog field or propagation step that failed validation.
        field: &'static str,
    },
    /// An off-epoch row has incomplete motion and the policy requires an error.
    #[error("missing proper motion in {side:?} catalog at row {index}")]
    MissingProperMotion {
        /// Input catalog containing the incomplete row.
        side: CrossMatchSide,
        /// Zero-based row index in that catalog.
        index: usize,
    },
}

impl LargeStarCatalog {
    /// Find geometric candidate pairs in two ICRS stellar catalogs.
    ///
    /// Evaluates each row at most once at the common epoch via
    /// [`CatalogRecord::propagate_to`]. Rows already at the requested epoch
    /// need no motion estimate. See [`MissingProperMotion`] for other rows.
    ///
    /// Results are ordered by left row index, then increasing separation,
    /// then right row index. Unmatched and skipped rows emit no candidates.
    /// Candidate counts describe geometric ambiguity. Physical association
    /// requires further evidence, and this search permits many-to-many pairs.
    ///
    /// The right catalog is sorted by declination to restrict the angular
    /// tests to a declination window. For L left rows, R right rows and P
    /// pairs in those windows, search time is `O(R log(R+1) + L log(R+2) + P)`.
    /// All-candidates mode also sorts each left row's candidates. Storage is
    /// `O(R)` plus the returned candidates. In the worst case, every left row's
    /// window contains all right rows.
    ///
    /// # Errors
    ///
    /// Rejects non-finite epochs, invalid radii, invalid positions or
    /// non-finite motion components (including in skipped rows), and
    /// non-finite propagation results. Accepts right ascension in [0, 2π]
    /// and declination in [−π/2, π/2]. With [`MissingProperMotion::Error`],
    /// an incomplete off-epoch row also causes an error. Input validation
    /// applies even when the opposite catalog is empty.
    pub fn cross_match(
        &self,
        other: &Self,
        options: CrossMatchOptions,
    ) -> Result<Vec<CrossMatchCandidate>, CrossMatchError> {
        if !options.epoch.value().is_finite() {
            return Err(CrossMatchError::InvalidEpoch);
        }
        if !(0.0..=PI).contains(&options.radius.value()) {
            return Err(CrossMatchError::InvalidRadius);
        }

        let mut right_rows = Vec::with_capacity(other.len());
        for (index, record) in other.records().iter().enumerate() {
            if let Some(row) = evaluate(record, index, CrossMatchSide::Right, options)? {
                right_rows.push(row);
            }
        }
        right_rows.sort_unstable_by(|a, b| a.dec_deg().total_cmp(&b.dec_deg()));

        // A great-circle distance is at least the absolute declination
        // difference. Pad only the window for rounding in the unit conversions
        // and separation calculation, then apply the exact radius test.
        let window = options.radius.to::<Degree>().value() + 16.0 * f64::EPSILON * 180.0;
        let mut right_counts = vec![0; other.len()];
        let mut result = Vec::new();
        for (index, record) in self.records().iter().enumerate() {
            let Some(left) = evaluate(record, index, CrossMatchSide::Left, options)? else {
                continue;
            };
            let start = result.len();
            let lower = left.dec_deg() - window;
            let upper = left.dec_deg() + window;
            let first = right_rows.partition_point(|row| row.dec_deg() < lower);
            let mut left_count = 0;
            let mut nearest = None;
            for right in right_rows[first..]
                .iter()
                .take_while(|row| row.dec_deg() <= upper)
            {
                let separation = left
                    .direction
                    .angular_separation(&right.direction)
                    .to::<Radian>();
                if separation > options.radius {
                    continue;
                }
                left_count += 1;
                right_counts[right.source.index] += 1;
                let candidate = CrossMatchCandidate {
                    left: left.source,
                    right: right.source,
                    separation,
                    left_candidate_count: 0,
                    right_candidate_count: 0,
                };
                match options.mode {
                    CrossMatchMode::AllCandidates => result.push(candidate),
                    CrossMatchMode::Nearest => {
                        if nearest
                            .as_ref()
                            .is_none_or(|best| by_separation(&candidate, best).is_lt())
                        {
                            nearest = Some(candidate);
                        }
                    }
                }
            }
            if let Some(best) = nearest {
                result.push(best);
            }
            let candidates = &mut result[start..];
            candidates.sort_unstable_by(by_separation);
            for candidate in candidates {
                candidate.left_candidate_count = left_count;
            }
        }
        for candidate in &mut result {
            candidate.right_candidate_count = right_counts[candidate.right.index];
        }
        Ok(result)
    }
}

fn by_separation(a: &CrossMatchCandidate, b: &CrossMatchCandidate) -> core::cmp::Ordering {
    a.separation
        .value()
        .total_cmp(&b.separation.value())
        .then_with(|| a.right.index.cmp(&b.right.index))
}

struct EvaluatedRow {
    source: CrossMatchSource,
    direction: ICRS,
}

impl EvaluatedRow {
    fn dec_deg(&self) -> f64 {
        self.direction.dec().value()
    }
}

fn evaluate(
    record: &CatalogRecord,
    index: usize,
    side: CrossMatchSide,
    options: CrossMatchOptions,
) -> Result<Option<EvaluatedRow>, CrossMatchError> {
    let invalid = |field| CrossMatchError::InvalidRecord { side, index, field };
    if !record.epoch_jyr.is_finite() {
        return Err(invalid("epoch_jyr"));
    }
    if !(0.0..=TAU).contains(&record.ra.value()) {
        return Err(invalid("ra"));
    }
    if !(-FRAC_PI_2..=FRAC_PI_2).contains(&record.dec.value()) {
        return Err(invalid("dec"));
    }
    for (field, value) in [
        ("pm_ra_cosdec", record.pm_ra_cosdec),
        ("pm_dec", record.pm_dec),
    ] {
        if value.is_some_and(|motion| !motion.is_finite()) {
            return Err(invalid(field));
        }
    }

    let off_epoch = record.epoch_jyr != options.epoch.value();
    let assumed_stationary =
        off_epoch && (record.pm_ra_cosdec.is_none() || record.pm_dec.is_none());
    if assumed_stationary {
        match options.missing_proper_motion {
            MissingProperMotion::Skip => return Ok(None),
            MissingProperMotion::Error => {
                return Err(CrossMatchError::MissingProperMotion { side, index });
            }
            MissingProperMotion::AssumeStationary => {}
        }
    }
    let (ra, dec) = if off_epoch && !assumed_stationary {
        record.propagate_to(options.epoch.value())
    } else {
        (record.ra, record.dec)
    };
    if !ra.value().is_finite() || !dec.value().is_finite() {
        return Err(invalid("epoch propagation"));
    }
    // All RAs represent the same direction at an exact pole. Keep the
    // recorded RA in propagation above, where it defines the tangent basis.
    let ra = if dec.value().abs() == FRAC_PI_2 || ra.value() == TAU {
        Radians::new(0.0)
    } else {
        ra
    };
    Ok(Some(EvaluatedRow {
        source: CrossMatchSource {
            index,
            source_id: record.source_id,
            assumed_stationary,
        },
        // The angles are already canonical. Avoid wrapping declination via
        // addition/subtraction of 180°, which rounds away tiny differences.
        direction: ICRS::new_unchecked(dec.to::<Degree>(), ra.to::<Degree>().normalize()),
    }))
}
