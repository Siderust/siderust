// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use alloc::vec;
use alloc::vec::Vec;

use crate::coordinates::spherical::direction::ICRS;
use crate::qtty::{Degree, Degrees, JulianYears, Radian, Radians};

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
///
/// The comparison epoch and search radius are required. By default the matcher
/// returns [all candidates](CrossMatchMode::AllCandidates) and
/// [returns an error](MissingProperMotion::Error) if an off-epoch row has
/// incomplete proper motion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrossMatchOptions {
    epoch: JulianYears,
    radius: Radians,
    mode: CrossMatchMode,
    missing_proper_motion: MissingProperMotion,
}

impl CrossMatchOptions {
    /// Create options for a common epoch and an inclusive search radius in [0, π].
    ///
    /// Supply the numeric Julian epoch (e.g. `JulianYears::new(2016.0)` for
    /// J2016.0), using the same convention as [`CatalogRecord::epoch_jyr`].
    /// The quantity carries the unit, and the caller supplies epochs for both
    /// catalogs and the comparison on a common time scale. Gaia DR3 uses TCB.
    ///
    /// [`LargeStarCatalog::cross_match`] validates the epoch and radius.
    pub const fn new(epoch: JulianYears, radius: Radians) -> Self {
        Self {
            epoch,
            radius,
            mode: CrossMatchMode::AllCandidates,
            missing_proper_motion: MissingProperMotion::Error,
        }
    }

    /// Set the candidate selection mode for each left row.
    #[must_use]
    pub const fn with_mode(mut self, mode: CrossMatchMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set the policy for incomplete motion on off-epoch rows.
    #[must_use]
    pub const fn with_missing_proper_motion(mut self, policy: MissingProperMotion) -> Self {
        self.missing_proper_motion = policy;
        self
    }

    /// Return the comparison epoch as a numeric Julian epoch.
    pub const fn epoch(&self) -> JulianYears {
        self.epoch
    }

    /// Return the inclusive angular search radius.
    pub const fn radius(&self) -> Radians {
        self.radius
    }

    /// Return the candidate selection mode.
    pub const fn mode(&self) -> CrossMatchMode {
        self.mode
    }

    /// Return the policy for incomplete motion on off-epoch rows.
    pub const fn missing_proper_motion(&self) -> MissingProperMotion {
        self.missing_proper_motion
    }
}

/// Identity and motion information for one endpoint of a candidate pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
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
#[non_exhaustive]
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
#[non_exhaustive]
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
        /// Diagnostic label for the catalog field or propagation step.
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
        if !options.epoch.is_finite() {
            return Err(CrossMatchError::InvalidEpoch);
        }
        if !(Radians::zero()..=Radians::HALF_TURN).contains(&options.radius) {
            return Err(CrossMatchError::InvalidRadius);
        }

        let mut right_rows = Vec::with_capacity(other.len());
        for (index, record) in other.records().iter().enumerate() {
            if let Some(row) = evaluate(record, index, CrossMatchSide::Right, options)? {
                right_rows.push(row);
            }
        }
        right_rows.sort_unstable_by(|a, b| {
            a.direction
                .dec()
                .value()
                .total_cmp(&b.direction.dec().value())
        });

        // A great-circle distance is at least the absolute declination
        // difference. Pad only the window for rounding in the unit conversions
        // and separation calculation, then apply the exact radius test.
        let window = options.radius.to::<Degree>() + Degrees::new(16.0 * f64::EPSILON * 180.0);
        let mut right_counts = vec![0; other.len()];
        let mut result = Vec::new();
        for (index, record) in self.records().iter().enumerate() {
            let Some(left) = evaluate(record, index, CrossMatchSide::Left, options)? else {
                continue;
            };
            let start = result.len();
            let lower = left.direction.dec() - window;
            let upper = left.direction.dec() + window;
            let first = right_rows.partition_point(|row| row.direction.dec() < lower);
            let mut left_count = 0;
            let mut nearest = None;
            for right in right_rows[first..]
                .iter()
                .take_while(|row| row.direction.dec() <= upper)
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
    if !(Radians::zero()..=Radians::FULL_TURN).contains(&record.ra) {
        return Err(invalid("ra"));
    }
    if !(-Radians::QUARTER_TURN..=Radians::QUARTER_TURN).contains(&record.dec) {
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

    let epoch_jyr = options.epoch.value();
    let off_epoch = record.epoch_jyr != epoch_jyr;
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
        record.propagate_to(epoch_jyr)
    } else {
        (record.ra, record.dec)
    };
    if !ra.is_finite() || !dec.is_finite() {
        return Err(invalid("epoch propagation"));
    }
    // All RAs represent the same direction at an exact pole. Keep the
    // recorded RA in propagation above, where it defines the tangent basis.
    let ra = if dec.abs() == Radians::QUARTER_TURN || ra == Radians::FULL_TURN {
        Radians::zero()
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
