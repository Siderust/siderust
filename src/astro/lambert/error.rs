//! Unified error type for the Lambert solver.
//!
//! All public APIs in this crate fail through a single [`LambertError`]
//! variant. The variants are stable and carry enough context to let
//! callers decide whether to retry with different inputs (e.g. relax
//! `revolutions`) or surface a domain error to the user.

/// Errors returned by the Lambert solver.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum LambertError {
    /// Gravitational parameter must be strictly positive (km³/s²).
    #[error("non-positive gravitational parameter ({0})")]
    NonPositiveMu(f64),
    /// Time of flight must be a finite, strictly positive duration.
    #[error("non-positive time of flight ({0} s)")]
    NonPositiveTof(f64),
    /// Initial or final position vector has effectively zero magnitude.
    #[error("zero-magnitude position vector")]
    ZeroPosition,
    /// `r1`, `r2` and the origin are collinear, so the chord/transfer
    /// plane is undefined.
    #[error("collinear positions: chord cannot be resolved unambiguously")]
    Collinear,
    /// The requested revolution count exceeds the maximum permitted by
    /// the supplied time of flight (Izzo Eq. 21).
    ///
    /// `requested` is what the caller asked for; `max` is the largest
    /// `N` for which a solution exists at the supplied `tof`.
    #[error("requested {requested} revolutions exceeds N_max = {max} for the given TOF")]
    RevolutionsExceedNMax {
        /// Number of revolutions requested by the caller.
        requested: u32,
        /// Largest physically admissible revolution count.
        max: u32,
    },
    /// Householder iteration failed to converge within the iteration cap.
    ///
    /// Carries the residual at the last iterate (non-dimensional time
    /// units).
    #[error("Householder iteration failed to converge (residual = {residual:e})")]
    DidNotConverge {
        /// Final residual `|T(x) − T*|` in non-dimensional units.
        residual: f64,
    },
}
