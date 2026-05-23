// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! [`SpiceEphemerisProvider`] — typed adapter implementing
//! [`siderust::pod::providers::EphemerisProvider`].
//!
//! Center selection is **explicit at construction time**: the caller
//! must pick a typed `affn` reference center (e.g.
//! `siderust::coordinates::centers::Barycentric`) and the matching NAIF
//! body id. Every state returned by the provider carries that center on
//! its typed [`affn::cartesian::Position`], so downstream POD code
//! cannot accidentally mix barycentric and heliocentric vectors.

use std::marker::PhantomData;
use std::sync::Arc;

use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::ICRS;
use crate::pod::providers::EphemerisProvider;
use affn::cartesian::Position;
use qtty::unit::Kilometer;
use qtty::Quantity;
use tempoch::{EncodedTime, J2000s, Time, TDB};

use crate::formats::spice::{SpiceError, SpkKernel};

/// One body state returned by [`SpiceEphemerisProvider`].
///
/// The position is typed by the construction-time reference center `C`,
/// so adding a barycentric position to a heliocentric position fails to
/// compile via `affn`'s affine algebra rules. Velocity is returned as a
/// raw `[f64; 3]` (km/s) for now; once the workspace `affn::Velocity`
/// or `qtty::Velocity` types stabilize for the J2000-reference-frame
/// case, this can be tightened without breaking the public surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpiceState<C: ReferenceCenter<Params = ()>> {
    /// Position in the kernel's frame (typically ICRS / J2000),
    /// expressed in kilometers from the construction-time center.
    pub position: Position<C, ICRS, Kilometer>,
    /// Velocity in km/s in the same frame.
    pub velocity_km_s: [f64; 3],
}

/// Typed `EphemerisProvider` backed by an in-memory [`SpkKernel`].
///
/// Construct with the NAIF id of the desired reference center and the
/// matching `affn` typed center marker. Cloning is cheap (`Arc`-shared
/// kernel).
///
/// # Examples
///
/// ```rust,no_run
/// use siderust::formats::spice::{SpkKernel, well_known};
/// use siderust::pod::spice::SpiceEphemerisProvider;
/// use siderust::coordinates::centers::Barycentric;
/// use siderust::pod::providers::EphemerisProvider;
///
/// let kernel = SpkKernel::open("de440.bsp")?;
/// let provider: SpiceEphemerisProvider<Barycentric> =
///     SpiceEphemerisProvider::new(kernel, well_known::SSB);
/// // Earth state at J2000:
/// let state = provider.state(well_known::EARTH, 0.0)?;
/// let _x_km = state.position.x();
/// # Ok::<_, siderust::formats::spice::SpiceError>(())
/// ```
pub struct SpiceEphemerisProvider<C: ReferenceCenter<Params = ()>> {
    kernel: Arc<SpkKernel>,
    center_naif_id: i32,
    _center: PhantomData<fn() -> C>,
}

impl<C: ReferenceCenter<Params = ()>> Clone for SpiceEphemerisProvider<C> {
    fn clone(&self) -> Self {
        Self {
            kernel: Arc::clone(&self.kernel),
            center_naif_id: self.center_naif_id,
            _center: PhantomData,
        }
    }
}

impl<C: ReferenceCenter<Params = ()>> std::fmt::Debug for SpiceEphemerisProvider<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpiceEphemerisProvider")
            .field("center_name", &C::center_name())
            .field("center_naif_id", &self.center_naif_id)
            .field("kernel", &self.kernel)
            .finish()
    }
}

impl<C: ReferenceCenter<Params = ()>> SpiceEphemerisProvider<C> {
    /// Construct a provider that returns positions relative to
    /// `center_naif_id`, typed at `C`.
    ///
    /// The caller is responsible for matching the typed marker `C`
    /// (e.g. `Barycentric`) with the NAIF id (e.g. `0`); this crate
    /// does not enforce that mapping at compile time because POD
    /// adapters may want to use custom centers.
    pub fn new(kernel: SpkKernel, center_naif_id: i32) -> Self {
        Self {
            kernel: Arc::new(kernel),
            center_naif_id,
            _center: PhantomData,
        }
    }

    /// Construct from an `Arc`-shared kernel (useful when several
    /// providers share one kernel with different centers).
    pub fn from_arc(kernel: Arc<SpkKernel>, center_naif_id: i32) -> Self {
        Self {
            kernel,
            center_naif_id,
            _center: PhantomData,
        }
    }

    /// NAIF id of the construction-time reference center.
    pub fn center_naif_id(&self) -> i32 {
        self.center_naif_id
    }

    /// Borrow the underlying kernel.
    pub fn kernel(&self) -> &SpkKernel {
        &self.kernel
    }

    /// Convenience: query state at a typed `tempoch` TDB epoch.
    ///
    /// Equivalent to [`EphemerisProvider::state`] but accepts a typed
    /// [`Time<TDB>`] instead of a raw `f64` so callers using the
    /// `tempoch` time-stack stay typed end-to-end.
    pub fn state_at(
        &self,
        body_naif_id: i32,
        epoch: Time<TDB>,
    ) -> Result<SpiceState<C>, SpiceError> {
        let encoded: EncodedTime<TDB, J2000s> = epoch.to::<J2000s>();
        let secs: Quantity<qtty::unit::Second> = encoded.raw();
        EphemerisProvider::state(self, body_naif_id, secs.value())
    }
}

impl<C: ReferenceCenter<Params = ()>> EphemerisProvider for SpiceEphemerisProvider<C> {
    type State = SpiceState<C>;
    type Error = SpiceError;

    fn state(&self, body_naif_id: i32, epoch_seconds_tdb: f64) -> Result<Self::State, Self::Error> {
        let s = self
            .kernel
            .state(body_naif_id, self.center_naif_id, epoch_seconds_tdb)?;
        let position = Position::<C, ICRS, Kilometer>::new(
            Quantity::<Kilometer>::new(s[0]),
            Quantity::<Kilometer>::new(s[1]),
            Quantity::<Kilometer>::new(s[2]),
        );
        Ok(SpiceState {
            position,
            velocity_km_s: [s[3], s[4], s[5]],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::{Barycentric, Geocentric};

    /// Same synthetic kernel used in `kernel::tests`.
    fn synthetic_kernel() -> SpkKernel {
        let mut buf = vec![0u8; 3 * 1024];
        buf[0..8].copy_from_slice(b"DAF/SPK ");
        buf[8..12].copy_from_slice(&2i32.to_le_bytes());
        buf[12..16].copy_from_slice(&6i32.to_le_bytes());
        buf[76..80].copy_from_slice(&2i32.to_le_bytes());
        let rec = &mut buf[1024..2048];
        rec[16..24].copy_from_slice(&2.0_f64.to_le_bytes());

        let seg1_start_word: i32 = 257;
        let seg1_end_word: i32 = seg1_start_word + 8;
        let seg2_start_word: i32 = seg1_end_word + 1;
        let seg2_end_word: i32 = seg2_start_word + 8;

        let s1 = &mut rec[24..64];
        s1[0..8].copy_from_slice(&0.0_f64.to_le_bytes());
        s1[8..16].copy_from_slice(&1.0_f64.to_le_bytes());
        s1[16..20].copy_from_slice(&399i32.to_le_bytes());
        s1[20..24].copy_from_slice(&0i32.to_le_bytes());
        s1[24..28].copy_from_slice(&1i32.to_le_bytes());
        s1[28..32].copy_from_slice(&2i32.to_le_bytes());
        s1[32..36].copy_from_slice(&seg1_start_word.to_le_bytes());
        s1[36..40].copy_from_slice(&seg1_end_word.to_le_bytes());

        let s2 = &mut rec[64..104];
        s2[0..8].copy_from_slice(&0.0_f64.to_le_bytes());
        s2[8..16].copy_from_slice(&1.0_f64.to_le_bytes());
        s2[16..20].copy_from_slice(&301i32.to_le_bytes());
        s2[20..24].copy_from_slice(&399i32.to_le_bytes());
        s2[24..28].copy_from_slice(&1i32.to_le_bytes());
        s2[28..32].copy_from_slice(&2i32.to_le_bytes());
        s2[32..36].copy_from_slice(&seg2_start_word.to_le_bytes());
        s2[36..40].copy_from_slice(&seg2_end_word.to_le_bytes());

        let wf = |word: i32, v: f64, b: &mut [u8]| {
            let off = (word as usize - 1) * 8;
            b[off..off + 8].copy_from_slice(&v.to_le_bytes());
        };
        wf(seg1_start_word, 0.5, &mut buf);
        wf(seg1_start_word + 1, 0.5, &mut buf);
        wf(seg1_start_word + 2, 1.0e8, &mut buf);
        wf(seg1_start_word + 3, 0.0, &mut buf);
        wf(seg1_start_word + 4, 0.0, &mut buf);
        wf(seg1_end_word - 3, 0.0, &mut buf);
        wf(seg1_end_word - 2, 1.0, &mut buf);
        wf(seg1_end_word - 1, 5.0, &mut buf);
        wf(seg1_end_word, 1.0, &mut buf);
        wf(seg2_start_word, 0.5, &mut buf);
        wf(seg2_start_word + 1, 0.5, &mut buf);
        wf(seg2_start_word + 2, 0.0, &mut buf);
        wf(seg2_start_word + 3, 3.84e5, &mut buf);
        wf(seg2_start_word + 4, 0.0, &mut buf);
        wf(seg2_end_word - 3, 0.0, &mut buf);
        wf(seg2_end_word - 2, 1.0, &mut buf);
        wf(seg2_end_word - 1, 5.0, &mut buf);
        wf(seg2_end_word, 1.0, &mut buf);

        SpkKernel::from_bytes(buf).unwrap()
    }

    #[test]
    fn provider_returns_typed_barycentric_position_for_earth() {
        let kernel = synthetic_kernel();
        let provider: SpiceEphemerisProvider<Barycentric> = SpiceEphemerisProvider::new(kernel, 0);
        let state = provider.state(399, 0.5).unwrap();
        assert!((state.position.x().value() - 1.0e8).abs() < 1e-6);
        assert!(state.position.y().value().abs() < 1e-6);
    }

    #[test]
    fn provider_returns_typed_geocentric_position_for_moon() {
        let kernel = synthetic_kernel();
        let provider: SpiceEphemerisProvider<Geocentric> = SpiceEphemerisProvider::new(kernel, 399);
        let state = provider.state(301, 0.5).unwrap();
        assert!((state.position.y().value() - 3.84e5).abs() < 1e-6);
        assert!(state.position.x().value().abs() < 1e-6);
    }

    #[test]
    fn cloned_provider_shares_kernel() {
        let kernel = synthetic_kernel();
        let p1: SpiceEphemerisProvider<Barycentric> = SpiceEphemerisProvider::new(kernel, 0);
        let p2 = p1.clone();
        let s1 = p1.state(399, 0.5).unwrap();
        let s2 = p2.state(399, 0.5).unwrap();
        assert_eq!(s1.position.x().value(), s2.position.x().value());
        assert_eq!(p1.center_naif_id(), 0);
        assert_eq!(p2.center_naif_id(), 0);
    }

    #[test]
    fn debug_repr_lists_center_name() {
        let kernel = synthetic_kernel();
        let provider: SpiceEphemerisProvider<Barycentric> = SpiceEphemerisProvider::new(kernel, 0);
        let dbg = format!("{provider:?}");
        assert!(dbg.contains("Barycentric"));
        assert!(dbg.contains("center_naif_id: 0"));
    }
}
