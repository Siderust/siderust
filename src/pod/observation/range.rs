// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Inter-satellite range helpers.

use affn::{ReferenceCenter, ReferenceFrame};
use qtty::unit::Meter;
use tempoch::{Time, TDB};

use crate::pod::observation::error::PodObservationsError;
use crate::pod::providers::EphemerisProvider;
use crate::pod::tabulated::TabulatedCartesianState;

const C_KM_S: f64 = 299_792.458;

/// Computes one-way geometric range with one light-time correction.
///
/// The receiver state is evaluated at `receive_epoch`; the transmitter state is
/// first evaluated at the receive epoch, then once more at
/// `receive_epoch - rho / c`.
pub fn one_way_light_time_range<P, C, F>(
    provider: &P,
    receiver_id: i32,
    transmitter_id: i32,
    receive_epoch: Time<TDB>,
) -> Result<qtty::Meter, PodObservationsError>
where
    P: EphemerisProvider<State = TabulatedCartesianState<C, F>>,
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let receive_seconds = receive_epoch.to::<tempoch::J2000s>().raw().value();
    let receiver = provider
        .state(receiver_id, receive_seconds)
        .map_err(|_| PodObservationsError::LightTimeNotConverged)?;
    let transmitter_at_receive = provider
        .state(transmitter_id, receive_seconds)
        .map_err(|_| PodObservationsError::LightTimeNotConverged)?;

    let initial_range_km = range_km(&receiver, &transmitter_at_receive);
    let emit_seconds = receive_seconds - initial_range_km / C_KM_S;
    let transmitter_at_emit = provider
        .state(transmitter_id, emit_seconds)
        .map_err(|_| PodObservationsError::LightTimeNotConverged)?;

    Ok(qtty::Meter::new(
        range_km(&receiver, &transmitter_at_emit) * 1_000.0,
    ))
}

fn range_km<C, F>(a: &TabulatedCartesianState<C, F>, b: &TabulatedCartesianState<C, F>) -> f64
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let displacement = a.position - b.position;
    displacement.magnitude().to::<Meter>().value() / 1_000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::Heliocentric;
    use crate::coordinates::frames::EME2000;
    use crate::pod::tabulated::{
        read_oem_tdb_ephemeris, TabulatedEphemeris, TabulatedEphemerisProvider,
    };

    const OEM_A: &str = "\
CCSDS_OEM_VERS = 2.0\n\
CREATION_DATE  = 2024-01-01T00:00:00\n\
ORIGINATOR     = TEST\n\
\n\
META_START\n\
OBJECT_NAME          = A\n\
OBJECT_ID            = -1001\n\
CENTER_NAME          = SUN\n\
REF_FRAME            = EME2000\n\
TIME_SYSTEM          = TDB\n\
START_TIME           = 2036-02-12T12:00:00.000\n\
STOP_TIME            = 2036-02-12T13:00:00.000\n\
META_STOP\n\
\n\
2036-02-12T12:00:00.000  0.0  0.0  0.0  0.0  0.0  0.0\n\
2036-02-12T13:00:00.000  0.0  0.0  0.0  0.0  0.0  0.0\n";

    const OEM_B: &str = "\
CCSDS_OEM_VERS = 2.0\n\
CREATION_DATE  = 2024-01-01T00:00:00\n\
ORIGINATOR     = TEST\n\
\n\
META_START\n\
OBJECT_NAME          = B\n\
OBJECT_ID            = -1002\n\
CENTER_NAME          = SUN\n\
REF_FRAME            = EME2000\n\
TIME_SYSTEM          = TDB\n\
START_TIME           = 2036-02-12T12:00:00.000\n\
STOP_TIME            = 2036-02-12T13:00:00.000\n\
META_STOP\n\
\n\
2036-02-12T12:00:00.000  1000.0  0.0  0.0  0.0  0.0  0.0\n\
2036-02-12T13:00:00.000  1000.0  0.0  0.0  0.0  0.0  0.0\n";

    type LisaEphemeris = TabulatedEphemeris<Heliocentric, EME2000>;
    type LisaProvider = TabulatedEphemerisProvider<Heliocentric, EME2000>;

    #[test]
    fn one_way_range_matches_static_separation() {
        let a: LisaEphemeris = read_oem_tdb_ephemeris(-1001, OEM_A.as_bytes()).unwrap();
        let b: LisaEphemeris = read_oem_tdb_ephemeris(-1002, OEM_B.as_bytes()).unwrap();
        let provider = LisaProvider::new(vec![a, b]).unwrap();
        let (coverage_start, _) = provider.ephemeris(-1001).unwrap().coverage().unwrap();
        let epoch = Time::<TDB>::from_raw_j2000_seconds(qtty::Second::new(
            coverage_start.to::<tempoch::J2000s>().raw().value() + 10.0,
        ))
        .unwrap();
        let range = one_way_light_time_range(&provider, -1001, -1002, epoch).unwrap();
        assert!((range.value() - 1_000_000.0).abs() < 1e-6);
    }
}
