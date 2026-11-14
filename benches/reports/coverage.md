# Workload coverage and scientific validation

The [catalog](../README.md#catalog) gives features and datasets for each target.
The cases below measure specific entry points. The linked tests and reference
data define numerical acceptance for their fixtures and model configurations.

## Coordinates and time

| Suite | Timed boundary | Validation and open coverage |
| --- | --- | --- |
| [converters](../02_converters.rs) | Convert a prebuilt barycentric ICRS position between Cartesian and spherical representations with AU distance. Each direction uses its own constant input. | [Coordinate tests](../../tests/test_coordinates.rs) cover coordinate behavior. Representation conversion leaves frame and center unchanged. Frame transforms, precession/nutation matrices and batch conversion lack isolated timings. |
| [sidereal](../08_sidereal.rs) | Evaluate GMST and GAST at JD TT 2460000.5 with precomputed UT1. The low-level GAST case also receives precomputed nutation. Despite its name, the current `gast_iau2006a` wrapper computes IAU 2000B internally. | See [sidereal unit tests](../../src/astro/sidereal.rs) and [SOFA-based orientation regression tests](../../tests/test_high_precision_earth_rotation_regression.rs). The observed difference between the GAST cases includes nutation evaluation. Time-scale conversion, EOP lookup and other nutation policies need separate workloads. |

## Ephemerides

| Suite | Timed boundary | Validation and open coverage |
| --- | --- | --- |
| [vsop87](../12_vsop87.rs) | One VSOP87A or E evaluation via a prebuilt trait object, for Mercury through Neptune at J2000 TT + 9,131.25 days. The A and E series have different origins and coefficient sets. | [Ephemeris tests](../../tests/test_ephemeris.rs) and [VSOP87E velocity checks](../../tests/test_vsop87e_vel.rs) provide numerical coverage. Batch throughput and epoch sensitivity are open. |
| [elp2000](../04_elp2000.rs) | Public Moon position in km at J2000 TT and offsets of +1,234.567, +18,262.5 and −3,652.5 days. | [Frame consistency tests](../../tests/ephemeris_frame_consistency.rs) exercise lunar conventions and external references. The four epochs sample evaluator cost. |
| [ephemeris_comparison](../05_ephemeris_comparison.rs) | The full backend's Sun barycentric, Earth heliocentric, Earth barycentric velocity and Moon geocentric queries at J2000 TT + 9,131.25 days. With `embedded-compact`, the three position groups also call the explicit compact backend at that same epoch. | Both paths expose canonical ecliptic positions with matching units and origins. The compact model truncates coefficients and checks TDB coverage. [Compact fixture tests](../../tests/compact_ephemeris.rs) and the [DE440s/ERFA validation workflow](../../tools/validation/README.md) define sampled accuracy criteria. The timing cases use one interior epoch and leave range edges, repeated epoch batches and the compact finite-difference velocity path open. |

```sh
cargo bench --locked -p siderust --bench ephemeris_comparison --features embedded-compact -- \
  'ephemeris/moon_geocentric'
cargo test --locked -p siderust --test compact_ephemeris --no-default-features --features embedded-compact
```

Host Criterion cases with both features enable full and compact tables in one
build. For flash size, RAM, stack high-water and device latency, measure a
compact-only application on its target. See the [embedded guide](../../doc/guides/embedded.md)
and the [ESP32 experiment](compact_ephemeris.md#esp32).

## Observation events

All four suites use Roque de los Muchachos. Window workloads start at midnight
UTC on 2026-01-01, converted to the library's MJD/TT time type before timing.
Single evaluations use MJD 51544.5 TT. The harnesses use the library
default time and Earth orientation policy and load no external EOP bundle. Each search includes dispatch, any internal
context construction, crossing discovery and result allocation/destruction.
The thresholds are those of the altitude API. See the
[observation model and examples](../../src/event/altitude/mod.rs) before comparing
with external tables that apply different refraction or horizon conventions.

| Suite | Workloads and models | Validation and open coverage |
| --- | --- | --- |
| [altitude_comparison](../01_altitude_comparison.rs) | Single Sun/Moon/Sirius evaluations; 7- and 365-day above 0°; 30-day below −18° for Sun and 0° for Moon/Sirius. | These are body-specific workloads with different evaluation chains and crossing counts. [Altitude API tests](../../tests/test_altitude_api.rs) cover dispatch and semantics. The 30-day cases also differ in threshold. |
| [solar_altitude](../09_solar_altitude.rs) | 30/184/365 days, daily predictor and internal alternatives. | [Solar report](solar_altitude_perf_report.md) links search and external checks. Polar, grazing and narrow-window timing remain open. |
| [moon_altitude](../06_moon_altitude.rs) | 30/184/365 days, cached lunar search and forced scan. | [Lunar report](moon_altitude_perf_report.md) identifies the shared model and cache. Direct uncached reference timing, reusable-query timing and a validated compact-only event corpus remain open. |
| [star_altitude](../11_star_altitude.rs) | Fixed ICRS Sirius (RA 101.287°, Dec −16.716°), 30/184/365-day above/below 0°, 10°–50° range and horizon crossings, plus 30/365-day above-horizon scan comparison. | [Stellar engine tests](../../src/event/stellar/altitude_periods.rs) compare analytical and scan paths for sampled sites/windows and cover high latitude and grazing. Moving targets, azimuth events and combined constraints remain open. |

## SPK and satellite propagation

| Suite | Timed boundary | Validation and open coverage |
| --- | --- | --- |
| [de441](../03_de441.rs) | Five loaded `RuntimeEphemeris` queries at J2000 TT, plus GMST/GAST at that epoch with precomputed UT1. BSP and context loading are outside timing. | [Runtime SPK tests](../../tests/test_spk_runtime.rs) and [CSPICE/ERFA frame tests](../../tests/ephemeris_frame_consistency.rs) cover segment selection, coverage errors and frame conversion. Fixed-epoch latency favors record locality. |
| [spk_perf](../10_spk_perf.rs) | Synthetic Type 2 Earth position/velocity and Moon position, Type 3 Sun position at J2000 + 500 days, optional parsing of resident BSP bytes including output destruction, and preloaded J2000 queries with `RuntimeEphemeris` or `SpkKernelSet`. | Synthetic data have one record and two coefficients per component. Type 2 and Type 3 cases use different outputs and chain lengths. Real parser cost depends on kernel structure. Cold file loading, epoch sweeps, multi-kernel resolution, binary PCK and lunar orientation are open. See [lunar SPICE reference tests](../../tests/spice_lunar.rs). |
| [sgp4](../07_sgp4.rs) | SGP4/SDP4 with WGS-72/AFSPC for Vallado TLEs 06251 (LEO), 24208 (GEO), 08195 (half-day resonance). Parsing occurs before timing. Initialization uses the parsed TLE. Single propagation uses +1440.25 minutes. The repeated case produces 100 TEME states at offsets `i × 60.25` minutes using the same propagator. | [Vallado verification](../../tests/sgp4_vallado.rs) and [fixture conventions](../../tests/fixtures/sgp4/README.md) describe the WGS-72 reference and tolerances. Catalog ingestion, long-term sweeps, and TEME-to-observer pipelines remain open. |

The highest-value extensions are a realistic bounded SPK epoch trace, an
end-to-end frame and time transform with pinned EOP inputs, and controlled
compact device latency and peak-memory measurements. POD, mission dynamics,
atmospheric models and other public modules have no registered Criterion suites.
