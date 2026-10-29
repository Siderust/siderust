# DE4xx lunar surface geometry

Enable `spice` and use `spice::LunarObserver` for geometric Sun/Earth centre
observations from a site in a **runtime NAIF lunar frame**. One immutable
`SpiceContext` supplies both SPK states and binary-PCK orientation. Nothing in
this path calls VSOP87, ELP2000 or the analytical body orientation providers.
It uses the existing SPK segment evaluator, rather than converting through
`RuntimeEphemeris`'s TT/ecliptic convenience interface. It reuses `SurfaceSite`,
ENU and `HorizontalObservation`; runtime frame IDs remain attached to results.

## Kernel configuration

The tested set is:

| Kernel | Function |
| --- | --- |
| `de440s.bsp` | J2000 Sun, EMB, Earth and Moon state chains |
| `moon_pa_de440_200625.bpc` | Type 2 lunar orientation, class ID 31008, base J2000 |
| `moon_de440_220930.tf` | PA/ME definitions and fixed frame offsets |
| `pck00011.tpc` | Physical constants and separate analytical orientation metadata |
| `naif0012.tls` | Leap-second metadata for SPICE consumers / reference generation |

Download these from NAIF (large files stay outside Git):

```sh
mkdir -p /tmp/lunar-kernels
cd /tmp/lunar-kernels
base=https://naif.jpl.nasa.gov/pub/naif/generic_kernels
curl -fLO "$base/spk/planets/de440s.bsp"
curl -fLO "$base/pck/moon_pa_de440_200625.bpc"
curl -fLO "$base/fk/satellites/a_old_versions/moon_de440_220930.tf"
curl -fLO "$base/pck/pck00011.tpc"
curl -fLO "$base/lsk/naif0012.tls"
```

The 2022 FK is archived; NAIF also distributes a newer DE440 FK. Frame
relationships are read from the supplied FK, never inferred from filenames.
The reference manifest in `tests/fixtures/lunar_de440/provenance.json` records
SHA-256 hashes. The example loads the files in the table's order. Later files
and later segments have priority when their coverage overlaps. Unsupported
selected segments return an error, rather than exposing older data. A newer
segment outside its coverage does not hide an older covering segment.

Kernel compatibility is an explicit caller responsibility: SPICE descriptors
do not authenticate a DE release. Use a documented, matching SPK/BPC/FK set;
do not combine DE441 states with DE440 orientation without assessing that
model change. The API cannot silently select a different ephemeris or analytical
orientation; its context is borrowed immutably for the observer's lifetime.
`kernels` records all loaded candidate aliases (including order within each
type), not cryptographic identity or a claim that every loaded file contributed
to a particular result. Retain the file hashes with scientific output.

## Frames and direction

For the archived DE440 FK:

```text
J2000 -- time-dependent Type 2 binary PCK --> MOON_PA_DE440 (31008)
MOON_PA_DE440 -- fixed FK rotation --> MOON_ME_DE440_ME421 (31009)
MOON_PA (31010) aliases MOON_PA_DE440
MOON_ME (31011) aliases MOON_ME_DE440_ME421
```

`MOON_ME_DE440_ME421` is deliberately aligned with the DE421 ME realization;
it is not an analytical IAU model. The FK uses 67.8526, 78.6944, 0.2785
arcseconds about axes 3, 2, 1. TK matrices map the defined frame **to its
relative frame**; the ordered coordinate rotations multiply as R3 R2 R1.
Binary PCK returns the base-to-fixed rotation R3(W) R1(delta) R3(phi).
`rotation_naif(from, to, et)` returns a matrix acting on column position
vectors. Its inverse is the transpose. All rotations use the observation ET.

`MoonIauMeanEarth` (including legacy `MoonPrincipalAxes`) retains its existing
analytical `IAU_MOON` meaning. It is never used to type a runtime DE frame.
`LunarObserver::new` checks Moon-centred FK chains and requires binary PCK,
with J2000 base, at the terminal class. Unknown frames, cycles, missing binary
data and unavailable coverage fail explicitly. It rejects an analytical-only
configuration even if text PCK has been loaded.

## Rust workflow and site conventions

Run the complete example from the repository root:

```sh
cargo run --features spice --example 24_spice_lunar -- /tmp/lunar-kernels
```

The core API is:

```rust,no_run
use siderust::bodies::solar_system::Moon;
use siderust::coordinates::surface::SurfaceSite;
use siderust::qtty::{DEG, M};
use siderust::spice::{SpiceContext, LunarObserver, LunarTarget};
use siderust::tempoch::{JulianDate, TDB, TT};

let mut context = SpiceContext::new();
context.load("de440s.bsp")?;
context.load("moon_pa_de440_200625.bpc")?;
context.load("moon_de440_220930.tf")?;
let site = SurfaceSite::<Moon>::new(32.3 * DEG, -89.0 * DEG, 600.0 * M)?;
let observer = LunarObserver::new(&context, "MOON_ME", site)?;
let epoch = JulianDate::<TT>::new(2461358.5).to::<TDB>().to_j2000s();
let sun = observer.observe(LunarTarget::Sun, epoch)?;
let earth = observer.observe(LunarTarget::Earth, epoch)?;
// sun.horizontal.alt(), .az(), .range(); earth.above_horizon()
# Ok::<(), Box<dyn std::error::Error>>(())
```

`SurfaceSite` supplies only the site's geometry here; its longitude/latitude
are **interpreted in the selected runtime frame**. Its analytical `FixedFrame`
is not used. Longitude is east-positive. Default Moon latitude is
planetocentric and default datum is the 1737.4 km mean sphere; height is along
the outward sphere normal in metres. Explicit `with_surface` sites can use
planetographic latitude and an oblate ellipsoid: height/local up then follow
the reference-surface normal, as in the existing surface API. The scalar
latitude convention is preserved; planetocentric and planetographic angles
are equal on a sphere.

For actual radial distances (e.g. from a terrain product), use:

```rust,ignore
let observer = LunarObserver::radial(
    &context, "MOON_ME", 32.3 * DEG, -89.0 * DEG, 1738000.0 * M,
)?;
```

This uses the actual centre distance, zero height on a sphere of that radius,
and radial up. It does not invent an ellipsoidal height or a terrain normal.
Azimuth is north through east, in [0,360) degrees; altitude is geometric above
the selected local horizon. At poles the longitude defines the local meridian;
at zenith/nadir the existing horizontal singularity convention reports zero
azimuth. `above_horizon()` means **target centre altitude > 0**, excluding disc
extent and terrain masking. Finite-distance parallax is retained. Earth/Moon
states accumulate at their common EMB centre, avoiding AU-scale cancellation.

## Time and corrections

`observe` requires `Time<TDB>`, in seconds past J2000 when evaluated (SPICE
ET). Passing TT directly is a type error. `JulianDate<TT>::to::<TDB>()` makes
conversion explicit. For civil UTC, construct `Time<UTC>` using tempoch's
`try_from_chrono_with` and a `TimeContext`, then call `.to::<TDB>()`.
UTC interpretation uses tempoch's historical time data, including leap seconds;
there is no hardcoded contemporary UTC/TAI offset in this pipeline. Loading
an LSK into `SpiceContext` does **not** replace tempoch's time data. Ensure the
civil time source and LSK are consistent when comparing implementations.

The time regression covers 1972, 1999, the 2016 leap second, 2017 and 2026.
TT agrees with CSPICE within 1 microsecond; the separately tested tempoch and
LSK DELTET TDB approximations agree within 50 microseconds. For stricter
same-epoch tests, pass identical ET as `Time::<TDB>::new(et)`. The geometry
fixtures do exactly that; time-conversion approximation error is not folded
into orientation tolerances. The legacy `tdb_minus_utc` context helper remains
an approximate offset diagnostic and is not used for observations.

All results are geometric (`abcorr="NONE"`): no light time, stellar aberration,
refraction, relativistic deflection or terrain masks. Do not compare these
numbers directly with apparent Horizons (`LT+S`) results.

## Supported operations and failure modes

Binary PCK supports Type 2 angle polynomials with finite coefficients,
validated directory/record metadata and inclusive segment coverage, in
little- or big-endian IEEE DAF. Type 3 and Type 20 are indexed but fail if
selected. No derivatives/state-frame transformations are exposed. The tested
lunar BPC contains Type 2, class 31008 in J2000, covering ET
[-14200747200,13447252800]; the SPK normally limits observations first.

FK supports proper orthonormal MATRIX and three-axis ANGLES definitions in
radians, degrees or arcseconds. Dynamic/two-vector specifications, CK/SCLK
composition, other inertial bases, DSK, EK, meta-kernel expansion and general
kernel-pool redefinition semantics are outside this API. The context state
path requires J2000 SPK edges and fails on other selected frames. It uses
supported existing SPK evaluators (DE440s uses Type 2). Missing/disconnected
states, centre cycles, out-of-coverage epochs and unsupported selected data
produce errors; no Earth/Moon mass-ratio or analytical fallback is used.

DAF record traversal is bounded by file length and rejects cycles, oversized
summary counts, invalid links and byte-order inconsistencies. BPC validates
addresses, finite ordered coverage, integral directory sizes/counts, checked
record-size products, midpoint/radius consistency and finite coefficients.
FK rejects improper/nonorthogonal/nonfinite matrices. No new unsafe code,
global kernel state, runtime CSPICE dependency or runtime dependency is added.
The low-level BPC parser works with `no_std + alloc`; the high-level `spice`
feature requires `std` as before.

## Independent CSPICE validation

`tests/spice_lunar.rs` uses 282 independently generated observations over
32 days at equator, mid-latitudes, both poles, near the south pole, a
CLPS-representative site, and both sides of the ±180-degree boundary. Sun
horizon crossings have samples one second above/below. PA, ME and a synthetic 30/50/20-degree FK frame,
rotation inverses/determinants, interpolation boundaries, temporal continuity,
file/segment precedence and invalid data are checked. Azimuth differences are
circular; there are no singular zenith/nadir azimuth comparisons.

Only 16 KiB of original SPK/BPC coefficient excerpts are committed, alongside
the FK, LSK, references and hashes. CSPICE `SPKSUB` and `PCKW02` preserve the
original coefficients. An extra BPC record preserves the original kernel's
selection at the excerpt's final endpoint. Reference values query the original
full kernels through SpiceyPy 8.3.0 / CSPICE N0067, never Siderust. CSPICE
`pxform`, `spkpos(...,"NONE",...)` and `latrec` supply matrices, positions and
sites; independent NumPy ENU dot products supply angles/ranges/centre visibility.

Observed maximum residuals on Linux x86-64 (both excerpts and full SPK/BPC):

| Quantity | Maximum residual | Regression bound |
| --- | --- | --- |
| Rotation matrix element | 2.23e-16 | 5e-14 |
| Fixed target/site-relative component | 2.98e-8 km (0.030 mm) | 2e-7 km |
| Slant range | 5.96e-8 km (0.060 mm) | 2e-7 km |
| Altitude/circular azimuth | 5.69e-14 degree | 1e-10 degree |

Bounds allow f64 evaluation/platform rounding, not model uncertainty. In
particular, angle reduction follows NAIF PCKE02's signed D_MOD subtraction.
These results establish same-kernel numerical agreement, not absolute lunar
pointing accuracy or mission certification. Full-BPC coverage-edge rotation
checks have a 5e-12 element bound, including epochs centuries from J2000.
The downstream 0.02-degree Horizons benchmark remains a separate comparison
requiring matched corrections, epochs, sites and frame models.

Regenerate offline reference data (development only; Rust builds/tests do not
require Python or CSPICE):

```sh
uv venv /tmp/lunar-reference-env
uv pip install --python /tmp/lunar-reference-env/bin/python spiceypy==8.3.0 numpy==2.5.4
KERNEL_DIR=/tmp/lunar-kernels /tmp/lunar-reference-env/bin/python \
  tools/cspice_refs/cspice_lunar_de440_references.py
cargo test --features spice --test spice_lunar
SIDERUST_LUNAR_KERNEL_DIR=/tmp/lunar-kernels \
  cargo test --features spice --test spice_lunar cspice_external -- --ignored --nocapture
```

The optional external test spans approximately 1850–2150 and tests the original
BPC's coverage endpoints. No large kernels or network access are needed for
ordinary CI.
