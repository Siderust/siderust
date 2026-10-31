# Frames and Centers Reference

This document is the authoritative guide to the **reference frames** and
**reference centers** used in siderust.  Read it before touching any coordinate
transformation or adding a new frame type to `affn`.

---

## 1. Reference Centers

A *center* is the origin of a coordinate system.  In siderust every
`Position` carries both a center and a frame.

| Center type | Struct | Where origin lives |
|---|---|---|
| `Barycentric` | `centers::Barycentric` | Solar-System barycentre |
| `Heliocentric` | `centers::Heliocentric` | Sun's centre of mass |
| `Geocentric` | `centers::Geocentric` | Earth's centre of mass |
| `Topocentric` | `centers::Topocentric` | Observer's location on Earth's surface |
| `Mercurycentric` | `centers::Mercurycentric` | Mercury's centre of mass |
| `Venuscentric` | `centers::Venuscentric` | Venus's centre of mass |
| `Marscentric` | `centers::Marscentric` | Mars's centre of mass |
| `Selenocentric` | `centers::Selenocentric` | Moon's centre of mass |
| `Jovicentric` | `centers::Jovicentric` | Jupiter's centre of mass |
| `Saturnocentric` | `centers::Saturnocentric` | Saturn's centre of mass |
| `Uranocentric` | `centers::Uranocentric` | Uranus's centre of mass |
| `Neptunocentric` | `centers::Neptunocentric` | Neptune's centre of mass |
| `Plutocentric` | `centers::Plutocentric` | Pluto's centre of mass |

`Topocentric` is *parameterised*: it carries a `Geodetic<ECEF>` so horizontal
coordinates know their own observation point without external context.

### Center-Shift Hub

All center-shift transformations route through the **Barycentric** hub:

```
Planetocentric → Barycentric → target center
```

Mercury–Neptune translations request the selected `Ephemeris` backend's
`try_planet_barycentric` capability. The default VSOP87 backend supplies
normalized VSOP87E positions. Pluto requests `try_pluto_barycentric`; the
default backend retains its coarse catalogue Keplerian orbit plus its Sun.
Other backends must opt in explicitly. Unsupported capabilities never select
VSOP87 or catalogue data implicitly. Selenocentric translations use the same
backend's `try_moon_geocentric` plus `try_earth_barycentric`.

`CenterShiftProvider::try_shift` returns `CenterShiftError`, preserving the
planetary-capability, ephemeris-evaluation or optional Lagrange-solver cause.
Reverse and composed shifts propagate the same errors using one provider path.
For positions with `Params = ()`, use the fallible coordinate entry point:

```rust
use siderust::coordinates::{cartesian::Position, centers::{Barycentric, Marscentric},
    frames::EclipticMeanJ2000, transform::{AstroContext, try_to_center_with}};
use siderust::qtty::AstronomicalUnit;
let context: AstroContext = AstroContext::default();
let input = Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(1., 0., 0.);
let result: Result<Position<Marscentric, EclipticMeanJ2000, AstronomicalUnit>, _> =
    try_to_center_with(&input, siderust::J2000, &context);
```

`shift`, `center_shift*`, `TransformCenter::to_center*` and combined `to*`
entry points remain infallible conveniences and panic when a required shift
is unavailable. Use `try_shift` or `try_to_center_with` when capability/coverage
must be checked. Custom `CenterShiftProvider` implementations implement
`try_shift` in 0.18; the default `shift` delegates to it. Custom Pluto backends
implement `try_pluto_barycentric`, not only the infallible wrapper.
`PlanetState<Body>` retains its heliocentric capability contract; implementors
may compose that state with their own Sun to supply a planetary-center shift.
`DynEphemeris`/`RuntimeEphemeris` retain their object-safe Sun/Earth/Moon API;
this does not add a runtime planetary-center transformation engine.

---

## 2. Earth-Fixed Frames

### 2.1 ECEF, Earth-Centred Earth-Fixed (WGS84 Earth-fixed convenience frame)

```rust
// affn::frames::ECEF
pub struct ECEF;
```

`ECEF` is the WGS84-tagged Earth-fixed frame used for geodetic site input and
for Earth-fixed Cartesian vectors in user-facing APIs.

**When to use `ECEF`:**
- Storing WGS84 geodetic coordinates (`Geodetic<ECEF>`).
- User-facing Earth-fixed Cartesian vectors when you want the familiar `ECEF`
  label instead of `ITRF`.

**Provider-layer note:**  frame rotations currently treat `ECEF` as
co-aligned with `ITRF`.  The distinction between the two frames is therefore
ellipsoidal / metadata-oriented in siderust today:
- `Geodetic<ECEF>` uses **WGS84**
- `Geodetic<ITRF>` uses **GRS80**
- `Position<Geocentric, ECEF, U>` and `Position<Geocentric, ITRF, U>` share
  the same Earth-orientation chain when rotated to inertial frames

### 2.2 ITRF, International Terrestrial Reference Frame (EOP-realised)

```rust
// affn::frames::ITRF
pub struct ITRF;
```

`ITRF` is the **physical** geocentric frame co-rotating with the solid Earth.
Its axes are realised through a global network of VLBI / SLR / GNSS stations.
The full IERS EOP chain applies:

```
ITRS → W⁻¹ → TIRS → ERA → CIRS → Q → GCRS/ICRS
```

**When to use `ITRF`:**
- Observatory geocentric coordinates sourced from ITRF2020 / VLBI solutions.
- Polar-motion-corrected baselines at the centimetre level.
- Any context where the phrase "ITRF coordinates" appears in the source data.

### 2.3 CIRS / TIRS, Intermediate Earth-orientation frames

```rust
pub struct CIRS;
pub struct TIRS;
```

These are the public intermediate frames of the IAU 2000/2006 Earth-rotation
chain:

```text
GCRS/ICRS ──Q──▶ CIRS ──ERA──▶ TIRS ──W──▶ ITRF/ECEF
```

Use them when you need to inspect or exchange states at a specific stage of the
Earth-orientation reduction instead of jumping directly between inertial and
Earth-fixed frames.

---

## 3. Celestial (Inertial) Frames

### 3.1 ICRS / ICRF

```rust
pub struct ICRS;
pub struct ICRF;
```

The **International Celestial Reference System** is the primary inertial frame.
Axes are fixed by extragalactic radio sources (quasars).  It is realised by
the ICRF catalogue.

ICRS is the system; ICRF is its observational realization. Both markers
support spherical and Cartesian coordinates. Siderust treats their axes as
identical to the accuracy of the source realization.

### 3.2 GCRS, Geocentric Celestial Reference System

```rust
pub struct GCRS;
```

For directions, Siderust treats GCRS axes as aligned with ICRS. This identity
is an orientation approximation, not a complete relativistic BCRS/GCRS
coordinate transformation. The IAU J2000 frame bias belongs between ICRS/GCRS
axes and the dynamical J2000 mean equator; do not apply it between GCRS and ICRS.

### 3.3 EquatorialMeanJ2000 / EME2000

```rust
pub struct EquatorialMeanJ2000;
pub struct EME2000;
```

IAU 2006 mean equator and equinox at J2000.0 (JD 2451545.0 TT).
The ICRS-to-frame mapping is SOFA/ERFA `bp06(2451545.0, 0.0)`'s fixed bias
matrix `rb` (about 23 mas). The public frame-bias extensions and transform
providers share this matrix. It is distinct from the empirical orientation of
the FK5 star catalogue returned by `fk5hip`; catalogue spin is not applied to
fixed J2000 axes.

`EME2000` is provided as an explicit CCSDS-facing synonym of the same axes.
Its transform to/from `EquatorialMeanJ2000` is the identity, but keeping the
name separate preserves frame provenance in exchanged data.

### 3.4 FK4 B1950

```rust
pub struct FK4B1950;
```

The **Fourth Fundamental Catalogue** equatorial frame, referred to the mean
equator and equinox of the Besselian epoch B1950.0.  The FK4 → ICRS
transformation uses the Standish (1982) frame-tie rotation matrix.

**When to use `FK4B1950`:**
- Converting positions from legacy FK4-epoch catalogues.
- Cross-matching with historical observing records.

**Accuracy note:**  The FK4 → FK5 (ICRS) transformation involves a frame
rotation, an equinox correction, and an E-term removal.  siderust implements
the constant rotation matrix; E-term removal is not yet included.

### 3.5 TEME, True Equator, Mean Equinox

```rust
pub struct TEME;
```

The **True Equator, Mean Equinox** frame is the de facto reference for
SGP4/SDP4 (TLE) satellite orbit propagation.  The TEME-to-TOD transformation
is a single rotation about the z-axis by the equation of the equinoxes:

```
TEME → Rz(Δψ · cos εA) → TOD (EquatorialTrueOfDate)
```

**When to use `TEME`:**
- Propagating TLE elements with SGP4/SDP4.
- Importing NORAD ephemerides.

### 3.6 Galactic

```rust
pub struct Galactic;
```

The **Galactic coordinate system** (IAU 1958), centred on the Sun, with the
fundamental plane parallel to the Galactic plane.  Galactic longitude *l*
increases towards the Galactic centre; Galactic latitude *b* is positive
towards the North Galactic Pole.

The Galactic → ICRS rotation is a pre-computed constant matrix derived from
the North Galactic Pole coordinates (α = 192.85948°, δ = 27.12825°) and the
Galactic centre position angle (l₀ = 32.93192°), per Murray (1989).

### 3.7 Ecliptic Frames

| Type | Description |
|---|---|
| `EclipticMeanJ2000` | IAU 2006 mean ecliptic/equinox at J2000 (`ecm06`) |
| `EclipticMeanOfDate` | Mean ecliptic and equinox of observation epoch |
| `EclipticTrueOfDate` | True ecliptic (adds nutation in obliquity and longitude) |

#### Canonical ecliptic and ephemeris boundaries

`EclipticMeanJ2000` has one backend-independent meaning: the IAU 2006 mean
ecliptic and equinox at J2000. For column vectors its ICRS mapping is
`C = Rx(-84381.406 arcsec) * rb`, exactly the convention of SOFA/ERFA
`ecm06(2451545.0, 0.0)`. It differs from SPICE `ECLIPJ2000` and from the
native dynamical ecliptics of the analytical theories.

| Source | Native axes | Normalization before public coordinate construction |
|---|---|---|
| VSOP87A/E | Inertial dynamical ecliptic/equinox J2000 | IMCCE published native→FK5 matrix `V`, fixed J2000 FK5→ICRS tie `H`, then `C` |
| ELP2000-82B | Moving inertial mean ecliptic and departure point | Existing Laskar rotation to its fixed native J2000 triad, ELP-specific native→FK5 matrix `L`, then `H` and `C` |
| Modern ICRF-realized DE kernels (e.g. DE440/441) | ICRF, labeled NAIF `J2000` (ID 1) | Existing `C`, after physical target/center-chain arithmetic |

`V` comes from IMCCE `vsop87.doc`, REFERENCE SYSTEM. `L` comes independently
from the ELP82B notice §8 p.12: obliquity 84381.40883 arcsec and equinox arc
-0.09845 arcsec. These native realizations are not assumed identical. The
Laskar step is retained exactly; it is not replaced by ordinary equinox-of-date
precession or repeated during normalization.

`H` is the fixed J2000 orientation from SOFA/ERFA `fk5hip` (orientation
components -19.9, -9.1, +22.9 mas), using Hipparcos as an optical realization
of ICRS. This conventional empirical tie does not remove catalogue zonal errors
or recover an exact DE200-to-DE440 alignment. The native published ties and
analytical trajectories retain their scientific uncertainties. All rotations
here are static: positions and VSOP velocities use the same `R`; there is no
catalogue-spin or `Rdot * position` term. The published IMCCE matrices have
12-decimal precision (orthogonality errors below 1e-12).

Public `vsop87a/e`, velocity and combined-state accessors normalize at their
source. `Vsop87Ephemeris` delegates to them without another correction.
`Moon::get_geo_position` converts TT input to TDB, evaluates ELP and its Laskar
rotation, then normalizes once. The legacy `Moon::vsop87a` accessor evaluates
the VSOP Earth–Moon barycenter coefficients, not the geocentric lunar state.

A shared `Position<Center, Frame, Unit>` describes physical coordinates, not
a calculation method. Different theories can estimate different trajectories
in the same axes; frame agreement is not ephemeris-accuracy agreement. Model
comparisons must match the physical target **and** reference center: e.g.
Mars `499` is not automatically system barycenter `4`. If a kernel lacks the
planet-center offset, omit the accuracy comparison rather than substitute the
system barycenter. No transformation is fitted to modern-DE residuals.
Keplerian/conic inputs are defined in the supplied frame; do not apply a model
correction to user-supplied canonical elements. Coarse built-in orbital elements
and the abbreviated Pluto series do not provide a precision realization tie.

#### ICRF-only SPK loader policy

`RuntimeEphemeris`, `SpkKernelSet` and the legacy body-chain parser require
producer-documented ICRF realization, frame ID 1 and Type 2/3 **for every segment**.
They reject an entire mixed/incompatible kernel, even for unrelated targets,
before returning or mutating a usable kernel set. Strict rejection prevents a
skipped higher-priority incompatible segment from exposing stale data.
`try_from_indexed_segments` applies the same policy to parsed input and checks
finite, ordered coverage. `RuntimeEphemeris` additionally requires its fixed
fundamental chains: Sun/SSB, EMB/SSB, Moon/EMB and optional Earth/EMB; alternate
chains for those targets are rejected. `SpkKernelSet` resolves general supported
target/center chains. General frame-aware SPICE context operations are separate.

Later segments/kernels win **within declared coverage**. Once selected, an
invalid record or unavailable center chain is an error, not permission to use
an older overlapping segment. An older segment may be used where a newer one
has no coverage. Legacy `parse_bsp` cannot represent a stack and rejects
multiple segments for a fundamental target; use the indexed/stacked loaders. Checked constructors return actionable errors; retained
infallible parsed constructors panic on rejected metadata. Metadata-free legacy
`BspSegments` explicitly assumes ICRF, and its checked constructor validates types.

ID 1 alone cannot distinguish historical FK5/DE200 data from modern ICRF
realizations. The caller must establish provenance from the kernel producer;
these ICRF APIs do not convert arbitrary historical J2000 realizations. Modern
DE states receive no extra classical FK5 bias. Unit conversions retain the
existing coefficient AU convention and use 149597870.700 km/AU and 86400 s/day.

Authoritative references:
[VSOP87 conventions](https://ftp.imcce.fr/pub/ephem/planets/vsop87/vsop87.doc),
[Bretagnon & Francou 1988](https://articles.adsabs.harvard.edu/pdf/1988A%26A...202..309B),
[ELP82B notice §§8–9](https://ftp.imcce.fr/pub/ephem/moon/elp82b/elp82b.ps),
[ERFA ecm06](https://github.com/liberfa/erfa/blob/v2.0.1/src/ecm06.c),
[ERFA fk5hip](https://github.com/liberfa/erfa/blob/v2.0.1/src/fk5hip.c),
[NAIF ICRF/J2000 conventions](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#ICRF%20vs%20J2000),
and [SPK precedence](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Data%20Precedence).
Permanent independent reference values, units, time conventions and tolerances
are recorded in `tests/data/ephemeris_frames.json` and exercised through Cargo
by `tests/ephemeris_frame_consistency.rs`. The ignored same-kernel DE440s test
requires `SIDERUST_BSP_PATH=/path/to/de440s.bsp`; normal reference tests need
no external kernel, Python or astronomy package.

---

## 4. Horizontal Frame

```rust
pub struct Horizontal;
```

Local horizontal frame centred on the observer.  Spherical coordinates are
**(altitude, azimuth)**, altitude measured from the horizon (+90° = zenith),
azimuth measured from North through East.

Every `Position<Topocentric, Horizontal, U>` implicitly carries a
`Geodetic<ECEF>` via the `Topocentric` center parameter (Earth-only path).

### 4.1 Generic body-surface observers

For observers on bodies other than Earth — and for a uniform Earth API —
use `coordinates::surface`:

| Type | Role |
|---|---|
| `SurfaceSite<B>` | Lat/lon/height on body `B` with an explicit `ReferenceSurface` |
| `SurfaceObserver<B>` | Parameterized centre carrying `SurfaceObserverParams` |
| `BodySurfaceModel` | Compile-time link: body → centre, fixed frame, default surface |

**Pipeline** (geometric Alt/Az):

```text
target → body-centred (Barycentric hub) → body-fixed (ICRS hub)
      → observer-relative LOS → local ENU → Horizontal
```

Earth sidereal-time (GMST/GAST/UT1) is **not** used for lunar or planetary
sites. Body orientation follows the corresponding IAU/WGCCRE model; the
Moon uses the complete analytical IAU_MOON periodic terms. The underlying
lunar arguments use TDB, converted from TT by tempoch. Other planetary
providers retain their existing linear elements and numerical TT approximation;
they do not implement all periodic terms of the current text PCK.


Lunar numerical validation uses all three axes of CSPICE `pxform` with
`pck00011.tpc` and `naif0012.tls`, at five TT epochs spanning 1900--2050.
Direction cosines agree within `1e-10` after TT-to-TDB conversion. This tests
implementation agreement, not absolute lunar pointing accuracy. NAIF describes
the analytical model's discrepancy from DE421 ME over about 2000--2020
as having an amplitude of 0.0051 degrees (155 m on the lunar surface);
this is not a universal error bound. Use binary PCKs and their
matching DE frame kernels for precision lunar cartography. The tests do not
establish a physical accuracy bound outside that epoch span.

Sources: [NAIF text PCK](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc),
[lunar frame definitions](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf),
[WGCCRE 2009](https://doi.org/10.1007/s10569-010-9320-4),
[WGCCRE 2015](https://doi.org/10.1007/s10569-017-9805-5).
`Moon::ROTATION` / `IauRotationParams` retain the secular linear elements for
compatibility; they are not the complete lunar orientation. Use the frame
provider for periodic lunar orientation.

**Support matrix**

| Body | Surface default | Fixed frame | Notes |
|---|---|---|---|
| Earth | WGS84 ellipsoid | `ECEF` | Coexists with `Topocentric` / GAST APIs; pgr lon east-positive |
| Moon | IAU sphere 1737.4 km | `MoonIauMeanEarth` (distinct from `MoonPrincipalAxes`) | NAIF IAU_MOON analytical Mean-Earth, **not** JPL DE PA / binary-PCK ME; pgr lon east-positive |
| Mars | IAU ellipsoid 3396.19 km | `MarsFixed` | Default latitude planetographic; **SPICE pgr lon west-positive** — convert explicitly |
| Mercury, Venus, Pluto | IAU mean spheres | `*Fixed` | Mercury/Pluto pgr lon west-positive; Venus east-positive |
| Gas giants | *(no `BodySurfaceModel`)* | — | Implement the trait only with a documented reference surface |

**Longitude:** Siderust stores east-positive longitudes for every body.
Planetocentric longitude (SPICE `latrec`/`reclat`) matches this directly.
Planetographic longitude (SPICE `pgrrec`/`recpgr`) defaults to west-positive
for prograde rotators such as Mars and east-positive for Earth/Moon/Venus;
convert with `BodySurfaceModel::east_positive_longitude` or
`SurfaceSite::from_planetographic_longitude`. No hidden sign flips.
SPICE kernel-pool overrides (`BODY*_PGR_POSITIVE_LON`) are **not** read by
Siderust — if an external product used a non-default sense, pass that sense
explicitly to the conversion helpers.

Planetocentric and planetographic **latitudes** are distinct on ellipsoids and
must not be mixed silently. See `coordinates::surface` rustdocs,
`examples/22_moon_surface_quickstart.rs`, and
`examples/21_body_surface_horizontal.rs`.

**Offline reference regeneration** (development only; SpiceyPy optional):
`python3 tools/cspice_refs/cspice_surface_references.py`.
CSPICE supplies lunar `pxform` / Mars `pgrrec` and TT→TDB; ENU Alt/Az/range
literals are independent sphere geometry, not CSPICE `azlcpo`.

---

## 5. Planetary Body-Fixed Frames

These frames rotate with the respective body, using IAU pole and prime-meridian
orientation parameters (α₀, δ₀, W). The active, column-vector body-fixed to
ICRS rotation implemented by the provider is:

```
R = Rz(α₀ + 90°) · Rx(90° − δ₀) · Rz(W)
```

where α₀ and δ₀ define the body's north pole direction in ICRS, and W is
the prime meridian angle (degrees). This order matters: `affn`'s ZXZ constructor
applies its first Z rotation first, then X, then the second Z. For the Moon,
α₀, δ₀ and W additionally include the 13 periodic arguments of the NAIF
IAU_MOON text-PCK and W includes a quadratic term. `MoonIauMeanEarth` is the
analytical marker. It is distinct from `MoonPrincipalAxes`, which has no
analytical provider. Migrate old analytical uses of the principal-axes name to
`MoonIauMeanEarth`; actual principal-axes data must not be relabeled. DE-specific
`MOON_PA` and `MOON_ME` require suitable binary PCK and frame kernels. See
[NAIF lunar frame definitions](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf).

Most frame markers are zero-sized structs from `affn::frames` behind its
`astro` feature, re-exported through `siderust::coordinates::frames` and
`planetary`. The analytical `MoonIauMeanEarth` marker is defined in Siderust;
geometry remains in `affn`, astronomical rotations in Siderust.

The **IAU rotation parameters** remain in `siderust::coordinates::frames::planetary`.
The corresponding `FrameRotationProvider` implementations are organized under
`siderust::coordinates::transform::providers::frames_planetary`.

| Frame | Struct | Body | Spherical coords |
|---|---|---|---|
| Mercury body-fixed | `MercuryFixed` | Mercury | (lat, lon, radius) |
| Venus body-fixed | `VenusFixed` | Venus | (lat, lon, radius) |
| Mars body-fixed | `MarsFixed` | Mars | (lat, lon, radius) |
| Moon analytical ME | `MoonIauMeanEarth` | Moon | (lat, lon, radius) |
| Moon principal axes | `MoonPrincipalAxes` (no analytical provider) | Moon | (lat, lon, radius) |

| Jupiter System III | `JupiterSystemIII` | Jupiter | (lat, lon, radius) |
| Saturn body-fixed | `SaturnFixed` | Saturn | (lat, lon, radius) |
| Uranus body-fixed | `UranusFixed` | Uranus | (lat, lon, radius) |
| Neptune body-fixed | `NeptuneFixed` | Neptune | (lat, lon, radius) |
| Pluto body-fixed | `PlutoFixed` | Pluto | (lat, lon, radius) |

The coordinate Serde representation remains component data plus any center
parameters; `affn` does not embed a frame-name discriminator. The caller must
retain frame provenance and deserialize into the intended typed frame. Old
unlabeled data do not prove analytical/PA equivalence and must not be relabeled
without establishing their original convention.

### IAU Rotation Parameters

Each body's rotation is defined by three functions of time:

- **α₀(T)**, right ascension of the north pole (degrees), T in Julian centuries
- **δ₀(T)**, declination of the north pole (degrees), T in Julian centuries
- **W(d)**, prime meridian angle (degrees), d in days from J2000.0

These are stored as `IauRotationParams` constants (e.g., `MARS_ROTATION`)
in `coordinates::frames::planetary`.

---

## 6. The Earth-Rotation Chain

The canonical ITRS → EquatorialMeanJ2000 rotation is implemented in
`crate::astro::earth_rotation_provider::itrs_to_equatorial_mean_j2000_rotation`.
The same chain is also exposed through public `FrameRotationProvider`
implementations for `GCRS`, `CIRS`, `TIRS`, `ITRF`, and `ECEF`.

```
ITRS  ──W⁻¹──▶  TIRS  ──ERA──▶  CIRS  ──Q──▶  GCRS/ICRS  ──rb──▶  EquatorialMeanJ2000
```

| Step | Matrix | Source |
|---|---|---|
| **W** | Polar motion (xₚ, yₚ, s′) | EOP → `polar_motion_matrix_from_eop` |
| **ERA** | Earth Rotation Angle | UT1 → `earth_rotation_angle` |
| **Q** | CIO/CIP (X, Y, s) | selected nutation model + EOP dX,dY (IAU 2006A default) |
| **rb** | Fixed IAU 2006 ICRS → MeanJ2000 bias | `FrameRotationProvider<ICRS, EquatorialMeanJ2000>` |

Public provider coverage now includes:
- `GCRS ↔ CIRS`
- `CIRS ↔ TIRS`
- `TIRS ↔ ITRF`
- `ICRS ↔ CIRS/TIRS/ITRF/ECEF`
- `EclipticMeanJ2000`, `EquatorialMeanJ2000`, `EME2000`, `ICRF`, `TEME`,
  `Galactic`, and `FK4B1950` ↔ `CIRS/TIRS/ITRF/ECEF`

**Time-scale contract:**  The `jd` argument to `itrs_to_equatorial_mean_j2000_rotation`
is a **Terrestrial Time (TT)** Julian Date.

---

## 7. Frame-Rotation Hub

ICRS is the orientation hub; direct frame-pair implementations may bypass it:

```
F1 → ICRS → F2
```

The following constant rotation matrices are pre-computed:

| Rotation | Implementation | Reference |
|---|---|---|
| FK4 B1950 → ICRS | `FK4_TO_FK5` composed with the shared inverse bias | Standish (1982) + frame bias |
| Galactic → ICRS | `GALACTIC_TO_ICRS` | Murray (1989), NGP coordinates |
| ICRS → EquatorialMeanJ2000 | `bias::frame_bias_icrs_to_j2000()` | IAU 2006 frame bias |

Time-dependent rotations:

| Pair | Method |
|---|---|
| GCRS ↔ CIRS | CIO/CIP matrix **Q** |
| CIRS ↔ TIRS | Earth Rotation Angle |
| TIRS ↔ ITRF | Polar motion matrix **W** |
| TEME ↔ TOD | Rz(equation of equinoxes) |
| EquatorialMeanJ2000 ↔ EquatorialMeanOfDate | IAU 2006 precession with fixed bias removed |
| EquatorialMeanOfDate ↔ EquatorialTrueOfDate | selected nutation model (IAU 2006A default) |
| Body-fixed ↔ ICRS | body-specific IAU rotation model (Moon: periodic IAU 2009) |

Provider implementations are grouped by responsibility:

- `providers::frames_inertial`
- `providers::frames_earth`
- `providers::frames_catalog`
- `providers::frames_teme`
- `providers::frames_planetary`
- `providers::centers_standard`
- `providers::centers_planetary`

Shared helper utilities in `providers::mod` handle inversion, hub composition,
and rotation of ecliptic ephemeris vectors into the requested destination frame.

---

## 8. Geodetic Coordinates and Ellipsoid Conversion

### 8.1 `Geodetic<F, U>` (siderust) / `ellipsoidal::Position<C, F, U>` (affn)

```rust
pub type Geodetic<F, U = Meter> = affn::ellipsoidal::Position<Geocentric, F, U>;
```

`Geodetic<F>` is the first-class type for Earth-surface positions.  The ellipsoid
is encoded in the frame `F` via `HasEllipsoid`:

- `Geodetic<ECEF>`, WGS84
- `Geodetic<ITRF>`, GRS80

### 8.2 Converting to/from Cartesian ECEF

```rust
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::ECEF;
use qtty::*;

let coord = Geodetic::<ECEF>::new(-17.8925 * DEG, 28.7543 * DEG, 2396.0 * M);
let ecef = coord.to_cartesian::<Meter>();
// Round-trip
let back = Geodetic::<ECEF>::from_cartesian(&ecef);
```

Named getters (`lat()`, `lon()`, `altitude()`) are generated automatically by the
`#[frame(inherent, ellipsoid = "...")]` derive attribute on each frame.

---

## 9. Summary Table

| Shorthand | Full name | Axes fixed by | Typical use |
|---|---|---|---|
| ICRS | Int'l Celestial Ref. System | Quasar positions | Star catalogues |
| GCRS | Geocentric CRS | treated as ICRS-aligned for directions | Aberration, parallax |
| MeanJ2000 | Equatorial Mean J2000 | IAU 2006 `bp06` J2000 bias | Internal working frame |
| EME2000 | Earth Mean Equator/Equinox J2000 | Same axes as MeanJ2000 | CCSDS / flight-dynamics exchange |
| FK4 B1950 | FK4 Equatorial B1950 | FK4 catalogue | Legacy catalogues |
| TEME | True Equator Mean Equinox | SGP4 convention | TLE/satellite orbits |
| Galactic | Galactic (IAU 1958) | Galactic plane + NGP | Milky Way structure |
| EclMeanJ2000 | Ecliptic Mean J2000 | IAU 2006 `ecm06` at J2000 | Planetary ephemerides |
| CIRS | Celestial Intermediate RS | CIP + CIO | Intermediate Earth-orientation diagnostics |
| TIRS | Terrestrial Intermediate RS | CIP + ERA | Intermediate Earth-orientation diagnostics |
| ECEF | Earth-Centred Earth-Fixed | Provider-aligned with ITRF; WGS84 geodetic tag | WGS84 site input / Earth-fixed vectors |
| ITRF | Int'l Terrestrial Ref. Frame | Full EOP chain | Geodetic / VLBI |
| MercuryFixed | Mercury body-fixed | IAU 2015 rotation | Mercury surface |
| VenusFixed | Venus body-fixed | IAU 2015 rotation | Venus surface |
| MarsFixed | Mars body-fixed | IAU 2015 rotation | Mars surface |
| MoonIauMeanEarth | Moon analytical Mean-Earth | NAIF IAU_MOON (2009), periodic | Lunar surface |
| JupiterSysIII | Jupiter System III | IAU 2015 rotation | Jupiter magnetosphere |
| SaturnFixed | Saturn body-fixed | IAU 2015 rotation | Saturn surface |
| UranusFixed | Uranus body-fixed | IAU 2015 rotation | Uranus surface |
| NeptuneFixed | Neptune body-fixed | IAU 2015 rotation | Neptune surface |
| PlutoFixed | Pluto body-fixed | IAU 2015 rotation | Pluto surface |
| Horizontal | Local horizon | Observer's local vertical | Telescope pointing |

---

## 10. Current Non-Goals

Local orbital frames (`RTN`, `LVLH`, `VNC`) and covariance transport between
them and the inertial Cartesian frame are now first-class through the
orbit-relative API in `astro::dynamics::frames::LocalOrbitalFrame<M>` and
`astro::dynamics::covariance::StateCovariance<F>`.  Each instance is built
from a concrete `OrbitState`, so the state-dependent nature of the frame is
preserved at the type level.

Other operational frames (e.g. body-fixed frames for minor planets) remain
out of scope for the typed frame system. Generic surface observers for the
major rocky bodies are in scope as of v0.16 via `coordinates::surface`.

---

*Updated for the v0.18.0 frame and fallible center-shift contracts.*
