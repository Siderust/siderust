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

For planets (Mercury–Pluto), barycentric positions are computed via Keplerian
orbit propagation (from `solar_system.rs`) plus Sun's barycentric offset.
For the Moon (Selenocentric), the position comes from the ephemeris
`moon_geocentric()` combined with Earth's barycentric position.

---

## 2. Earth-Fixed Frames

### 2.1 ECEF — Earth-Centred Earth-Fixed (mathematical placeholder)

```rust
// affn::frames::ECEF
pub struct ECEF;
```

`ECEF` is a **generic** Earth-fixed label.  It rotates with the Earth using
ERA / GMST but intentionally **does not** apply the IERS polar-motion matrix
**W**(xₚ, yₚ, s′).

**When to use `ECEF`:**
- First-order geodetic → topocentric conversions where ≤ 10 m accuracy suffices.
- Internal bookkeeping when a labelled Earth-fixed frame is needed before a
  full EOP chain is available.

**Accuracy note:**  Omitting polar motion introduces an error of roughly
±10 m (up to ~30 m at solar maximum) in geocentric Cartesian coordinates.
For observatory positioning at the metre level or better, store observatory
data in `Geodetic<ECEF>` and call `.to_cartesian::<Meter>()`, which uses the
WGS84 ellipsoid encoded in the `ECEF` frame.

### 2.2 ITRF — International Terrestrial Reference Frame (EOP-realised)

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

---

## 3. Celestial (Inertial) Frames

### 3.1 ICRS / ICRF

```rust
pub struct ICRS;  // Spherical  (α, δ)
pub struct ICRF;  // Cartesian unit vector
```

The **International Celestial Reference System** is the primary inertial frame.
Axes are fixed by extragalactic radio sources (quasars).  It is realised by
the ICRF catalogue.

**Note:** `ICRS` and `ICRF` differ only in coordinate style (spherical vs.
Cartesian direction); they share the same axes.

### 3.2 GCRS — Geocentric Celestial Reference System

```rust
pub struct GCRS;
```

**Approximation in siderust:** `GCRS` currently uses `frames::ICRS` axes
under the hood.  The true GCRS is rotated from ICRS by < 1 mas (the
frame-bias matrix **B**), which is below the accuracy floor of most
calculations.  If sub-mas accuracy is required, apply the frame-bias
explicitly.

### 3.3 EquatorialMeanJ2000

```rust
pub struct EquatorialMeanJ2000;
```

Mean equatorial frame referred to the standard epoch J2000.0
(JD 2 451 545.0 TT).  Axes are defined by the mean ecliptic and equinox of
J2000, using the IAU 1976/1980 precession-nutation constants.

This is the **working frame** for most internal transforms in siderust
because precession-only rotations are cheap and widely tabulated.

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

### 3.5 TEME — True Equator, Mean Equinox

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
| `EclipticMeanJ2000` | Mean ecliptic and equinox of J2000 |
| `EclipticMeanOfDate` | Mean ecliptic and equinox of observation epoch |
| `EclipticTrueOfDate` | True ecliptic (adds nutation in obliquity and longitude) |

---

## 4. Horizontal Frame

```rust
pub struct Horizontal;
```

Local horizontal frame centred on the observer.  Spherical coordinates are
**(altitude, azimuth)** — altitude measured from the horizon (+90° = zenith),
azimuth measured from North through East.

Every `Position<Topocentric, Horizontal, U>` implicitly carries a
`Geodetic<ECEF>` via the `Topocentric` center parameter.

---

## 5. Planetary Body-Fixed Frames

These frames rotate with the respective body, defined by IAU 2015 rotation
parameters (α₀, δ₀, W).  The transformation from body-fixed to ICRS is:

```
R = Rz(−(α₀ + 90°)) · Rx(−(90° − δ₀)) · Rz(−W)
```

where α₀ and δ₀ define the body's north pole direction in ICRS, and W is
the prime meridian angle (degrees).

The **frame marker types** (zero-sized structs) are defined in `affn::frames`
behind `feature = "astro"` so that `affn`'s derive macros can generate
inherent `Direction`/`Position` constructors and getters (`lat()`, `lon()`,
`radius()`).  They are re-exported from `siderust::coordinates::frames`
and `siderust::coordinates::frames::planetary` for convenience.

The **IAU rotation parameters** and all `FrameRotationProvider` implementations
remain in `siderust::coordinates::frames::planetary` and
`siderust::coordinates::transform::providers`.

| Frame | Struct | Body | Spherical coords |
|---|---|---|---|
| Mercury body-fixed | `MercuryFixed` | Mercury | (lat, lon, radius) |
| Venus body-fixed | `VenusFixed` | Venus | (lat, lon, radius) |
| Mars body-fixed | `MarsFixed` | Mars | (lat, lon, radius) |
| Moon PA | `MoonPrincipalAxes` | Moon | (lat, lon, radius) |
| Jupiter System III | `JupiterSystemIII` | Jupiter | (lat, lon, radius) |
| Saturn body-fixed | `SaturnFixed` | Saturn | (lat, lon, radius) |
| Uranus body-fixed | `UranusFixed` | Uranus | (lat, lon, radius) |
| Neptune body-fixed | `NeptuneFixed` | Neptune | (lat, lon, radius) |
| Pluto body-fixed | `PlutoFixed` | Pluto | (lat, lon, radius) |

### IAU Rotation Parameters

Each body's rotation is defined by three functions of time:

- **α₀(T)** — right ascension of the north pole (degrees), T in Julian centuries
- **δ₀(T)** — declination of the north pole (degrees), T in Julian centuries
- **W(d)** — prime meridian angle (degrees), d in days from J2000.0

These are stored as `IauRotationParams` constants (e.g., `MARS_ROTATION`)
in `coordinates::frames::planetary`.

---

## 6. The Earth-Rotation Chain

The canonical ITRS → EquatorialMeanJ2000 rotation is implemented in
`crate::astro::earth_rotation_provider::itrs_to_equatorial_mean_j2000_rotation`.

```
ITRS  ──W⁻¹──▶  TIRS  ──ERA──▶  CIRS  ──Q──▶  GCRS/ICRS  ──P──▶  EquatorialMeanJ2000
```

| Step | Matrix | Source |
|---|---|---|
| **W** | Polar motion (xₚ, yₚ, s′) | EOP → `polar_motion_matrix_from_eop` |
| **ERA** | Earth Rotation Angle | UT1 → `earth_rotation_angle` |
| **Q** | CIO/CIP (X, Y, s) | nutation IAU 2000B ± EOP dX,dY |
| **P** | Precession ICRS → MeanJ2000 | `FrameRotationProvider<ICRS, EquatorialMeanJ2000>` |

**Time-scale contract:**  The `jd` argument to `itrs_to_equatorial_mean_j2000_rotation`
is a **Terrestrial Time (TT)** Julian Date.

---

## 7. Frame-Rotation Hub

All frame rotations route through the **ICRS** hub frame:

```
F1 → ICRS → F2
```

The following constant rotation matrices are pre-computed:

| Rotation | Matrix constant | Reference |
|---|---|---|
| FK4 B1950 → ICRS | `FK4_TO_ICRS` | Standish (1982) + frame bias |
| Galactic → ICRS | `GALACTIC_TO_ICRS` | Murray (1989), NGP coordinates |
| ICRS → EquatorialMeanJ2000 | `FRAME_BIAS_ICRS_TO_J2000` | IAU 2006 frame bias |

Time-dependent rotations:

| Pair | Method |
|---|---|
| TEME ↔ TOD | Rz(equation of equinoxes) |
| EquatorialMean ↔ EquatorialTrueOfDate | IAU 2006 precession + IAU 2000B nutation |
| Body-fixed ↔ ICRS | IAU 2015 rotation parameters |

---

## 8. Geodetic Coordinates and Ellipsoid Conversion

### 8.1 `Geodetic<F, U>` (siderust) / `ellipsoidal::Position<C, F, U>` (affn)

```rust
pub type Geodetic<F, U = Meter> = affn::ellipsoidal::Position<Geocentric, F, U>;
```

`Geodetic<F>` is the first-class type for Earth-surface positions.  The ellipsoid
is encoded in the frame `F` via `HasEllipsoid`:

- `Geodetic<ECEF>` — WGS84
- `Geodetic<ITRF>` — GRS80

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
| GCRS | Geocentric CRS | ≈ ICRS (< 1 mas offset) | Aberration, parallax |
| MeanJ2000 | Equatorial Mean J2000 | IAU 1976 precession | Internal working frame |
| FK4 B1950 | FK4 Equatorial B1950 | FK4 catalogue | Legacy catalogues |
| TEME | True Equator Mean Equinox | SGP4 convention | TLE/satellite orbits |
| Galactic | Galactic (IAU 1958) | Galactic plane + NGP | Milky Way structure |
| EclMeanJ2000 | Ecliptic Mean J2000 | Ecliptic plane of J2000 | Planetary ephemerides |
| ECEF | Earth-Centred Earth-Fixed | ERA only (no polar motion) | First-order site positions |
| ITRF | Int'l Terrestrial Ref. Frame | Full EOP chain | Geodetic / VLBI |
| MercuryFixed | Mercury body-fixed | IAU 2015 rotation | Mercury surface |
| VenusFixed | Venus body-fixed | IAU 2015 rotation | Venus surface |
| MarsFixed | Mars body-fixed | IAU 2015 rotation | Mars surface |
| MoonPA | Moon Principal Axes | IAU 2015 rotation | Lunar surface |
| JupiterSysIII | Jupiter System III | IAU 2015 rotation | Jupiter magnetosphere |
| SaturnFixed | Saturn body-fixed | IAU 2015 rotation | Saturn surface |
| UranusFixed | Uranus body-fixed | IAU 2015 rotation | Uranus surface |
| NeptuneFixed | Neptune body-fixed | IAU 2015 rotation | Neptune surface |
| PlutoFixed | Pluto body-fixed | IAU 2015 rotation | Pluto surface |
| Horizontal | Local horizon | Observer's local vertical | Telescope pointing |

---

*Last updated by the P0–P2 coordinate model refactor.*
