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
| `Geocentric` | `centers::Geocentric` | Earth's centre of mass |
| `Topocentric` | `centers::Topocentric` | Observer's location on Earth's surface |
| `Heliocentric` | `centers::Heliocentric` | Sun's centre of mass |
| `Bodycentric` | `centers::Bodycentric` | Centre of a named Solar-System body |

`Topocentric` is *parameterised*: it carries a `Geodetic<ECEF>` so horizontal
coordinates know their own observation point without external context.

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

### 3.4 Ecliptic Frames

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

## 5. The Earth-Rotation Chain

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

## 6. Geodetic Coordinates and Ellipsoid Conversion

### 6.1 `Geodetic<F, U>` (siderust) / `ellipsoidal::Position<C, F, U>` (affn)

```rust
pub type Geodetic<F, U = Meter> = affn::ellipsoidal::Position<Geocentric, F, U>;
```

`Geodetic<F>` is the first-class type for Earth-surface positions.  The ellipsoid
is encoded in the frame `F` via `HasEllipsoid`:

- `Geodetic<ECEF>` — WGS84
- `Geodetic<ITRF>` — GRS80

### 6.2 Converting to/from Cartesian ECEF

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

## 7. Summary Table

| Shorthand | Full name | Axes fixed by | Typical use |
|---|---|---|---|
| ICRS | Int'l Celestial Ref. System | Quasar positions | Star catalogues |
| GCRS | Geocentric CRS | ≈ ICRS (< 1 mas offset) | Aberration, parallax |
| MeanJ2000 | Equatorial Mean J2000 | IAU 1976 precession | Internal working frame |
| EclMeanJ2000 | Ecliptic Mean J2000 | Ecliptic plane of J2000 | Planetary ephemerides |
| ECEF | Earth-Centred Earth-Fixed | ERA only (no polar motion) | First-order site positions |
| ITRF | Int'l Terrestrial Ref. Frame | Full EOP chain | Geodetic / VLBI |
| Horizontal | Local horizon | Observer's local vertical | Telescope pointing |

---

*Last updated by the P0–P2 coordinate model refactor.*
