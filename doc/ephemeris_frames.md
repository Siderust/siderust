# Ephemeris reference frames

`EclipticMeanJ2000` in Siderust means the **IAU 2006 mean ecliptic and
equinox at J2000.0**. For column vectors its ICRS mapping is
`C = Rx(-84381.406 arcsec) * rb`, where `rb` is SOFA `bp06` at J2000.
This is SOFA/ERFA `ecm06(2451545.0, 0.0)`. `EquatorialMeanJ2000` and
`EME2000` use that same `rb` mean equator/equinox. They do not implement
the empirical FK5 *catalogue* orientation. ICRF is treated as an ICRS
realization with identical axes, to the accuracy of the source ephemeris.

`Position<Center, Frame, Unit>` describes physical coordinates. Changing
the calculation method does not change this meaning. Different theories
estimate different trajectories, with different errors, in these axes.

## Native analytical realizations

The numerical VSOP kernels remain native. Public planetary `vsop87a`,
`vsop87e`, velocity, and combined-state methods normalize **before**
constructing a typed coordinate. `Vsop87Ephemeris` delegates to them; no
second rotation occurs in that backend or in center transformations. The legacy
`Moon::vsop87a` method evaluates the VSOP Earth–Moon barycenter coefficients;
it is not used for the Moon's geocentric backend state.

For VSOP87A/E use the matrix `V` from the IMCCE `vsop87.doc` **REFERENCE
SYSTEM** section (rows shown, multiplying column vectors):

```text
 1               +0.000000440360  -0.000000190919
-0.000000479966   +0.917482137087  -0.397776982902
 0               +0.397776982902  +0.917482137087
```

It maps the dynamical inertial ecliptic/equinox J2000 of VSOP87A/E to FK5
J2000 equatorial coordinates. The original paper distinguishes this
inertial dynamical equinox from the rotational dynamical equinox and
from DE200. Rotating only by a standard obliquity omits the equinox tie.

For ELP2000-82B the series longitude is measured from a **departure point
on the moving inertial mean ecliptic**, not the equinox of date. The
existing Laskar matrix maps that triad to ELP's inertial mean ecliptic
and equinox J2000. Preserve that step. The ELP82B notice §8, p.12 gives
its **separate** native-to-FK5 matrix `L`:

```text
 1               +0.000000437913  -0.000000189859
-0.000000477299   +0.917482137607  -0.397776981701
 0               +0.397776981701  +0.917482137607
```

The notice derives this from its fitted obliquity 84381.40883 arcsec
(±0.00006 arcsec) and equinox arc -0.09845 arcsec (±0.00016 arcsec).
The ELP and VSOP native triads are therefore **not assumed identical**.
`Moon::get_geo_position` converts TT to TDB as required by the ELP notice,
keeps the Laskar rotation, then applies the ELP normalization once.

For both paths the final mapping is `C * H * V` or `C * H * L`. `H` is
the **fixed J2000 orientation** from SOFA/ERFA `fk5hip` (Mignard & Froeschlé
2000; orientation components -19.9, -9.1, +22.9 mas). The Hipparcos system
is used as the optical realization of ICRS. This is a conventional,
empirical FK5-to-ICRS tie, not a new fit to a modern planetary ephemeris.
Catalogue zonal errors and the uncertainty of the historical dynamical/FK5
tie are not eliminated. Do not interpret the result as an exact recovered
DE200-to-DE440 alignment or milliarcsecond-accurate analytical ephemeris.

The catalogue **spin** from `fk5hip` is deliberately excluded: these
are fixed J2000 axes, not FK5 catalogue positions/proper motions propagated
to the evaluation date. Thus `v_canonical = R * v_native`, without an
`Rdot * position` term. VSOP velocities retain their AU/TDB-day convention.
The two fixed matrices are composed once with `spin::Once`. Their published
12-decimal rounding gives orthogonality errors below 1e-12; tests account
for that precision rather than claiming machine-exact inverse matrices.

## JPL and SPK

DE440/441 are ICRF realizations. NAIF labels their axes `J2000`, frame
ID 1. This does **not** mean classical FK5/EME2000. The runtime path
already correctly computes `C * position_icrf` and `C * velocity_icrf`
after target/center chain arithmetic. No additional FK5 bias is applied.
The public frame-bias extension traits now share the pipeline's `bp06`
matrix instead of independently constructing a nearly equal matrix.

The ICRF-only SPK parser and legacy body-chain parser reject unsupported
frame IDs. `RuntimeEphemeris::try_from_indexed_segments` and
`SpkKernelSet::try_from_indexed_segments` validate metadata too; their
existing infallible constructors panic on invalid IDs. A non-J2000 segment
cannot silently disappear and expose an older overlapping segment.

**Frame ID 1 is necessary, not sufficient, to establish ICRF provenance.**
The loader contract requires the producer's documented ICRF realization,
as for DE440/441. Historical DE200/FK5 kernels are not supported by these
ICRF APIs. Neither a filename nor an ID proves physical realization.
Legacy `from_segments(BspSegments)` has no frame metadata; its caller
explicitly supplies ICRF data. General frame-aware SPICE operations remain
separate from these ICRF-only ephemeris loaders.

## Consumers and compatibility

- Sun/Earth and selenocentric translations, observer states, aberration,
  physical orientation, and spherical conversions consume normalized
  positions/velocities through their existing paths.
- Mercury–Neptune center shifts now use
  `Ephemeris::try_planet_barycentric`, including reverse and composed shifts.
  The default reports unsupported data. A custom static backend that uses
  planetary-center transformations must implement this capability with its
  own data (it may compose its `PlanetState<Body>` and Sun state). It can no
  longer inherit a silent VSOP87E substitution. `PlanetState` stays generic.
- `RuntimeEphemeris` continues to provide Sun/Earth/Moon through
  `DynEphemeris`; this change does not add runtime planetary-center support.
- Keplerian/conic propagation operates on supplied canonical-frame elements;
  it must not apply a backend correction to user-supplied elements. The
  built-in coarse element catalogue and the abbreviated Meeus Pluto series
  do not carry a precision realization tie. Their subarcsecond orientation
  uncertainty remains a limitation; no VSOP or ELP correction was guessed
  for them. The default backend's Pluto capability retains the catalogue Kepler
  approximation plus its Sun, not a high-accuracy DE Pluto state. A custom
  backend must explicitly implement `pluto_barycentric`; its default panics
  rather than silently substituting the catalogue model.
- `MoonIauMeanEarth` is now a **distinct marker** for analytical `IAU_MOON`.
  `MoonPrincipalAxes` remains the `affn` principal-axes marker and has no
  analytical provider. Migrate analytical lunar coordinates and transforms
  to `MoonIauMeanEarth`. DE-specific `MOON_PA`/`MOON_ME` needs binary PCK and
  frame kernels; it cannot be obtained by relabeling analytical coordinates.
- FFI/Python callers of Rust analytical ephemerides inherit the normalized
  numbers without changing coordinate layouts. Historical numerical snapshots
  must be updated. No ephemeris-specific frame generics or dependency changes
  are required; `no_std + alloc` is preserved.

## Independent validation

`tests/data/ephemeris_frames.json` is generated independently by
`scripts/validate_ephemeris_frames.py`, using IMCCE check data, ELP Table H,
ERFA 2.0.1.5 and CSPICE N0067 with the official `de440s.bsp`.

- VSOP A/E Mercury, Earth, Mars and Neptune check states: J2000, 1900 and
  1100; compare canonical position and analytical velocity. The printed
  IMCCE precision is 1e-10 AU or AU/day; test tolerance is 1e-10.
- Direct scalar evaluation of the original IMCCE coefficients at ±4000 years
  (Mercury/Earth/Mars) and ±6000 years (Neptune), for both A/E, checks the
  nominal-span endpoints to 1e-9 AU (scalar versus SIMD/FMA roundoff).
- All five ELP Table H full-series examples (1828–2047), checked to 0.1 m,
  including signs independently verified on the rendered notice. Table H
  has 1e-5 km quantization.
- Independent ERFA basis vectors test matrix direction and handedness;
  inverse tests allow 2e-12 for the rounded source matrices. Public ICRS
  round trips and analytic-velocity finite differences cover ±1000 years.
- Runtime DE440s Sun, Earth barycentric/heliocentric, Moon geocentric, and
  Earth velocity at 1900/J2000/2026/2100: CSPICE `spkez(...,'J2000','NONE',...)`
  plus ERFA `ecm06`, 1 mm position and 1e-7 km/day velocity tolerances.
  The kernel-dependent test is explicitly ignored unless requested.
- Backend-dispatch tests deliberately supply unrelated Mars coordinates;
  metadata tests exercise parsed segments and byte/legacy loaders.
  Compile-fail documentation prevents analytical/principal-axes interchange.

See [numerical comparison](ephemeris_frame_comparison.csv) for the **pre-fix
native output** comparison against DE440s at six epochs over 1850–2149,
covering Mercury, Earth, Mars, Neptune and Moon. Positions are km, angular
separations arcseconds, planetary velocities m/s. Lunar velocities are
omitted because ELP's public API supplies positions only. Each row separates
un-normalized residual, orientation-only displacement, and residual after
orientation normalization. This is not a claimed error bound for either model.

At J2000, Earth changes by 78.014 km (0.109390 arcsec), with its DE440s
position residual changing from 49.737 to 29.267 km. Its velocity correction
is 0.016275 m/s, residual 0.009585 → 0.006760 m/s. Moon's orientation-only
change is 215.860 m (0.110633 arcsec), residual 233.688 → 68.384 m. The
separate TT→TDB correction is about 0.08 m there (about 1.59 m in 2026).
Some residuals increase: e.g. Mars in 2100, 31.178 → 151.830 km. Orientation
must follow references regardless of which direction a model residual moves.

The CSV's TT/TDB column estimates the displacement from confusing the two
time scales using ERFA `dtdb` and the CSPICE speed. Planetary evaluation
already used TDB. SPK evaluation retains its existing NAIF DELTET TT→ET
conversion, matched in the reference script. Units use the exact SI AU
149597870.700 km and 86400 s/day. These conversions add no axis rotation;
this work preserves the existing AU interpretation of the coefficients and does
not refit their historical length scale. No rotation is fitted to DE440s.

IMCCE describes VSOP's nominal 1 arcsec spans as ±4000 years for inner
planets, ±2000 for Jupiter/Saturn, ±6000 for Uranus/Neptune; these concern
agreement with its historical fit, not a DE440 accuracy guarantee. ELP's
historical DE200 fit also leaves secular/model residuals, which grow away
from the fit interval. DE440s covers only 1849–2150: this comparison does
not validate modern-kernel residuals over the full analytical validity span.

Reproduce with an external kernel (not downloaded by tests):

```sh
uv run --with pyerfa==2.0.1.5 --with spiceypy==8.3.0 python \
  scripts/validate_ephemeris_frames.py --kernel /path/to/de440s.bsp \
  --fixtures /tmp/ephemeris_frames.json
SIDERUST_BSP_PATH=/path/to/de440s.bsp cargo test \
  --test ephemeris_frame_consistency runtime_de440s_matches_cspice_erfa_reference -- --ignored
```

For before/after output, run the `frame_consistency` example against the
parent revision and this revision in separate worktrees, saving CSVs.
Copy the example into the parent worktree first (it did not exist there).
Pass `--native-csv before.csv --canonical-csv after.csv --kernel ...` to the
script. It independently checks planetary normalization to 1e-5 km and
1e-8 m/s; Moon position allows its separate ≤2 ms time-scale correction.

## References

- [Bretagnon & Francou (1988), A&A 202, 309–315](https://articles.adsabs.harvard.edu/pdf/1988A%26A...202..309B), especially p.309.
- [IMCCE VSOP87 conventions and precision](https://ftp.imcce.fr/pub/ephem/planets/vsop87/vsop87.doc) and [check states](https://ftp.imcce.fr/pub/ephem/planets/vsop87/vsop87.chk).
- [Chapront-Touzé & Chapront, ELP2000-82B notice](https://ftp.imcce.fr/pub/ephem/moon/elp82b/elp82b.ps), §§8–9, and [reference code](https://ftp.imcce.fr/pub/ephem/moon/elp82b/elp82b_2).
- [ERFA `fk5hip`](https://github.com/liberfa/erfa/blob/v2.0.1/src/fk5hip.c) and [`ecm06`](https://github.com/liberfa/erfa/blob/v2.0.1/src/ecm06.c).
- [NAIF frame conventions, ICRF vs J2000](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#ICRF%20vs%20J2000).
- [DE440/441, Park et al. (2021)](https://doi.org/10.3847/1538-3881/abd414).
- [NAIF lunar frame kernel documentation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf), distinguishing IAU_MOON, PA and ME.

## Performance checks

Criterion on the same host, parent `087d7cf` versus this change (20 samples,
1 s warmup, 2 s measurement): Earth VSOP87A 25.61–26.06 µs versus
20.82–20.93 µs; VSOP87E 39.75–39.99 µs versus 38.00–39.05 µs; ELP offset
372.84–375.26 µs versus 370.48–373.66 µs. Concurrent test load makes the
apparent VSOP improvements unsuitable for a speedup claim. The ELP comparison
detected no significant change. These measurements show no observed regression,
not a precise bound on the small extra fixed matrix application.
