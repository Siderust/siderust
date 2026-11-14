# Observe from a lunar site using SPICE kernels

Enable `spice` and use `LunarObserver` when a site is defined in a
kernel-based lunar frame. One borrowed `SpiceContext` supplies both
Solar System states and lunar orientation. The observer uses the kernel
states directly and computes geometric Sun and Earth center directions,
including finite-distance parallax.

## Kernel set

The [runnable example](../../examples/24_spice_lunar.rs) loads these files
in order. Obtain them from the [NAIF generic kernel archive](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/)
and place them in one directory. The archived FK is under
`fk/satellites/a_old_versions/`. Compare the files with the
[reference hashes](../../tests/fixtures/lunar_de440/provenance.json)
when reproducing the regression.

| File | Role |
| --- | --- |
| `de440s.bsp` | SPK states for Sun, Earth–Moon barycenter, Earth, and Moon in NAIF J2000. |
| `moon_pa_de440_200625.bpc` | Type 2 binary PCK for the lunar principal-axes orientation, class ID 31008, based on J2000. |
| `moon_de440_220930.tf` | Frame kernel defining the PA and ME frames and their relationship. |
| `pck00011.tpc` | Text PCK constants and analytical orientation metadata. The lunar observer uses binary orientation. |
| `naif0012.tls` | Leap-second kernel used by external SPICE consumers and reference generation. Tempoch supplies time conversion in the observation path. |

```sh
cargo run --features spice --example 24_spice_lunar -- /path/to/kernels
```

Use a compatible SPK, binary PCK, and frame kernel. The implementation
checks the frame chain and selected segments, but descriptors do not
authenticate a DE release. Retain file hashes and load order with results.
Result metadata lists loaded candidates; it does not prove which file
supplied each segment.

## Site and epoch

In this FK, `MOON_PA` aliases `MOON_PA_DE440` (31008).
`MOON_ME` aliases `MOON_ME_DE440_ME421` (31009), a fixed offset
from the DE440 principal axes aligned to the DE421 Mean-Earth realization.
Site longitude and latitude are interpreted in the selected runtime frame.
The analytical `MoonIauMeanEarth` frame is a separate realization.
See the [DE440 frame kernel](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/a_old_versions/moon_de440_220930.tf)
and the [NAIF frame conventions](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html).
For direct comparisons, `rotation_naif(from, to, et)` maps column vectors
from source to destination axes. Its inverse is the transpose. A fixed
TK definition maps the defined frame to its relative frame, while a binary
PCK rotation maps base axes to body-fixed axes.

The example uses `LunarObserver::radial` with east-positive longitude,
planetocentric latitude, and an actual distance from the Moon's center
in metres. Local up is radial. To use height above a datum instead, pass
a `SurfaceSite<Moon>` to `LunarObserver::new`. Its default datum is the
1737.4 km sphere. For a custom ellipsoid, local up follows the surface normal.
At a pole, site longitude defines the local meridian. The horizontal API
returns zero azimuth at the zenith/nadir singularity.

`observe` accepts `Time<TDB>`, evaluated as seconds past J2000 (SPICE
ephemeris time, ET). Convert TT explicitly with Tempoch. For UTC, use the
appropriate civil time-data context before converting. Loading an LSK
does not swap Tempoch's leap-second data. Pass identical ET to both
implementations when isolating geometry from time conversion.

The output includes altitude, azimuth north through east, and slant
range. Visibility means target center altitude strictly above zero. The
geometry corresponds to SPICE `abcorr="NONE"`. Light-time, stellar
aberration, relativistic deflection, refraction, disc extent, and terrain
are outside this calculation.

## Diagnose a failed observation

Construction requires a Moon-centered FK chain terminating at a
J2000-based binary PCK. Evaluation requires covering, supported SPK and
PCK data at the same ET. BPC coverage is inclusive, and both state and
orientation must cover the epoch. For the archived BPC in the example, class 31008 spans
ET −14200747200 through 13447252800 seconds. The SPK chain can restrict
the usable interval further. Later eligible files and segments have
priority within coverage. A selected unsupported segment produces an error
with no fallback to older overlapping data or analytical models. Check the
[loader contract](../reference/guarantees.md#runtime-kernel-selection) for
the supported subset. The [CSPICE regression](../../tests/fixtures/lunar_de440/README.md)
records the kernel excerpts, reference generation, and measured agreement.
