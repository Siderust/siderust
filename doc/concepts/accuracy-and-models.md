# Model accuracy and reproducibility

Analytical ephemerides evaluate series of orbital terms. Siderust bundles
VSOP87 for the Sun and planets and ELP2000-82B for the Moon. Compact
ephemerides use truncated series selected for a bounded interval. JPL
Development Ephemeris (DE) kernels supply piecewise polynomial states from
a numerical ephemeris. The accuracy of these approaches depends on the
body, epoch, data, and quantity being compared. The ability to evaluate an
analytical series at an epoch does not establish its physical validity there.
See the [DE440 and DE441 model description](https://ssd.jpl.nasa.gov/doc/de440_de441.html).

## What the frame label establishes

Public ephemeris coordinates use the IAU 2006 `EclipticMeanJ2000` convention.
Full and compact VSOP use the published VSOP87 native-to-FK5 tie. ELP first
uses its Laskar step from its moving native axes to the native J2000 triad,
then applies the separate ELP native-to-FK5 tie. Both paths use the fixed J2000
orientation from ERFA `fk5hip` and then the canonical ICRS-to-ecliptic
rotation. The tie is empirical. It does not remove catalog zonal errors
or calibrate old analytical theories to DE440. The rotations are static,
so the VSOP position and velocity share the same mapping without a frame-rate
term. User-supplied Keplerian elements are already defined in their declared
frame and should keep that frame convention. The legacy full `Moon::vsop87a` accessor evaluates the Earth–Moon
barycenter coefficients. Use the Moon ephemeris for a geocentric lunar state.
See the [VSOP87 conventions](https://ftp.imcce.fr/pub/ephem/planets/vsop87/vsop87.doc),
[ELP82B archive catalog and notice](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/VI/79?format=html&tex=true)
and [ERFA FK5 orientation](https://github.com/liberfa/erfa/blob/v2.0.1/src/fk5hip.c).

For producer-documented ICRF-realized DE kernels, NAIF frame ID 1
(`J2000`) denotes orientation consistent with ICRF. The runtime ephemeris
applies the canonical ecliptic mapping after computing the physical chain.
The literal ID cannot establish a historical kernel's realization.
See the [NAIF frame documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html).

## Compare like observables

A geometric position is evaluated at the chosen epoch. An apparent
observation may also use light-time (retarded emission epoch), stellar
aberration (observer motion), gravitational deflection, and atmospheric
refraction. Observer parallax follows from subtracting the site position
from a finite target. A local equatorial or horizontal frame alone does
not specify which corrections an API applies. Check the particular
observation function and reference output before comparing.

Record the Siderust revision and Cargo features, model/backend, input
epoch and scale, units, target and origin, frame realization, site datum,
correction settings, and identities of time-data and ephemeris files.
For state comparisons, retain both the velocity convention and the
position convention. For near-zenith cases use an angular separation
between directions alongside any circular azimuth difference.

Independent regressions cover
[analytical frame normalization](../../tests/ephemeris_frame_consistency.rs),
[Earth orientation](../../tests/test_high_precision_earth_rotation_regression.rs),
[body-surface geometry](../../tests/body_surface_horizontal.rs),
[compact pointing](../../tools/validation/README.md), and
[same-kernel lunar geometry](../../tests/fixtures/lunar_de440/README.md).
The [analytical lunar axis test](../../tests/body_surface_horizontal.rs)
compares all three axes with CSPICE and `pck00011.tpc` at five TT epochs from
1900 through 2050, with a direction-cosine tolerance of 10⁻¹⁰. Regeneration
uses the [surface reference script](../../tools/cspice_refs/cspice_surface_references.py).

These comparisons establish agreement only for their tested data, epochs,
operations, and reference models. Separate approximation error, data
uncertainty, and numerical residuals when setting an application error budget.
