# Crate boundaries and transform invariants

Siderust owns astronomical bodies, provider integration, frame and center
transforms, observing and event workflows, and satellite and orbit
determination applications. Its dependencies provide reusable primitives
and datasets.

| Crate | Responsibility |
| --- | --- |
| `qtty` | Typed physical quantities, conversions, and arithmetic. |
| `affn` | Affine geometry, coordinate and frame types, geometric operators, and conic geometry. |
| `tempoch` | Time instants, scales, formats, intervals, and the interpretation of active time data. |
| `siderust-archive` | Scientific datasets, manifests, provenance, coefficient accessors, and optional acquisition. |
| `keplerian` | Two-body propagation, Kepler equation and anomaly solvers, and Lambert transfers. |
| `cheby` | Reusable Chebyshev approximation, evaluation, and calculus. |
| `principia` | Numerical dynamics, force-model primitives, integrators, and covariance machinery. |
| `optica` | Typed optics, spectra, and radiative-transfer primitives. |

## Composition

Keep astronomical frame rotations, origin translations, and observation
corrections separate. Rotation providers usually compose through ICRS and
center shifts through the Solar System barycenter. Add direct primitives
in the appropriate [provider module](../../src/coordinates/transform/providers)
and compose inverses and indirect paths from them. Check planetary ZXZ
composition against independent basis vectors because the order of
coordinate rotations matters. Rotate ephemeris
displacements into the requested frame, using backend-specific
normalization exactly once at the ephemeris boundary. Preserve coverage
and capability errors through reverse and composed shifts.

A surface observation uses the body's center, orientation, and reference
surface. The body marker carries those relationships through
`BodySurfaceModel`. Add a default surface only when its datum, latitude
convention, and longitude sense are specified. The Earth-specific
sidereal paths retain their own epoch and EOP requirements.

A `CoordinateWithPM` stores a coordinate sample, epoch, and optional
proper motion. `Trackable` is the abstraction consumers use to request
coordinates at an epoch. A trackable fixed direction and a propagating
body use different physical models. Altitude and azimuth provider traits
supply observation-specific behavior independently.

## Binding boundaries

The [C ABI](../../siderust-ffi/README.md) owns astronomy handles, coordinate
tags, and observation semantics. Reuse quantity carriers from `qtty-ffi`
and time and interval carriers from `tempoch-ffi`. Preserve affine
distinctions, frames, centers, and units when crossing the ABI.
Keep handle ownership and destruction explicit. Check adapter exposure
in the relevant downstream repository before documenting it as supported.

## Public astronomy APIs

| Kind | Examples | When to expose it |
| --- | --- | --- |
| Composable primitive | `TransformFrame`, `PositionAstroExt`, `TransformCenter`, `Ephemeris`, `project_position_with` | The operation has a single astronomical meaning, with independently useful inputs and output. |
| Scientific algorithm | Lunar phase geometry, rise and set searches, light-time iteration | The algorithm computes a domain result or enforces an order of operations needed for physical correctness. |
| Convenience wrapper | `Sun::get_horizontal`, `Moon::get_horizontal`, default-context transforms | The entry point represents a common workflow and delegates to shared primitives. Document the model, time scale, origin, and any physical effects. |

Use frame, center, and unit types for compile-time distinctions. Use a trait
when callers need to substitute a provider or support several coordinate
types through the same contract. Use a function for a single scientific
operation with a concrete input and result. Avoid encoding choices of frame,
center, target, and model as a family of function names. Accept the existing
context when an operation needs time-dependent models. Evaluate each provider
at the appropriate epoch and pass the same model context through related
operations.

### Earth horizon coordinates

A geometric sightline from a finite-distance target starts with an ephemeris
position at the reception epoch. Translate it to the observer with
`TransformCenter::to_center_with`, rotate it to `EquatorialTrueOfDate` with
`PositionAstroExt::to_frame_with`, and call
`coordinates::transform::horizontal::project_position_with`. The position
carries the observer site. The frame provider uses UT1 for Earth rotation,
TT for the selected bias-precession-nutation model, and EOP for terrestrial
polar motion and observed celestial pole offsets. The offsets correct the
CIP in the CIO-based Earth orientation chain. The true equator and equinox
of date follow the selected nutation model.

For free directions, pass the site explicitly to `project_direction_with`
and `unproject_direction_with`. The position inverse
`unproject_position_with` carries the observer and range back to the
true-of-date equatorial frame. Use a subsequent typed frame rotation to
express the recovered sightline in other axes. To change a finite-distance
observer, translate the original target to the new observer first.

The explicit UT1/TT APIs use the equinox-based sidereal convention with
zero polar motion and celestial pole offsets. The mean-of-date transform
uses mean sidereal time. Body-fixed surface projection uses the body's
orientation model and local tangent basis. These paths share the relevant
local geometry and preserve their time and frame contracts.

Lunar altitude searches interpolate only geocentric Cartesian position.
They reuse the site's WGS84 position and local basis, and evaluate the
default terrestrial rotation at each TT query. A crossing solve refines
the interpolated altitude signal. Ephemeris and interpolation error must
be assessed separately from the root solver's time tolerance.

### Physical observation effects

Build the line of sight using the target state and the observer state.
Evaluate target motion when the model requires it. For light-time, solve for
the emission epoch using the target trajectory and the reception observer
state. Apply deflection and aberration in the order and frame defined by the
selected observation algorithm. Use the final horizon direction for any
atmospheric refraction model. Frame bias, precession, and nutation describe
axes. They apply to geometric sightlines as well as to physically corrected
sightlines.

Use the typed observation wrappers or the specific effect functions for
explicit corrections. If a policy's required data or effect is unsupported,
return an error. If the contract permits partial results, expose the applied
stages through `AppliedCorrections`. Existing target and event interfaces can
provide different levels of fidelity. In particular, the legacy stellar
horizontal policy selects how input angles are interpreted and uses a fixed
nutation model. A solver that needs full apparent observations should accept
the appropriate observation evaluator at its boundary.

When changing numerical paths, verify the applicable independent references
and compare implementations at matching epochs, sites, and models. For a
convenience API migration, publish a working composition with equivalent
semantics before deprecating the wrapper.
