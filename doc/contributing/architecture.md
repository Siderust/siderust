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
