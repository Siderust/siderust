# `atmosphere`

## Scientific scope

This module models line-of-sight effects introduced by the terrestrial
atmosphere in optical astronomy. The covered phenomena are the standard ones
used in extinction and sky-brightness work: geometric airmass, Rayleigh and
Mie optical depth, ozone transmission, Van Rhijn airglow geometry, and
Beer-Lambert attenuation. In other words, this is the layer that turns a
top-of-atmosphere radiometric quantity into a site-dependent ground-level
transmission factor.

The implemented formulas are literature-backed approximations rather than a
full radiative-transfer package. Their natural regime is ground-based optical
and near-infrared observing, with validity limited by the cited models and by
the assumptions documented in each submodule.

## Technical scope

What you will find here:

- `airmass.rs`: compile-time-selectable airmass formulas such as
  `PlaneParallel`, `Young1994`, `Rozenberg1966`, and
  `KrisciunasSchaefer1991`.
- `rayleigh.rs`: Rayleigh optical depth and phase-function helpers.
- `mie.rs`: aerosol / Mie optical-depth parametrization and site presets.
- `airglow.rs`: Van Rhijn geometry for emission-layer enhancement.
- `extinction.rs`: Beer-Lambert transmission `exp(-X * tau)`.
- `profile.rs`: `AtmosphereProfile`, the typed bundle of site atmospheric
  parameters used by higher-level callers.
- `ozone.rs`: wavelength-dependent ozone transmission tables when the
  `photometry` feature is enabled.
- `scattering.rs`: shared scattering traits and typed factors.

The public surface uses typed `qtty` values throughout: zenith angles,
wavelengths, pressures, heights, airmasses, optical depths, and
transmittances are never exposed as bare `f64`.

## Scattering phase conventions

`NormalizedPhaseTable` is the strict atmospheric option. Construct it with
positive ascending typed `Nanometers`, ascending typed `Radians` covering
exactly `[0, PI]`, and typed Optica `ScatteringFactor` values in angular-row-major
order (`values[iy * nx + ix]`). Optica checks finite monotone axes, dimensions,
and the two-sample minimum per axis. Siderust checks the physical domain,
finite nonnegative values, and each slice's unity normalization:

```text
2π ∫[0,π] P(θ) sin(θ) dθ = 1
```

The angle is between photon propagation directions. Values from external data
normalized to a 4π integral must explicitly be divided by 4π before construction.
Construction rejects unnormalized data; it never guesses or corrects a convention.

Normalization integrates the actual linear angular segments analytically.
A midpoint formulation, a small-interval Taylor series, and compensated summation
avoid cancellation in narrow segments and accumulation error in large tables.
The absolute integral tolerance is `1e-6`: it allows small discretization errors
in sufficiently resolved analytic samples and exceeds ordinary f64 roundoff.
It does not establish agreement with an unresolved analytic function, especially
near HG's singular limits. Linear wavelength interpolation preserves normalization
because it forms convex combinations of normalized angular slices, subject to
floating-point rounding. See the compiling example on `NormalizedPhaseTable`.

`try_phase` accepts exact boundaries and rejects nonfinite coordinates and all
out-of-domain queries. It uses Optica `PhaseTable::try_interp_at` with
`OutOfRange::Error`; no queries allocate. Neither infallible phase trait is
implemented. Optica's infallible `interp_at` and `PhaseModel::Tabulated.phase`
can clamp even with an Error policy, so they do not provide this contract.

`TabulatedPhaseFunction` remains the generic correction-table API: it accepts
incomplete or unnormalized data and clamps queries. Its Siderust
`PhaseFunction` trait and `ScatteringFactor` marker are distinct Rust types from
Optica's similarly named trait and marker. The strict adapter uses the Optica
marker consistently. Both markers are dimensionless numerical representations;
physical phase densities are per solid angle (`sr^-1`). Consumers needing typed
reciprocal-solid-angle densities must explicitly convert at that seam.

Use the existing `siderust::optica::phase` analytic models directly:
`RayleighPhaseFunction`, `HenyeyGreensteinPhaseFunction`, and
`DoubleHenyeyGreensteinPhaseFunction`. Their densities integrate to unity.
HG's `try_new` validates parameters; positive `g` means forward scattering,
negative `g` backward scattering, and zero isotropic scattering. Analytic models
assume wavelength independence and do not validate query coordinates. Very strong
peaks require adequate numerical resolution. The existing Siderust Rayleigh API
now delegates to Optica's kernel through an explicit marker conversion.

## References

- Bodhaine, B. A., Wood, N. B., Dutton, E. G., and Slusser, J. R. (1999).
  "On Rayleigh optical depth calculations." Journal of Atmospheric and Oceanic
  Technology 16, 1854-1861.
- Patat, F., Moehler, S., O'Brien, K., et al. (2011). "Optical atmospheric
  extinction over Cerro Paranal." Astronomy and Astrophysics 527, A91.
- Krisciunas, K., and Schaefer, B. E. (1991). "A model of the brightness of
  moonlight." Publications of the Astronomical Society of the Pacific 103,
  1033-1039.
- European Southern Observatory. SkyCalc and Advanced Cerro Paranal Sky Model.
  <https://www.eso.org/sci/software/pipelines/skytools/skycalc>

## Layered optical columns

`LayeredAtmosphere` adds finite, horizontally uniform, plane-parallel columns
without changing `AtmosphereProfile`. `AtmosphericOptics` groups existing Optica
coefficients into Rayleigh, aerosol and molecular absorption. Literal layers
are wavelength independent over an explicit caller-declared domain; `from_profile`
uses the existing Bodhaine/Patat wavelength laws over 360–570 nm and **explicit**
constituent layer fractions and aerosol SSA. Fractions must sum to one within
`8*n*f64::EPSILON`; only that rounding allowance is normalized. No measured
vertical distribution or aerosol scale height is inferred from a site preset.

Edges are absolute meters in one declared datum (MSL for site profiles), ordered
and contiguous by construction. Ground and observer are the first edge; the top
is finite. Profile construction compares that first edge with the observer's
kilometer value after conversion to meters using a finite, four-machine-epsilon
relative/absolute tolerance. This accepts representation-only conversion
rounding but rejects real height differences; it never moves an edge. Internal
edges select the layer above; the final edge belongs to the last layer. Partial
vertical optical depths use exact constant-layer products, with checked sums and
no query allocation, quadrature, extrapolation or clamping.
See the compiling example on `LayeredAtmosphere`.

Aerosol extinction is split by explicit SSA; vacuum is valid and its stored Optica
SSA=0 remains an undefined-ratio sentinel. Complete constituent states and totals
are validated, including public-field mutation and sum overflow. Profile columns
use local pressure plus zero additional altitude correction. Unsupported domains
and invalid numerical states return `LayeredAtmosphereError`.

This height API does not implement an infallible `Medium`. `at_origin` explicitly
binds AFFN center parameters and an absolute origin height to a borrowed Optica
`TryMedium` implementation in the local Horizontal frame. It rejects center
mismatches, nonfinite positions and out-of-domain heights/wavelengths. It never
guesses a geoid conversion from ellipsoidal to MSL height. Optica rays and path
integrators remain available; its infallible-medium optical-depth integrator is
not used to discard errors from this fallible column. These primitives do not perform angular integration,
scattering transport, surface reflection or empirical site calibration.
