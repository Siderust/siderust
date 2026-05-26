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
