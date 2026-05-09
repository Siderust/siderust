# `spectra`

## Scientific scope

Many astronomical response functions and atmospheric quantities are naturally
functions of wavelength but are distributed in practice as sampled tables:
filter passbands, throughput curves, irradiance spectra, ozone transmission,
and sky radiance. This module provides the typed container and interpolation
rules for those one-dimensional sampled spectral functions.

The scientific assumption is intentionally simple and explicit: a spectrum is
represented by strictly monotonic samples and is interpreted between samples by
piecewise-linear interpolation. That choice matches common synthetic-photometry
practice and keeps data provenance attached to the exact sampled curve.

## Technical scope

What you will find here:

- `sampled.rs`: the `SampledSpectrum<X, Y, S>` container.
- `interp.rs`: interpolation and out-of-range policies.
- `integrate.rs`: typed integration helpers.
- `loaders/`: file-ingestion utilities for sampled spectra.
- `passbands/`: curated built-in filter datasets such as Bessell 1990.
- `algo.rs`: lower-level untyped kernels for parity with existing pipelines.
- `provenance.rs`, `error.rs`: provenance and error types re-exported by the
  top-level module.

This module is available behind the `spectra` feature and is the 1-D spectral
counterpart to the higher-rank `tables` module.

## References

- Bessell, M. S. (1990). "UBVRI passbands." Publications of the Astronomical
  Society of the Pacific 102, 1181-1199.
- SVO Filter Profile Service. Spanish Virtual Observatory.
  <https://svo2.cab.inta-csic.es/theory/fps/>
- European Southern Observatory. SkyCalc and Advanced Cerro Paranal Sky Model.
  <https://www.eso.org/sci/software/pipelines/skytools/skycalc>
