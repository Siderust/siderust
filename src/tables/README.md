# `tables`

## Scientific scope

Some astronomical reference products are not naturally one-dimensional spectra
but multidimensional lookup grids: sky-brightness intensity as a function of
wavelength and geometry, atmospheric transmission as a function of airmass and
wavelength, or any other tabulated model whose authoritative form is a sampled
grid. This module generalizes the sampled-spectrum idea to one-, two-, and
three-dimensional typed tables.

The implemented interpolants are piecewise multilinear. The module is designed
to preserve the sampled-domain interpretation of the source data rather than to
impose a global analytic fit.

## Technical scope

What you will find here:

- `grid1d.rs`, `grid2d.rs`, `grid3d.rs`: typed linear, bilinear, and trilinear
  grid containers.
- `algo.rs`: lower-level interpolation kernels and axis-validation helpers.
- `TableError`: the shared error taxonomy.
- Re-exports of `OutOfRange`, `Provenance`, and `DataSource` so that grid
  interpolation and reproducibility metadata remain aligned with `spectra`.

This module is available behind the `tables` feature and is intended for
authoritative tabulated products whose original structure should remain visible
to downstream scientific code.

## References

- Leinert, Ch., Bowyer, S., Haikala, L. K., et al. (1998). "The 1997 reference
  of diffuse night sky brightness." Astronomy and Astrophysics Supplement
  Series 127, 1-99.
- Noll, S., Kausch, W., Barden, M., et al. (2012). "An atmospheric radiation
  model for Cerro Paranal. I. The optical spectral range." Astronomy and
  Astrophysics 543, A92.
- European Southern Observatory. SkyCalc - ESO Sky Model Calculator.
  <https://www.eso.org/sci/software/pipelines/skytools/skycalc>
