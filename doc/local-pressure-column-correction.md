# Observer-local Rayleigh column correction

AtmosphereProfile.surface_pressure is documented as observer-local pressure.
Hydrostatic column mass p_local/g already accounts for altitude. The wrapper
now passes local pressure and zero altitude to the existing Bodhaine helper.
The low-level helper retains its numerical behavior and signature: its
pressure is a reference at height zero of the supplied exponential column.
Either reference pressure plus altitude or local pressure plus zero is valid.

At h=2.635 km, H=8 km the previous profile introduced an extra exp(-h/H)
factor (~0.71937). Tests compare barometric parameterizations and altitude
invariance at fixed local pressure. Generic profile users see increased
Rayleigh depth; NSB direct transport and airglow already avoid this error,
and the separate Jones correction reuses their shared helper.

Reference: Bodhaine et al. (1999), On Rayleigh Optical Depth Calculations,
J. Atmos. Oceanic Technol. 16, 1854-1861. The rational wavelength law is
unchanged; this patch fixes only the column semantics in the profile wrapper.
