# Lunar beam irradiance correction

Jones et al. (2013), Eq. 1 (https://doi.org/10.1051/0004-6361/201322433)
uses the reference lunar solid angle and one inverse-square correction.
The old helper used the actual solid angle and another distance correction,
yielding D^-4. Conservation of radiance instead gives
L_disc = E_sun A/pi and E_moon = L_disc Omega = E_sun A (R/D)^2.
The downstream phase function (sr^-1, integral one) converts this beam
irradiance to sky radiance; no further solid-angle conversion belongs there.

The new reflected_lunar_spectral_irradiance_jones2013 helper returns a typed
structural spectral irradiance. The legacy radiance-named helper retains its
Rust signature and returns the same corrected numeric beam irradiance, with
its historical unit-marker error explicitly documented. No lunar albedo or
phase calibration changes, and the mean-distance normalization is unchanged.
This avoids an unannounced public return-type break; new code should use the
typed irradiance helper, and legacy users must audit their physical conversion.

At 350000/384400 km the old response is 1.454994458; the corrected ratio is
(384400/350000)^2 = 1.2062315102. Independent tests construct disc radiance
multiplied by apparent solid angle, check distance ratios, zero illumination,
and invalid/nonfinite inputs. Downstream NSB evidence lives in
nsb-cross-validation/investigation. Shipping requires an upstream release and
a downstream lockfile update; no unpublished dependency is committed.
