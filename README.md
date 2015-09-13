# atca-1934-modelling
Spectral modelling of the source PKS 1934-638

This code is used to derive the "new" flux density model for the ATCA primary
flux density calibrator PKS 1934-638 at high frequency. It does this by adding
measurements of 1934-638 in the 3mm band.

The 3mm band measurements of 1934-638 were performed on 2012 August 12 with the ATCA
in its H75 configuration. Observations were made between 15:22 UTC and 17:59 UTC.
Being night-time observations (local time is 10 hours ahead of UTC), the weather conditions
over this period were stable, with ambient temperature varying between 6.3
and 7.5 degrees Celcius, relative humidity between 80% and 84%, and air pressure between 975
and 977 mB. Atmospheric phase stability, as monitored by an independent interferometer
looking at a 30 GHz geostationary satellite signal, never exceeded 130 microns.

During this time, 1934-638 was observed for a total of 63 minutes, between 15:40 UTC
and 17:30 UTC, and the observations covered the elevation range 32 - 44 degrees.
The flux density of 1934-638 was determined by comparing it to the planet Uranus, using the
de Pater (1990) model for its flux density. Uranus was itself observed for a total 320
seconds, at an elevation of approximately 55.2 degrees.
System temperature was monitored by regular observations of an absorbing material, kept
at ambient temperature. This absorber moves directly over the feed horn (there is one
in each antenna, all observed simultaneously), and by comparing the power seen while
observing the absorber to that seen while observing empty sky, we can calculate the
system temperature, including the atmospheric contribution. The system temperatures
slowly increased from approximately 700 K to 880 K over the duration of the measurements.

The measurements were made in two independent 2048 MHz bands, centred at 93000 MHz and
95000 MHz. All data reduction was performed on each band independently. No opacity
corrections are made to this data, nor are elevation-dependent gain changes corrected
for. The observed visibility amplitudes are scaled to compensate for changes in the
system temperature, which are linearly interpolated in time between observations of
the absorber.

Bandpass calibration was performed with the source PKS 1921-293, which had a flux
density (as measured against Uranus) of 10.73 Jy at 93 GHz, and 10.17 Jy at 95 GHz.
A phase reference calibrator - PKS 1758-651 - was observed intermittently with
1934-638, and a total of 35 minutes of integration time was obtained. Its flux density
was measured to be 0.48 +/- 0.02 Jy at 93 GHz, and 0.47 +/- 0.02 Jy at 95 GHz. The
phases measured from this calibrator were interpolated over the observations of 1934-638.

An image was made for each of the bands, and the Stokes I flux density was measured
for 1934-638 in each, assuming it to be a point source.
At 93 GHz, its flux density was measured to be 0.11 +/- 0.01 Jy,
and at 95 GHz, 0.11 +/- 0.01 Jy. For both the 93 GHz and 95 GHz images, the synthesised beam size was
measured to be 7.70 x 5.17 arcseconds.

These flux densities, measured at such high frequency, are an excellent constraint on
the flux density model for 1934-638. Our aim is to modify the Sault (2003) flux density
model to incorporate these 3mm band flux densities, while assuming that it generates correct
flux densities in the ATCA 15mm band. We make this assumption in order that flux densities referenced
against the modified model will closely match those made against the Sault (2003) model, which
has been in use for ATCA data reduction since 2003.

To do this, we first evaluated the Sault (2003) model to get a list of the Stokes I flux density of
1934-638 every 128 MHz between 10 GHz and 24 GHz. We assume a very conservative uncertainty
of 0.1 Jy for each of these flux density values (this being slightly less than 10% of the
flux density of 1934-638 between 17 and 24 GHz). The measurements of 1934-638 at 93 GHz
and 95 GHz, as listed above, are added to this list, and a first-order linear least-squares fit is
made. The resulting flux density model fit is log(S / 1 Jy) = 5.8870 - 1.3763 log (f / 1 MHz), where
S is the flux density (in Jy), f is the frequency (in MHz), and the logarithms are base 10.
