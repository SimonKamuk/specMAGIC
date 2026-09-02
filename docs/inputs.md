# Inputs

* **Aerosol climatology.** This is used for calculating the effective cloud albedo. At present the code is using one from ECMWF.
* **Water vapour climatology.** See above.
* **Ozone climatology.** This is used for calculating the effective cloud albedo. At present the code is using one from ERA-Interim.
* **Land use map.** This is used for calculating the reflection from land. Currently a coarse one is being used from IGBPa.
* **Spectral aerosol lookup table.** RTM output. See the DWD article for information on how this is used and calculated.
* **Spectral water vapour lookup table.** As above.
* **Spectral ozone lookup table.** As above.
* **Spectral correction for all-weather irradiance calculation.** Used to correct from clear- to the all-sky case. Source unnown.

Be aware that if running with `--albedo MODIS`, other land maps are being downloaded and used.