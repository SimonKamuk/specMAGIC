# Inputs

Every static table is named in [`magic-config.asc`](../magic-config.asc), which is read
**positionally** by `loadConfig` in `src/read.cpp`. The order of the value lines in that
file is therefore significant: inserting or removing a line shifts every subsequent field.
The tables themselves are loaded in `src/main.cpp`, in the block marked
`Read lookup tables`, and in the `Tables::Climate` constructor.

All climatologies are long-term monthly means. Only the month of the satellite image being
processed is retained in memory; the other eleven columns are parsed and discarded.

## Climatologies

### Aerosol optical depth — `climatologies/macc-ecwmf-clim03-06rich.dat`

MACC/ECMWF aerosol climatology for 2003-06, smoothed. Read by `Tables::Climatology` into
`Tables::Aerosol::aod`, and used both for the clear-sky irradiance and (via the surface
albedo) for the effective cloud albedo.

* Grid: 720 x 360, 0.5 degrees, cell-centred (-179.75 to 179.75).
* Format: header line, then `lon lat` followed by 12 monthly AOD values, then a further
  12 columns which the code **skips** (`extra_cols_to_skip = 12`).
* Units: AOD is dimensionless, nominally at 550 nm.

The skipped columns are all `0.950` and appear to be a single-scattering albedo. The code
does not use them: single-scattering albedo and asymmetry parameter are hardcoded as
`ssa = 0.95` and `gg = 0.7` in `Tables::AerosolOptics`, on the grounds that the existing
climatology has no spatial variation in either.

### Water vapour — `climatologies/h2o_ecmwf.dat`

ERA-Interim total-column water vapour. Read by `Tables::Climatology` into
`Tables::Climate::water`, and passed to the water vapour absorption LUT.

* Grid: 1441 x 721, 0.25 degrees, node-centred (-180 to 180 inclusive).
* Units: **kg m^-2**, equivalently mm of precipitable water. Observed range in the file is
  roughly 0.2 to 73, which matches the 0-70 axis of `luts/wvcorr-l.lut`.
* This file is tracked with Git LFS. Without `git lfs pull` it is a ~130 byte pointer file,
  the read fails, and the first daytime pixel dereferences a null field pointer.

Note that the older `h2oclim.dat` (see *Unused files* below) is in **cm**, not kg m^-2, and
so is not a drop-in replacement — substituting it would put the LUT lookup out by a factor
of ten.

### Ozone — `climatologies/ERAint-o3-clim.asc`

ERA-Interim total-column ozone. Read by `Tables::Climatology` into
`Tables::Climate::ozone`, and passed to the ozone absorption LUT.

* Grid: 361 x 181, 1 degree, node-centred (-180 to 180 inclusive). Note that the config
  file records this as 360 x 181; the recorded dimensions are not used, since
  `Tables::Climatology::read` infers the grid by scanning the file.
* Units: **Dobson Units**. Observed range in the file is 157 to 457.

The low end of that range sits **below the 210 DU floor** of the `luts/o3corr-l.lut`
concentration axis. Every such cell is Antarctic (1.7 percent of the file, all between 66
and 90 degrees south), and the default output grid starts at the equator, whose minimum is
234 DU. But `magicHelpers::interpolate` does not clamp: given a concentration below the
axis floor it computes a weight greater than one, prints a `Dodgy interpolation` warning
and returns `-1000`, which is then added to the irradiance. Moving the output grid south of
about 66 degrees would trigger this.

### Land use — `climatologies/IGBPa_2006.map`

IGBP land-use classification, used to select a spectral ground albedo and a zenith-angle
reflectance correction. Read by `Tables::LandUse`.

* Grid: 2160 x 1080, 6 cells per degree, indexed by `Tables::LandUse::index`.
* Format: raw bytes, one `unsigned char` per cell, no header. The file must be exactly
  2332800 bytes or the read is rejected.
* Units: land class index, 1 to 20, dimensionless. The class table is in
  [`physics.md`](physics.md).

### MODIS BRDF — `climatologies/modis-brdf/`

Only read when the driver script is given `--albedo MODIS`; otherwise
`ModisBrdf::ModisBrdfAlbedo` stays empty and the land-use albedo is used. The ten NetCDF
files are fetched by `py_utils/download_MODIS_maps.py` and are named
`Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_<source>.nc`, where `<source>` is
`Band1` to `Band7`, `vis`, `nir` or `shortwave`.

* Grid: 7200 x 3600, 0.05 degrees, by 12 months. Only the region covered by the output
  grid is read into memory.
* Variables: `FISO_MEAN`, `FVOL_MEAN`, `FGEO_MEAN`.
* Units: stored as `short` and scaled by 0.001 to give dimensionless BRDF kernel weights;
  32767 is the fill value and marks the pixel invalid. The three weights are combined into
  a black-sky and a white-sky albedo and mixed by diffuse fraction, giving a surface albedo
  in 0-1. Anything outside that range falls back to the land-use albedo.

## Lookup tables

### Clear-sky aerosol RTM — `luts/magic-clear_spectral.lut`

Radiative transfer model output, read by `Tables::RTM`. This is the core clear-sky table.
See the DWD article for how it is generated and used.

* Shape: 2 `gg` x 3 `ssa` x 11 `aod` x 32 Kato bands. The file is six blocks, each opened
  by a `ssa= <value> gg= <value>` line, and each block is 11 x 32 rows of seven columns:
  `aod lambda Im gtau ag btau ab`.
* Axes: `ssa` in {0.70, 0.85, 1.00}, `gg` in {0.60, 0.78}, `aod` in {0.00, 0.10, 0.20,
  0.30, 0.45, 0.60, 0.80, 1.00, 1.20, 1.50, 2.00}.
* Units: `aod` dimensionless; `lambda` in nm; `Im` is a fitted Lambert-Beer prefactor with
  dimensions of W m^-2 per Kato band; `gtau` and `btau` are dimensionless
  optical-depth-like terms (negative) for the global and beam components; `ag` and `ab` are
  dimensionless exponents applied to the cosine of the solar zenith angle.

Interpolation is linear in `aod` and `ssa` but **nearest-neighbour in `gg`**.

### Absorber corrections — `luts/wvcorr-l.lut` and `luts/o3corr-l.lut`

Water vapour and ozone absorption corrections. Both are read by
`Tables::Correction::Absorber` and share a format: a header line, then `nrows x ncols`
records of six columns, `concen lambda delta_Ig ag delta_Ib ab`.

* Shape: 18 x 32 for water vapour, 8 x 32 for ozone.
* Axes: water vapour concentration 0, 2.5, ... 15 then 20, 25, ... 70 kg m^-2; ozone
  concentration 210, 255, ... 525 DU.
* Units: `lambda` in nm; `delta_Ig` and `delta_Ib` are **negative** irradiance losses in
  W m^-2 per Kato band, for the global and beam components, added to the clear-sky result;
  `ag` and `ab` are dimensionless exponents applied to the cosine of the solar zenith
  angle.

The finer sampling of the water vapour axis at low concentrations reflects the logarithmic
saturation of water vapour absorption.

### Cloud spectral correction — `luts/lambdacor.lut`

Corrects the broadband cloud transmission to a wavelength-dependent transmission, i.e.
from the clear-sky to the all-sky case. Read by `Tables::Correction::Spectral`.

* Format: header line; then one row of six clear-sky index `k*` values corresponding to six
  cloud optical depths; then 32 rows of `lambda` plus six correction factors.
* Units: `k*` and the correction factors are dimensionless, the latter being multiplicative
  on the per-band GHI; `lambda` in nm.
* Coverage: the file's own header notes that the corrections are only calculated over
  317.3-2638.5 nm. The rows outside that range repeat the nearest computed value.

The file has 34 non-empty lines, all of which are read; the remainder is trailing blank
lines.

### Spectral ground albedo — `luts/reflec.lut`

Surface albedo per land class, used whenever MODIS BRDF albedo is unavailable or
out of range. Read by `Tables::GroundAlbedo`.

* Shape: header line, then 20 rows (one per land class) of 32 columns (one per Kato band).
* Units: albedo, dimensionless, 0 to 1.

The value is multiplied by a land-class zenith-angle correction from
`Reflectivity::zenithCorrection` before use.

## Unused files

The following are present in the repository but are not referenced by the config and are
never read. They are retained from the legacy DWD code.

* `climatologies/h2oclim.dat` — NCEP long-term monthly water vapour, 144 x 73, in **cm**.
  Superseded by `h2o_ecmwf.dat`; see the units warning above.
* `climatologies/landuse.dat` — superseded by `IGBPa_2006.map`.
* `climatologies/CVS/` — leftover CVS metadata.
* `climatologies/readme_aerosol_data.txt` — documentation of alternative aerosol
  climatologies and their grid dimensions.

Several config fields are likewise parsed but unused: the recorded climatology dimensions
`xadim`, `yadim`, `xhdim`, `yhdim`, `xo3dim` and `yo3dim` (only `xadim` and `yadim` are
sanity-checked for positivity), and `iconflag` and `iconres`.

## A note on the climatology grids

`Tables::Climatology::read` derives the grid spacing as `ddeg = 360 / nlon`. This is exact
for the aerosol climatology, whose longitudes are cell-centred and so do not repeat the
antimeridian. It is slightly wrong for the water vapour and ozone climatologies, whose
longitudes run from -180 to 180 inclusive: `nlon` is 1441 and 361 rather than 1440 and 360,
giving a `ddeg` around 0.07 and 0.3 percent too small respectively. The resulting
interpolation weights are correspondingly off, which is small but not zero.
