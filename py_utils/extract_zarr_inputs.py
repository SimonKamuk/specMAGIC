"""
Extracts the clear-sky inputs for one satellite image from the zarr archives, and writes them
on the specMAGIC output grid, so that magic can use them instead of the climatologies.

    python extract_zarr_inputs.py <image_list> <config> <satellite_toml> <out_dir> <list_out>
        [--cams PATH] [--dini PATH] [--lsasaf PATH]

Every source is optional. Whatever cannot be extracted is left out with a warning, and magic
falls back to the climatology (or the land-use albedo) for it. This script therefore exits
with 0 unless its own arguments are wrong.

One NetCDF file is written per source and source time, and their paths are listed in
<list_out>, which magic reads through the ZARR_INPUTS environment variable.
"""
import os
os.environ["DASK_NUM_WORKERS"] = "1"
import sys
import argparse
import tomllib
from datetime import datetime, timezone

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

# 1 Dobson Unit of ozone in kg m-2
DU_IN_KG_M2 = 2.1415e-5

KG_M2 = ("kg m**-2", "kg m-2")
DIMENSIONLESS = ("~", "1", "")

# How far a source time may lie from the image time. The nearest available time wins,
# and a later one is only used if it is already in the store.
HOURLY_MAX_BEFORE = np.timedelta64(3, "h")
HOURLY_MAX_AFTER = np.timedelta64(30, "m")
DAILY_MAX_BEFORE = np.timedelta64(10, "D")


def warn(message):
    print(f"[WARNING] zarr inputs: {message}", file=sys.stderr, flush=True)


## ------ What the model grid and image look like ------

def config_values(config_path):
    """The value lines of magic-config.asc, as read positionally by loadConfig in src/read.cpp."""
    values = []
    with open(config_path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                values.append(line)
    return values


def output_grid(config_path, subsatellite_longitude_deg):
    """The lat/lon axes that Geography::makeCoordinates builds in headers/navigation.hpp."""
    v = config_values(config_path)
    latdim, londim = int(v[27]), int(v[28])

    # The config is parsed into floats, then widened, so mirror that rounding
    latbegin, lonbegin, dxy = (float(np.float32(x)) for x in v[29:32])

    lat = latbegin + dxy * np.arange(latdim)
    lon = lonbegin + dxy * np.arange(londim) + subsatellite_longitude_deg
    lon = np.where(lon > 180.0, lon - 360.0, lon)

    attrs = {"latbegin": latbegin, "lonbegin": lonbegin, "dxy": dxy,
             "subsatellite_longitude_deg": subsatellite_longitude_deg}
    return lat, lon, attrs


def image_time(image_list):
    """Start time of the image, parsed from its file name as Image::read_from_list does."""
    with open(image_list) as f:
        name = next(line.strip() for line in f if line.strip())

    stamp = os.path.basename(name).split("-")[2][:12]
    return np.datetime64(datetime.strptime(stamp, "%Y%m%d%H%M"), "s")


def subsatellite_longitude(satellite_toml):
    with open(satellite_toml, "rb") as f:
        return float(tomllib.load(f)["navigation"]["subsatellite_longitude_deg"])


## ------ Choosing the source time ------

def preferred_order(times, target, max_before, max_after):
    """Indices of the usable times, best first: nearest, and the earlier of two equally near."""
    times = np.asarray(times).astype("datetime64[s]")
    diff = times - target
    idx = np.flatnonzero((diff >= -max_before) & (diff <= max_after))
    return idx[np.lexsort((times[idx], np.abs(diff[idx])))]


def first_with_data(order, load, label):
    """Loads candidates in order and returns the first one with any data on the output grid."""
    for i in order:
        field = load(i)
        if np.isfinite(field).any():
            return i, field
    warn(f"{label}: no usable time within the allowed window")
    return None, None


def check_units(da, accepted, label):
    units = da.attrs.get("units", "")
    if units not in accepted:
        warn(f"{label}: unexpected units '{units}', expected one of {accepted}; skipping")
        return False
    return True


## ------ Regridding ------

def axis_weights(src, dst):
    """Lower neighbour and upper weight of each dst point on an ascending src axis."""
    i = np.clip(np.searchsorted(src, dst, side="right") - 1, 0, len(src) - 2)
    w = (dst - src[i]) / (src[i + 1] - src[i])
    inside = (dst >= src[0]) & (dst <= src[-1])
    return i, w, inside


def bilinear_regular(src_lat, src_lon, values, lat, lon):
    """Bilinear interpolation from a global lat/lon grid with 0..360 longitudes."""

    # Ascending latitude, and a wrapped column so the antimeridian interpolates
    order = np.argsort(src_lat)
    src_lat, values = src_lat[order], values[order]
    src_lon = np.append(src_lon, src_lon[0] + 360.0)
    values = np.concatenate([values, values[:, :1]], axis=1)

    iy, wy, oky = axis_weights(src_lat, lat)
    ix, wx, okx = axis_weights(src_lon, np.mod(lon, 360.0))
    wy, wx = wy[:, None], wx[None, :]

    out = ((1 - wy) * (1 - wx) * values[np.ix_(iy, ix)]
           + (1 - wy) * wx * values[np.ix_(iy, ix + 1)]
           + wy * (1 - wx) * values[np.ix_(iy + 1, ix)]
           + wy * wx * values[np.ix_(iy + 1, ix + 1)])

    out[~oky, :] = np.nan
    out[:, ~okx] = np.nan
    return out.astype(np.float32)


def bilinear_points(src_y, src_x, values, y, x):
    """Bilinear interpolation of scattered points on a regular, ascending y/x grid."""
    iy, wy, oky = axis_weights(src_y, y)
    ix, wx, okx = axis_weights(src_x, x)

    out = ((1 - wy) * (1 - wx) * values[iy, ix]
           + (1 - wy) * wx * values[iy, ix + 1]
           + wy * (1 - wx) * values[iy + 1, ix]
           + wy * wx * values[iy + 1, ix + 1])

    out[~(oky & okx)] = np.nan
    return out


## ------ The sources ------

def extract_cams(path, target, lat, lon):
    """Aerosol and ozone, plus water vapour as the backup for DINI."""
    ds = xr.open_zarr(path)
    fields = [("aod", "aod550", DIMENSIONLESS, 1.0),
              ("o3", "gtco3", KG_M2, 1.0 / DU_IN_KG_M2),
              ("h2o_cams", "tcwv", KG_M2, 1.0)]

    order = preferred_order(ds.valid_time.values, target, HOURLY_MAX_BEFORE, HOURLY_MAX_AFTER)

    results = []
    for name, var, units, scale in fields:
        label = f"cams {var}"
        if var not in ds:
            warn(f"{label}: not in {path}; skipping")
            continue
        if not check_units(ds[var], units, label):
            continue

        def load(i):
            da = ds[var].isel(valid_time=i)
            values = da.values.astype(np.float64) * scale
            return bilinear_regular(da.latitude.values, da.longitude.values, values, lat, lon)

        i, field = first_with_data(order, load, label)
        if i is None:
            continue

        results.append((name, field, ds.valid_time.values[i],
                        {"source_store": path, "source_variable": var}))
    return results


def extract_dini(path, target, lat, lon):
    """Water vapour from the most recent DINI run valid at the chosen time."""
    ds = xr.open_zarr(path)
    label = "dini wvint"
    if "wvint" not in ds:
        warn(f"{label}: not in {path}; skipping")
        return []
    if not check_units(ds.wvint, KG_M2, label):
        return []

    # Only transform the output points that can lie inside the DINI domain
    dini_lon = np.mod(ds.longitude.values + 180.0, 360.0) - 180.0
    dini_lat = ds.latitude.values
    in_box = ((lat[:, None] >= dini_lat.min()) & (lat[:, None] <= dini_lat.max())
              & (lon[None, :] >= dini_lon.min()) & (lon[None, :] <= dini_lon.max()))
    if not in_box.any():
        warn(f"{label}: the output grid does not overlap the DINI domain; skipping")
        return []

    lat2d, lon2d = np.broadcast_arrays(lat[:, None], lon[None, :])
    crs = CRS.from_wkt(ds.dini_projection.attrs["crs_wkt"])
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(
        lon2d[in_box], lat2d[in_box])

    # Every (run, step) pair, and the time it is valid at
    runs = ds.forecast_reference_time.values
    steps = ds.step.values
    if not np.issubdtype(steps.dtype, np.timedelta64):
        steps = steps.astype("timedelta64[h]")
    run_idx, step_idx = np.meshgrid(np.arange(len(runs)), np.arange(len(steps)), indexing="ij")
    run_idx, step_idx = run_idx.ravel(), step_idx.ravel()
    valid = (runs[run_idx] + steps[step_idx]).astype("datetime64[s]")

    # Nearest valid time first, and for the same valid time the most recent run
    diff = valid - target
    ok = np.flatnonzero((diff >= -HOURLY_MAX_BEFORE) & (diff <= HOURLY_MAX_AFTER))
    order = ok[np.lexsort((-runs[run_idx[ok]].astype("datetime64[s]").astype(np.int64),
                           valid[ok], np.abs(diff[ok])))]

    def load(k):
        values = ds.wvint.isel(forecast_reference_time=run_idx[k], step=step_idx[k]).values
        field = np.full(lat2d.shape, np.nan, dtype=np.float32)
        field[in_box] = bilinear_points(ds.y.values, ds.x.values, values, y, x)
        return field

    k, field = first_with_data(order, load, label)
    if k is None:
        return []

    run = runs[run_idx[k]]
    return [("h2o", field, valid[k],
             {"source_store": path, "source_variable": "wvint",
              "forecast_reference_time": str(np.datetime_as_string(run, unit="m")),
              "step_hours": int(steps[step_idx[k]] / np.timedelta64(1, "h"))})]


def extract_lsasaf(path, target, lat, lon):
    """Visible and near-infrared black-sky albedo, used as they are."""
    ds = xr.open_zarr(path)
    fields = [("albedo_vis", "AL-VI-DH"), ("albedo_nir", "AL-NI-DH")]

    # Geostationary scan angles, in radians, of every output point
    crs = CRS.from_cf(ds.geostationary.attrs)
    height = ds.geostationary.attrs["perspective_point_height"]
    lat2d, lon2d = np.broadcast_arrays(lat[:, None], lon[None, :])
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon2d, lat2d)
    x, y = x / height, y / height

    # Nearest neighbour; points off the disc come back non-finite
    sx, sy = ds.x.values, ds.y.values
    ix = np.rint((x - sx[0]) / ((sx[-1] - sx[0]) / (len(sx) - 1)))
    iy = np.rint((y - sy[0]) / ((sy[-1] - sy[0]) / (len(sy) - 1)))
    on_disc = (np.isfinite(ix) & np.isfinite(iy)
               & (ix >= 0) & (ix < len(sx)) & (iy >= 0) & (iy < len(sy)))
    if not on_disc.any():
        warn("lsasaf: the output grid is not on the LSA SAF disc; skipping")
        return []
    ix, iy = ix[on_disc].astype(int), iy[on_disc].astype(int)

    # Only read the window that the output grid covers
    x0, x1, y0, y1 = ix.min(), ix.max() + 1, iy.min(), iy.max() + 1

    # Daily product: the image's own day, or the most recent day before it
    day = target.astype("datetime64[D]").astype("datetime64[s]")
    order = preferred_order(ds.time.values, day, DAILY_MAX_BEFORE, np.timedelta64(0, "s"))

    results = []
    for name, var in fields:
        label = f"lsasaf {var}"
        if var not in ds:
            warn(f"{label}: not in {path}; skipping")
            continue

        def load(i):
            window = ds[var].isel(time=i, y=slice(y0, y1), x=slice(x0, x1)).values
            field = np.full(lat2d.shape, np.nan, dtype=np.float32)
            field[on_disc] = window[iy - y0, ix - x0]
            field[(field < 0) | (field > 1)] = np.nan
            return field

        i, field = first_with_data(order, load, label)
        if i is None:
            continue

        when = ds.time.values[i]
        if when.astype("datetime64[D]") != target.astype("datetime64[D]"):
            warn(f"{label}: no product for the image's day, using {str(when)[:10]}")

        results.append((name, field, when, {"source_store": path, "source_variable": var}))
    return results


## ------ Output ------

def write_source_files(source, results, lat, lon, grid_attrs, target, out_dir):
    """One file per source time, written atomically. Returns the paths written."""
    written = []
    for when in sorted({r[2] for r in results}):
        stamp = str(np.datetime_as_string(when, unit="m")).replace("-", "").replace("T", "").replace(":", "")
        path = os.path.join(out_dir, f"{source}-{stamp}.nc")

        data_vars = {}
        for name, field, t, attrs in results:
            if t != when:
                continue
            data_vars[name] = xr.Variable(("lat", "lon"), field, attrs=attrs)

        ds = xr.Dataset(data_vars, coords={"lat": ("lat", lat), "lon": ("lon", lon)})
        ds.attrs = {**grid_attrs, "source": source,
                    "source_time": str(np.datetime_as_string(when, unit="m")),
                    "image_time": str(np.datetime_as_string(target, unit="m")),
                    "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}

        encoding = {name: {"zlib": True, "complevel": 1, "dtype": "float32"} for name in data_vars}
        tmp = path + ".tmp"
        ds.to_netcdf(tmp, engine="netcdf4", encoding=encoding)
        os.replace(tmp, path)

        print(f" Wrote {', '.join(data_vars)} valid {ds.attrs['source_time']} to {path}")
        written.append(path)
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("image_list")
    parser.add_argument("config")
    parser.add_argument("satellite_toml")
    parser.add_argument("out_dir")
    parser.add_argument("list_out")
    parser.add_argument("--cams", default="")
    parser.add_argument("--dini", default="")
    parser.add_argument("--lsasaf", default="")
    args = parser.parse_args()

    target = image_time(args.image_list)
    lat, lon, grid_attrs = output_grid(args.config, subsatellite_longitude(args.satellite_toml))
    os.makedirs(args.out_dir, exist_ok=True)

    print(f" Extracting zarr inputs for image time {target}...")

    written = []
    for source, path, extract in (("cams", args.cams, extract_cams),
                                  ("dini", args.dini, extract_dini),
                                  ("lsasaf", args.lsasaf, extract_lsasaf)):
        if not path:
            continue
        try:
            results = extract(path, target, lat, lon)
        except Exception as e:
            warn(f"{source}: could not extract from {path} ({type(e).__name__}: {e}); "
                 "magic will use the climatology instead")
            continue
        written += write_source_files(source, results, lat, lon, grid_attrs, target, args.out_dir)

    with open(args.list_out, "w") as f:
        f.writelines(p + "\n" for p in written)


if __name__ == "__main__":
    main()
