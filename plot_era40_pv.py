#!/usr/bin/env python3
"""
Visualise one ERA-40 PV-surface file (e4oper.an.pv.YYYYMM from NCAR GDEX d117004).

The file holds, for every 6-hourly analysis in the month, seven fields on the
+-2 PVU surface (reduced Gaussian N80 grid, GRIB1):
    pt   potential temperature [K]      pres  pressure [Pa]     z  geopotential
    u    eastward wind [m/s]            v     northward wind    q  humidity    o3

This script picks one analysis time, interpolates the reduced Gaussian points
onto a regular lat/lon grid, and draws potential temperature (filled) with
wind vectors on top.  Optionally it also writes the regridded fields to netCDF.

Usage
    python plot_era40_pv.py e4oper.an.pv.195709                    # first time in file
    python plot_era40_pv.py e4oper.an.pv.195709 --time 1957-09-15T12
    python plot_era40_pv.py e4oper.an.pv.195709 --time 1957-09-15T12 --proj npstereo --out sep15.png
    python plot_era40_pv.py e4oper.an.pv.195709 --list              # show available times
    python plot_era40_pv.py e4oper.an.pv.195709 --time 1957-09-15T12 --netcdf out.nc

Requirements:  pip install eccodes numpy scipy matplotlib   (+ cartopy, optional)
"""
import argparse
import datetime as dt
import sys

import numpy as np

# cartopy (pyproj) must be imported BEFORE eccodes: with the pip wheels, the
# reverse order makes Python abort with "double free or corruption" at exit.
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.util import add_cyclic_point
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False

import eccodes as ec
from scipy.interpolate import griddata

FIELDS = ("pt", "pres", "u", "v")  # what we read; everything else is skipped


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------
def scan_times(path):
    """Return the sorted list of analysis datetimes present in the file."""
    times = set()
    with open(path, "rb") as f:
        while True:
            g = ec.codes_grib_new_from_file(f)
            if g is None:
                break
            times.add(_valid_time(g))
            ec.codes_release(g)
    return sorted(times)


def _valid_time(g):
    d = ec.codes_get(g, "validityDate")
    t = ec.codes_get(g, "validityTime")
    return dt.datetime(d // 10000, d // 100 % 100, d % 100, t // 100, t % 100)


def read_fields(path, when):
    """Read pt/pres/u/v at datetime `when`. Returns (lats, lons, {name: values})."""
    out, lats, lons = {}, None, None
    with open(path, "rb") as f:
        while True:
            g = ec.codes_grib_new_from_file(f)
            if g is None:
                break
            try:
                name = ec.codes_get(g, "shortName")
                if name in FIELDS and _valid_time(g) == when:
                    if lats is None:  # geometry is identical for every message
                        lats = ec.codes_get_array(g, "latitudes")
                        lons = ec.codes_get_array(g, "longitudes")
                        grid = ec.codes_get(g, "gridType")
                        n = ec.codes_get(g, "N")
                        print(f"grid: {grid}, N={n}, {lats.size} points")
                    out[name] = ec.codes_get_values(g)
            finally:
                ec.codes_release(g)
    missing = [k for k in FIELDS if k not in out]
    if missing:
        sys.exit(f"fields {missing} not found at {when:%Y-%m-%d %H}Z (use --list)")
    return lats, lons, out


# ----------------------------------------------------------------------------
# Regridding reduced Gaussian -> regular lat/lon
# ----------------------------------------------------------------------------
def regrid(lats, lons, fields, dlat=1.125, dlon=1.125):
    """Linear interpolation of scattered Gaussian points onto a regular grid.

    The longitude array is padded by one period on each side so the
    interpolation wraps cleanly across 0/360.
    """
    lon_p = np.concatenate([lons - 360, lons, lons + 360])
    lat_p = np.tile(lats, 3)
    glat = np.arange(90, -90 - 1e-6, -dlat)
    glon = np.arange(0, 360, dlon)
    GLON, GLAT = np.meshgrid(glon, glat)
    pts = np.column_stack([lon_p, lat_p])
    out = {}
    for k, v in fields.items():
        zi = griddata(pts, np.tile(v, 3), (GLON, GLAT), method="linear")
        # poles lie outside the convex hull of the Gaussian points -> fill nearest
        nan = np.isnan(zi)
        if nan.any():
            zi[nan] = griddata(pts, np.tile(v, 3), (GLON[nan], GLAT[nan]), method="nearest")
        out[k] = zi
    return glat, glon, out


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def plot(glat, glon, f, when, proj="global", out=None, stride=None):
    import matplotlib
    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    speed = np.hypot(f["u"], f["v"])
    if not HAVE_CARTOPY and proj != "global":
        print("cartopy not installed; falling back to a plain lat/lon map")
        proj = "global"

    # --- restrict to the region drawn and close the longitude seam -----------
    if proj == "npstereo":
        keep = (glat > 15) & (glat < 89.9)      # pole row makes contourf degenerate
    elif proj == "spstereo":
        keep = (glat < -15) & (glat > -89.9)
    else:
        keep = np.ones_like(glat, bool)
    lat = glat[keep]
    fld = {k: v[keep] for k, v in f.items()}
    if HAVE_CARTOPY:
        lon = None
        for k in fld:
            fld[k], lon = add_cyclic_point(fld[k], coord=glon)
    else:
        lon = glon

    # --- axes ----------------------------------------------------------------
    if HAVE_CARTOPY:
        data_crs = ccrs.PlateCarree()
        if proj == "npstereo":
            crs, edge_lat = ccrs.NorthPolarStereo(central_longitude=-90), 20
        elif proj == "spstereo":
            crs, edge_lat = ccrs.SouthPolarStereo(central_longitude=0), -20
        else:
            crs, edge_lat = ccrs.PlateCarree(central_longitude=180), None
        fig = plt.figure(figsize=(12, 7 if proj == "global" else 9))
        ax = plt.axes(projection=crs)
        if edge_lat is not None:
            # circular polar cap down to edge_lat (set in projected coords: a
            # set_extent over the full longitude range degenerates in cartopy)
            import matplotlib.path as mpath
            r = np.hypot(*crs.transform_point(0, edge_lat, data_crs))
            ax.set_xlim(-r, r)
            ax.set_ylim(-r, r)
            t = np.linspace(0, 2 * np.pi, 200)
            ax.set_boundary(mpath.Path(np.column_stack([np.sin(t), np.cos(t)]) * 0.5 + 0.5),
                            transform=ax.transAxes)
        try:  # cartopy downloads Natural Earth coastlines on first use
            from cartopy.io import shapereader
            shapereader.natural_earth(resolution="110m", category="physical", name="coastline")
            ax.add_feature(cfeature.COASTLINE.with_scale("110m"), lw=0.6, color="0.25")
        except Exception as e:
            print(f"coastlines unavailable ({e.__class__.__name__}); drawing without them")
        gl = ax.gridlines(draw_labels=(proj == "global"), lw=0.3, color="0.6", alpha=0.6)
        if proj == "global":
            gl.top_labels = gl.right_labels = False
        tr = dict(transform=data_crs)
    else:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        tr = {}

    # --- potential temperature: one sequential, perceptually uniform ramp ----
    # Typical 2-PVU theta range is ~280-380 K; the levels adapt to the data.
    lo, hi = np.nanpercentile(fld["pt"], [1, 99])
    levels = np.arange(np.floor(lo / 5) * 5, np.ceil(hi / 5) * 5 + 1, 5)
    cf = ax.contourf(lon, lat, fld["pt"], levels=levels, cmap="magma", extend="both", **tr)
    cb = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05, shrink=0.7, aspect=40)
    cb.set_label("potential temperature on the 2 PVU surface [K]")

    # --- pressure of the surface (hPa) as thin contours = tropopause height --
    cs = ax.contour(lon, lat, fld["pres"] / 100.0, levels=np.arange(100, 600, 50),
                    colors="white", linewidths=0.5, alpha=0.7, **tr)
    if proj == "global":  # clabel misplaces text on the polar projections
        ax.clabel(cs, fmt="%d", fontsize=7, inline=True)

    # --- wind vectors, subsampled so they stay legible -----------------------
    if stride is None:
        stride = 6 if proj == "global" else 2
    s = slice(None, None, stride)
    qkw = dict(color="black", edgecolor="white", linewidth=0.4,
               scale=1500, width=0.0018, headwidth=4, **tr)
    if HAVE_CARTOPY and proj != "global":
        # resample onto an even grid in map coordinates so the arrows do not
        # bunch up towards the pole
        qkw["regrid_shape"] = 35
    q = ax.quiver(lon[s], lat[s], fld["u"][s, s], fld["v"][s, s], **qkw)
    ax.quiverkey(q, 0.9, 1.03, 50, "50 m/s", labelpos="E", coordinates="axes",
                 fontproperties={"size": 8})

    ax.set_title(f"ERA-40  2 PVU surface   {when:%Y-%m-%d %H}Z     "
                 f"max wind {np.nanmax(speed):.0f} m/s", loc="left", fontsize=11)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print("wrote", out)
    else:
        plt.show()


# ----------------------------------------------------------------------------
def write_netcdf(path, glat, glon, f, when):
    import xarray as xr
    ds = xr.Dataset(
        {k: (("lat", "lon"), v.astype("f4")) for k, v in f.items()},
        coords={"lat": glat, "lon": glon, "time": np.datetime64(when)},
    )
    ds["pt"].attrs.update(long_name="potential temperature on 2 PVU surface", units="K")
    ds["pres"].attrs.update(long_name="pressure of 2 PVU surface", units="Pa")
    ds["u"].attrs.update(long_name="eastward wind on 2 PVU surface", units="m s-1")
    ds["v"].attrs.update(long_name="northward wind on 2 PVU surface", units="m s-1")
    ds.to_netcdf(path)
    print("wrote", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", help="e4oper.an.pv.YYYYMM GRIB file")
    p.add_argument("--time", help="analysis time, e.g. 1957-09-15T12 (default: first in file)")
    p.add_argument("--list", action="store_true", help="list analysis times in the file and exit")
    p.add_argument("--proj", choices=["global", "npstereo", "spstereo"], default="global")
    p.add_argument("--res", type=float, default=1.125, help="output grid spacing in degrees (default 1.125 ~ N80)")
    p.add_argument("--stride", type=int, help="plot every Nth wind vector")
    p.add_argument("--out", help="save figure to this PNG/PDF instead of showing it")
    p.add_argument("--netcdf", help="also write the regridded fields to this netCDF file")
    a = p.parse_args()

    times = scan_times(a.file)
    if a.list:
        print(f"{len(times)} analysis times, {times[0]:%Y-%m-%d %H}Z .. {times[-1]:%Y-%m-%d %H}Z (6-hourly)")
        return
    when = dt.datetime.fromisoformat(a.time) if a.time else times[0]
    if when not in times:
        sys.exit(f"{when} not in file; nearest: {min(times, key=lambda t: abs(t - when))}")

    lats, lons, fields = read_fields(a.file, when)
    print(f"{when:%Y-%m-%d %H}Z  theta {fields['pt'].min():.0f}-{fields['pt'].max():.0f} K, "
          f"max |V| {np.hypot(fields['u'], fields['v']).max():.0f} m/s")
    glat, glon, reg = regrid(lats, lons, fields, a.res, a.res)
    if a.netcdf:
        write_netcdf(a.netcdf, glat, glon, reg, when)
    plot(glat, glon, reg, when, proj=a.proj, out=a.out, stride=a.stride)


if __name__ == "__main__":
    main()