#!/usr/bin/env python3
"""
Plot ERA-40 fields on the full Gaussian N80 grid (320 x 160):

  theta  potential temperature on the 2 PVU surface   theta2pvu_YYYYMM.grib
  u, v   wind at 300 hPa                               uv300_YYYYMM.grib

Give one or both files; what is drawn follows from what you give unless you
narrow it with --show.

    theta only   -> filled contours of theta
    uv only      -> filled contours of wind speed + wind vectors
    both         -> filled contours of theta + wind vectors (speed in title)

Usage
    python plot_era40.py --theta theta2pvu_199812.grib --list
    python plot_era40.py --theta theta2pvu_199812.grib --time 1998-12-15T12
    python plot_era40.py --uv uv300_199812.grib --time 1998-12-15T12 --proj npstereo
    python plot_era40.py --theta theta2pvu_199812.grib --uv uv300_199812.grib --time 1998-12-15T12 --out dec15.png
    python plot_era40.py --theta ... --uv ... --time ... --show theta      # ignore the wind file for drawing

Inputs may be GRIB (regular or reduced Gaussian grid; a reduced grid is
interpolated to the full grid on the fly) or netCDF written by
reduced_to_full.py.  The time must exist in every file given.

Requirements: python-eccodes, numpy, matplotlib; cartopy optional (maps,
polar projections); xarray + netcdf4 only for .nc inputs.
"""
import argparse
import datetime as dt
import sys

import numpy as np

# cartopy (pyproj) must be imported BEFORE eccodes with the pip wheels,
# otherwise Python aborts at exit.
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.util import add_cyclic_point
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False

import eccodes as ec

NAMES = {"theta": ("pt",), "uv": ("u", "v")}


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def _valid_time(gid):
    d, t = ec.codes_get(gid, "validityDate"), ec.codes_get(gid, "validityTime")
    return dt.datetime(d // 10000, d // 100 % 100, d % 100, t // 100, t % 100)


def _reduced_to_full(values, pl):
    nlon = int(pl.max())
    target = np.arange(nlon) * 360.0 / nlon
    out, pos = np.empty((pl.size, nlon)), 0
    for i, n in enumerate(pl):
        row = values[pos:pos + n]
        pos += n
        if n == nlon:
            out[i] = row
            continue
        src = np.arange(n) * 360.0 / n
        out[i] = np.interp(target, np.r_[src[-1] - 360, src, src[0] + 360], np.r_[row[-1], row, row[0]])
    return out


def read_grib(path, names, when=None):
    """Return lats, lons, times, {name: 2-D array} (for `when`; if None, only times)."""
    lats = lons = None
    times, fields = set(), {}
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                name = ec.codes_get(gid, "shortName")
                if name not in names:
                    continue
                t = _valid_time(gid)
                times.add(t)
                if when is None or t != when:
                    continue
                gt = ec.codes_get(gid, "gridType")
                if gt == "regular_gg":
                    nj, ni = ec.codes_get(gid, "Nj"), ec.codes_get(gid, "Ni")
                    arr = ec.codes_get_values(gid).reshape(nj, ni)
                    lo = ec.codes_get_array(gid, "distinctLongitudes")
                elif gt == "reduced_gg":
                    pl = ec.codes_get_array(gid, "pl").astype(int)
                    arr = _reduced_to_full(ec.codes_get_values(gid), pl)
                    lo = np.arange(arr.shape[1]) * 360.0 / arr.shape[1]
                else:
                    sys.exit(f"{path}: unsupported gridType {gt}")
                la = ec.codes_get_array(gid, "distinctLatitudes")
                if lats is None:
                    lats, lons = la, lo
                fields[name] = arr
            finally:
                ec.codes_release(gid)
    return lats, lons, sorted(times), fields


def read_netcdf(path, names, when=None):
    import xarray as xr
    ds = xr.open_dataset(path)
    times = sorted(dt.datetime.utcfromtimestamp(int(t) / 1e9) for t in ds.time.values.astype("datetime64[ns]").astype("int64"))
    fields = {}
    if when is not None:
        sel = ds.sel(time=np.datetime64(when))
        fields = {n: sel[n].values for n in names if n in ds}
    return ds.lat.values, ds.lon.values, times, fields


def read(path, names, when=None):
    if path.endswith(".nc"):
        return read_netcdf(path, names, when)
    return read_grib(path, names, when)


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
def make_axes(plt, proj):
    """Return fig, ax, transform-kwargs for the requested projection."""
    if not HAVE_CARTOPY:
        if proj != "global":
            print("cartopy not installed; falling back to a plain lat/lon plot")
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        return fig, ax, {}
    data_crs = ccrs.PlateCarree()
    if proj == "npstereo":
        crs, edge = ccrs.NorthPolarStereo(central_longitude=-90), 20
    elif proj == "spstereo":
        crs, edge = ccrs.SouthPolarStereo(central_longitude=0), -20
    else:
        crs, edge = ccrs.PlateCarree(central_longitude=180), None
    fig = plt.figure(figsize=(12, 7 if proj == "global" else 9))
    ax = plt.axes(projection=crs)
    if edge is not None:
        import matplotlib.path as mpath
        r = np.hypot(*crs.transform_point(0, edge, data_crs))
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)
        t = np.linspace(0, 2 * np.pi, 200)
        ax.set_boundary(mpath.Path(np.column_stack([np.sin(t), np.cos(t)]) * 0.5 + 0.5), transform=ax.transAxes)
    try:
        from cartopy.io import shapereader
        shapereader.natural_earth(resolution="110m", category="physical", name="coastline")
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), lw=0.6, color="0.25")
    except Exception as e:
        print(f"coastlines unavailable ({e.__class__.__name__}); drawing without them")
    gl = ax.gridlines(draw_labels=(proj == "global"), lw=0.3, color="0.6", alpha=0.6)
    if proj == "global":
        gl.top_labels = gl.right_labels = False
    return fig, ax, dict(transform=data_crs)


def nice_levels(arr, step):
    lo, hi = np.nanpercentile(arr, [1, 99])
    return np.arange(np.floor(lo / step) * step, np.ceil(hi / step) * step + step / 2, step)


def plot(lats, lons, f, when, show, proj="global", out=None, stride=None):
    import matplotlib
    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # region + cyclic point
    if proj == "npstereo":
        keep = (lats > 15) & (lats < 89.9)
    elif proj == "spstereo":
        keep = (lats < -15) & (lats > -89.9)
    else:
        keep = np.ones_like(lats, bool)
    lat = lats[keep]
    fld = {k: v[keep] for k, v in f.items()}
    lon = lons
    if HAVE_CARTOPY:
        for k in fld:
            fld[k], lon = add_cyclic_point(fld[k], coord=lons)

    fig, ax, tr = make_axes(plt, proj)
    title = [f"ERA-40  {when:%Y-%m-%d %H}Z"]

    # --- filled field: theta if requested, else wind speed ---------------------
    if "theta" in show:
        levels = nice_levels(fld["pt"], 5)
        cf = ax.contourf(lon, lat, fld["pt"], levels=levels, cmap="magma", extend="both", **tr)
        label = "potential temperature on the 2 PVU surface [K]"
        title.append("theta(2 PVU)")
    else:
        speed = np.hypot(fld["u"], fld["v"])
        levels = nice_levels(speed, 5)
        levels = levels[levels >= 0]
        cf = ax.contourf(lon, lat, speed, levels=levels, cmap="Blues", extend="max", **tr)
        label = "wind speed at 300 hPa [m/s]"
    cb = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05, shrink=0.7, aspect=40)
    cb.set_label(label)

    # --- wind vectors ---------------------------------------------------------
    if "uv" in show:
        if stride is None:
            stride = 6 if proj == "global" else 2
        s = slice(None, None, stride)
        qkw = dict(color="black", edgecolor="white", linewidth=0.4,
                   scale=1500, width=0.0018, headwidth=4, **tr)
        if HAVE_CARTOPY and proj != "global":
            qkw["regrid_shape"] = 35
        q = ax.quiver(lon[s], lat[s], fld["u"][s, s], fld["v"][s, s], **qkw)
        ax.quiverkey(q, 0.9, 1.03, 50, "50 m/s", labelpos="E", coordinates="axes",
                     fontproperties={"size": 8})
        title.append(f"wind 300 hPa (max {np.hypot(fld['u'], fld['v']).max():.0f} m/s)")

    ax.set_title("   ".join(title), loc="left", fontsize=11)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print("wrote", out)
    else:
        plt.show()


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--theta", help="theta2pvu_YYYYMM.grib (or .nc)")
    p.add_argument("--uv", help="uv300_YYYYMM.grib (or .nc)")
    p.add_argument("--show", nargs="+", choices=["theta", "uv"],
                   help="what to draw (default: everything the given files allow)")
    p.add_argument("--time", help="analysis time, e.g. 1998-12-15T12 (default: first in file)")
    p.add_argument("--list", action="store_true", help="list analysis times and exit")
    p.add_argument("--proj", choices=["global", "npstereo", "spstereo"], default="global")
    p.add_argument("--stride", type=int, help="plot every Nth wind vector")
    p.add_argument("--out", help="save figure to this PNG/PDF instead of showing it")
    a = p.parse_args()

    files = {k: v for k, v in (("theta", a.theta), ("uv", a.uv)) if v}
    if not files:
        p.error("give --theta and/or --uv")
    show = a.show or list(files)
    missing = [k for k in show if k not in files]
    if missing:
        p.error(f"--show {missing} but no file given for it")

    # times common to the files we will draw from
    times = None
    for k in show:
        _, _, t, _ = read(files[k], NAMES[k])
        times = set(t) if times is None else times & set(t)
    times = sorted(times)
    if not times:
        sys.exit("the files have no analysis time in common")
    if a.list:
        print(f"{len(times)} times, {times[0]:%Y-%m-%d %H}Z .. {times[-1]:%Y-%m-%d %H}Z")
        return
    when = dt.datetime.fromisoformat(a.time) if a.time else times[0]
    if when not in times:
        sys.exit(f"{when} not available; nearest: {min(times, key=lambda t: abs(t - when))}")

    lats = lons = None
    fields = {}
    for k in show:
        la, lo, _, fl = read(files[k], NAMES[k], when)
        if lats is not None and not (np.allclose(la, lats) and np.allclose(lo, lons)):
            sys.exit("theta and uv files are on different grids")
        lats, lons = la, lo
        fields.update(fl)
    need = [n for k in show for n in NAMES[k] if n not in fields]
    if need:
        sys.exit(f"fields {need} not found at {when:%Y-%m-%d %H}Z")
    print(f"grid {lons.size} x {lats.size}, {when:%Y-%m-%d %H}Z, fields {sorted(fields)}")
    plot(lats, lons, fields, when, show, proj=a.proj, out=a.out, stride=a.stride)


if __name__ == "__main__":
    main()