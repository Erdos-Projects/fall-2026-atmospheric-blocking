#!/usr/bin/env python3
"""
Visualise the blocking-index files written by blocking_index.py (B_YYYYMM.nc).

Default: a time-longitude heatmap (Hovmoller diagram) of B, one row per day,
diverging colours centred on B = 0, with the B = 0 contour drawn.  Several
monthly files can be given; they are concatenated in time (a whole winter,
say).

--day YYYY-MM-DD : instead, show B2d for that day as a longitude x central-
                   latitude heatmap, with the fixed CBL curve overlaid.
--cbl C          : Hovmoller of B2d at central latitude C (40..64) instead of
                   the fixed-CBL B.

Usage
    python plot_blocking.py B_199812.nc
    python plot_blocking.py B_199812.nc B_199901.nc B_199902.nc --out winter9899.png
    python plot_blocking.py B_199812.nc --day 1998-12-15
    python plot_blocking.py B_199812.nc --cbl 50

Requirements: xarray, netcdf4, matplotlib
"""
import argparse

import numpy as np
import xarray as xr


def load(files):
    """Open one or more monthly files and concatenate them along time.
    (Done by hand rather than with open_mfdataset, which needs dask.)"""
    parts = [xr.open_dataset(f).load() for f in files]
    if len(parts) == 1:
        return parts[0]
    base = parts[0]
    timed = [v for v in base.data_vars if "time" in base[v].dims]       # B, B2d, nsteps
    static = [v for v in base.data_vars if "time" not in base[v].dims]  # cbl_fixed
    ds = xr.concat([p[timed] for p in parts], dim="time").sortby("time")
    for v in static:
        ds[v] = base[v]
    ds.attrs = dict(base.attrs, source=", ".join(p.attrs.get("source", "?") for p in parts))
    return ds


def sym_limits(arr, pct=99):
    v = float(np.nanpercentile(np.abs(arr), pct))
    return -v, v


def hovmoller(ds, cbl=None, out=None):
    import matplotlib
    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    if cbl is None:
        data, label = ds.B, "B at the fixed central blocking latitude [K]"
    else:
        data, label = ds.B2d.sel(cbl=cbl), f"B at central latitude {cbl:g}N [K]"
    arr = data.values
    vmin, vmax = sym_limits(arr)
    lon = ds.lon.values
    t = ds.time.values

    # fixed, screen-friendly height: ~4.5 in for one month, capped at 7 in for a winter
    fig, ax = plt.subplots(figsize=(11, min(7.0, 3.5 + 0.035 * t.size)))
    # pcolormesh with cell edges so each day/longitude is one cell
    lon_e = np.r_[lon - 2.5, lon[-1] + 2.5]
    t_e = np.r_[t, t[-1] + np.timedelta64(1, "D")] - np.timedelta64(12, "h")
    pm = ax.pcolormesh(lon_e, t_e, arr, cmap="RdBu_r", vmin=vmin, vmax=vmax, shading="flat")
    ax.contour(lon, t, arr, levels=[0], colors="0.2", linewidths=0.6)
    # thin line at each month boundary
    months = t.astype("datetime64[M]")
    for i in np.flatnonzero(months[1:] != months[:-1]) + 1:
        ax.axhline(t_e[i], color="0.3", lw=0.8)
    cb = fig.colorbar(pm, ax=ax, pad=0.02)
    cb.set_label(label)
    ax.set_xlim(-2.5, 357.5)
    ax.set_xticks(np.arange(0, 360, 30))
    ax.set_xlabel("longitude [deg E]")
    # date ticks: every 5 days for one month, every 10 for longer spans
    step = 5 if t.size <= 40 else 10
    ax.yaxis.set_major_locator(mdates.DayLocator(bymonthday=range(1, 30, step)))   # 1, 11, 21 (never 31)
    ax.yaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.set_ylim(t_e[-1], t_e[0])           # time runs downward
    ax.set_title(f"Blocking index  {np.datetime_as_string(t[0], 'D')} .. {np.datetime_as_string(t[-1], 'D')}"
                 f"   (red: theta higher poleward of the CBL)", loc="left", fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        print("wrote", out)
    else:
        plt.show()


def day_map(ds, day, out=None):
    import matplotlib
    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sel = ds.sel(time=np.datetime64(day))
    arr = sel.B2d.values.T                 # (cbl, lon)
    vmin, vmax = sym_limits(ds.B2d.values)  # same scale across the file(s), not just this day
    lon, cbl = ds.lon.values, ds.cbl.values

    fig, ax = plt.subplots(figsize=(11, 4.5))
    lon_e = np.r_[lon - 2.5, lon[-1] + 2.5]
    cbl_e = np.r_[cbl - 0.5, cbl[-1] + 0.5]
    pm = ax.pcolormesh(lon_e, cbl_e, arr, cmap="RdBu_r", vmin=vmin, vmax=vmax, shading="flat")
    ax.contour(lon, cbl, arr, levels=[0], colors="0.2", linewidths=0.6)
    ax.plot(lon, sel.cbl_fixed.values, color="black", lw=1.2, ls="--", label="fixed CBL")
    cb = fig.colorbar(pm, ax=ax, pad=0.02)
    cb.set_label("B [K]")
    ax.set_xlim(-2.5, 357.5)
    ax.set_xticks(np.arange(0, 360, 30))
    ax.set_xlabel("longitude [deg E]")
    ax.set_ylabel("central latitude [deg N]")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_title(f"Blocking index as a function of central latitude   {day}", loc="left", fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
        print("wrote", out)
    else:
        plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="B_YYYYMM.nc file(s)")
    p.add_argument("--day", help="YYYY-MM-DD: longitude x central-latitude map of B2d for this day")
    p.add_argument("--cbl", type=float, help="Hovmoller of B2d at this central latitude (40..64)")
    p.add_argument("--out", help="save figure instead of showing it")
    a = p.parse_args()
    ds = load(a.files)
    if a.day:
        day_map(ds, a.day, a.out)
    else:
        hovmoller(ds, a.cbl, a.out)


if __name__ == "__main__":
    main()