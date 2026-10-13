#!/usr/bin/env python3
"""
Frequency of sector blocking episodes by longitude, across all winters.

For each longitude, plots

    100 * (number of SBE days, summed over winters)
        / (number of days on which SBE is possible, summed over winters)

where a winter of N days allows SBE on max(N - 4, 0) of them (the first four
days of a winter can never qualify).

Every run also stores the length-72 vectors LIB, SB, SBE (this work, %),
MHW_LIB, MHW_SB, MHW_SBE (reference, %) and CBL (deg N) in
EPISODES_DIR/longitude_vectors.nc, for use by correlate_vectors.py.

Usage
    python plot_sbe_frequency.py EPISODES_DIR                 # reads SBE_*.nc
    python plot_sbe_frequency.py EPISODES_DIR --out sbe.png
    python plot_sbe_frequency.py EPISODES_DIR --all           # also LIB and SB (denominator: all N days)
    python plot_sbe_frequency.py EPISODES_DIR --csv sbe.csv   # also write the numbers
    python plot_sbe_frequency.py EPISODES_DIR --all --cbl     # add a panel with the central blocking latitude
    python plot_sbe_frequency.py EPISODES_DIR --all --mhw     # overlay the MHW reference curves (dashed)

Requirements: numpy, xarray, netcdf4, matplotlib
"""
import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

SBE_LEAD = 4     # days at the start of a winter on which SBE cannot occur

# ---------------------------------------------------------------------------
# Reference frequencies (fraction of days) by longitude 0, 5, ..., 355,
# digitised from the solid black curves of the MHW figure (LIB / SB / SBE
# panels).  Reading accuracy is roughly +-0.005 (one pixel).
# ---------------------------------------------------------------------------
MHW = {
"LIB": np.array([
    0.307, 0.322, 0.332, 0.337, 0.332, 0.327, 0.323, 0.317, 0.316, 0.301, 0.298, 0.291,
    0.276, 0.260, 0.239, 0.237, 0.237, 0.238, 0.252, 0.251, 0.262, 0.254, 0.216, 0.168,
    0.132, 0.095, 0.080, 0.082, 0.077, 0.072, 0.077, 0.105, 0.144, 0.156, 0.156, 0.159,
    0.159, 0.161, 0.160, 0.166, 0.177, 0.167, 0.150, 0.147, 0.134, 0.129, 0.131, 0.127,
    0.128, 0.133, 0.120, 0.098, 0.075, 0.063, 0.060, 0.062, 0.082, 0.088, 0.100, 0.106,
    0.122, 0.156, 0.179, 0.191, 0.201, 0.204, 0.214, 0.216, 0.216, 0.236, 0.257, 0.275,
]),
"SB": np.array([
    0.411, 0.421, 0.430, 0.437, 0.444, 0.445, 0.441, 0.438, 0.428, 0.416, 0.410, 0.394,
    0.382, 0.372, 0.365, 0.360, 0.357, 0.355, 0.351, 0.348, 0.343, 0.323, 0.301, 0.273,
    0.230, 0.185, 0.148, 0.124, 0.122, 0.139, 0.157, 0.170, 0.184, 0.203, 0.223, 0.237,
    0.239, 0.239, 0.237, 0.234, 0.225, 0.228, 0.221, 0.213, 0.198, 0.188, 0.183, 0.171,
    0.161, 0.151, 0.142, 0.132, 0.122, 0.113, 0.108, 0.108, 0.116, 0.126, 0.150, 0.180,
    0.202, 0.222, 0.246, 0.269, 0.289, 0.300, 0.306, 0.317, 0.336, 0.356, 0.367, 0.379,
]),
"SBE": np.array([
    0.120, 0.129, 0.133, 0.136, 0.143, 0.150, 0.141, 0.138, 0.130, 0.118, 0.118, 0.110,
    0.106, 0.100, 0.098, 0.098, 0.098, 0.097, 0.092, 0.085, 0.075, 0.069, 0.057, 0.047,
    0.029, 0.023, 0.020, 0.018, 0.021, 0.026, 0.029, 0.032, 0.039, 0.040, 0.042, 0.045,
    0.042, 0.040, 0.040, 0.040, 0.037, 0.037, 0.035, 0.031, 0.024, 0.021, 0.020, 0.014,
    0.012, 0.011, 0.010, 0.006, 0.004, 0.003, 0.003, 0.006, 0.008, 0.014, 0.018, 0.027,
    0.035, 0.038, 0.044, 0.051, 0.056, 0.057, 0.055, 0.065, 0.077, 0.092, 0.097, 0.101,
]),
}


def central_blocking_latitude(lon, b_dir=None):
    """The fixed CBL per longitude: from a B_*.nc file if b_dir is given,
    otherwise from blocking_index.py (same folder as this script)."""
    if b_dir:
        f = sorted(glob.glob(os.path.join(b_dir, "B_*.nc")))
        if not f:
            sys.exit(f"no B_*.nc in {b_dir}")
        ds = xr.open_dataset(f[0])
        return ds.cbl_fixed.sel(lon=lon).values
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from blocking_index import CBL_FIXED
    except ImportError:
        sys.exit("cannot find blocking_index.py next to this script; pass --b-dir instead")
    return np.asarray(CBL_FIXED)


def tally(files, var, lead):
    """Return lon, total count per longitude, total possible days, n_winters."""
    count = possible = None
    lon = None
    for f in files:
        ds = xr.open_dataset(f)
        v = ds[var].values.astype(int)
        n = v.shape[0]
        c = v.sum(axis=0)
        p = max(n - lead, 0)
        if count is None:
            lon, count, possible = ds.lon.values, c, p
        else:
            if not np.array_equal(ds.lon.values, lon):
                sys.exit(f"{f}: longitude axis differs")
            count, possible = count + c, possible + p
    return lon, count, possible, len(files)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("episodes_dir")
    p.add_argument("--all", action="store_true", help="also plot LIB and SB frequencies")
    p.add_argument("--out", help="save figure instead of showing it")
    p.add_argument("--csv", help="write the per-longitude percentages to this CSV")
    p.add_argument("--mhw", action="store_true",
                   help="overlay the MHW reference frequencies (dashed, same colour as the matching curve)")
    p.add_argument("--cbl", action="store_true",
                   help="add a panel showing the central blocking latitude against longitude")
    p.add_argument("--b-dir", help="folder with B_*.nc, used to read the CBL for --cbl "
                                   "(default: import it from blocking_index.py)")
    a = p.parse_args()

    series = {}
    for var, lead in (("SBE", SBE_LEAD), ("SB", 0), ("LIB", 0)):
        files = sorted(glob.glob(os.path.join(a.episodes_dir, f"{var}_*.nc")))
        if not files:
            sys.exit(f"no {var}_*.nc files in {a.episodes_dir}")
        lon, count, possible, nw = tally(files, var, lead)
        if possible == 0:
            sys.exit(f"{var}: no days on which it is possible")
        series[var] = (100.0 * count / possible, count, possible, nw)
        print(f"{var}: {nw} winters, {possible} possible days per longitude, "
              f"mean {series[var][0].mean():.2f}%, max {series[var][0].max():.2f}% at {lon[series[var][0].argmax()]:.0f}E")

    # CBL (needed for --cbl and for the stored vectors); skip quietly if unavailable
    try:
        cbl = central_blocking_latitude(lon, a.b_dir)
    except SystemExit as e:
        if a.cbl:
            raise
        print(f"note: CBL not stored ({e})")
        cbl = None

    # always store the longitude vectors next to the episode files
    vec_path = os.path.join(a.episodes_dir, "longitude_vectors.nc")
    data = {v: ("lon", series[v][0].astype("f8"),
                {"units": "percent of possible days", "days": int(series[v][1].sum()),
                 "possible_days_per_lon": int(series[v][2]), "winters": int(series[v][3])})
            for v in ("LIB", "SB", "SBE")}
    for v in ("LIB", "SB", "SBE"):
        data[f"MHW_{v}"] = ("lon", 100.0 * MHW[v], {"units": "percent of days", "source": "digitised MHW figure"})
    if cbl is not None:
        data["CBL"] = ("lon", np.asarray(cbl, "f8"), {"units": "degrees_north"})
    xr.Dataset(data, coords={"lon": ("lon", lon, {"units": "degrees_east"})}).to_netcdf(vec_path)
    print("stored longitude vectors in", vec_path)

    if not a.all:
        series = {"SBE": series["SBE"]}

    if a.csv:
        cols = ["lon"] + [f"{v}_pct" for v in series] + [f"{v}_days" for v in series]
        rows = np.column_stack([lon] + [series[v][0] for v in series] + [series[v][1] for v in series])
        np.savetxt(a.csv, rows, delimiter=",", header=",".join(cols), comments="",
                   fmt=["%.0f"] + ["%.3f"] * len(series) + ["%d"] * len(series))
        print("wrote", a.csv)

    import matplotlib
    if a.out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if a.cbl:
        fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True,
                                      gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08})
    else:
        fig, ax = plt.subplots(figsize=(11, 4.5))
    colours = {"SBE": "#1f5fa8", "SB": "#e07b39", "LIB": "#6a6a6a"}
    for var, (pct, count, possible, nw) in series.items():
        # close the curve around the globe by repeating the first point at 360
        x = np.r_[lon, 360.0]
        y = np.r_[pct, pct[0]]
        ax.plot(x, y, color=colours[var], lw=2 if var == "SBE" else 1.4,
                label=f"{var} (ERA-40, this work)" if a.mhw else var)
        if a.mhw:
            ref = 100.0 * MHW[var]
            ax.plot(x, np.r_[ref, ref[0]], color=colours[var], lw=1.4, ls="--", label=f"{var} (MHW)")
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 30))
    ax.set_ylabel("% of possible days")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="0.85", lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if a.cbl:
        ax2.plot(np.r_[lon, 360.0], np.r_[cbl, cbl[0]], color="black", lw=1.4)
        ax2.set_ylabel("CBL [deg N]")
        ax2.set_ylim(np.floor(cbl.min() / 5) * 5, np.ceil(cbl.max() / 5) * 5)
        ax2.grid(axis="y", color="0.85", lw=0.6)
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
        ax2.set_xlabel("longitude [deg E]")
    else:
        ax.set_xlabel("longitude [deg E]")
    nw = series["SBE"][3]
    poss = series["SBE"][2]
    ax.set_title(f"Sector blocking episode frequency, {nw} winter{'s' if nw != 1 else ''} "
                 f"({poss} possible days per longitude)", loc="left", fontsize=10)
    if len(series) > 1 or a.mhw:
        ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=3 if a.mhw else 1)
    if not a.cbl:
        fig.tight_layout()
    if a.out:
        fig.savefig(a.out, dpi=150, bbox_inches="tight")
        print("wrote", a.out)
    else:
        plt.show()


if __name__ == "__main__":
    main()