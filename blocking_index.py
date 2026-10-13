#!/usr/bin/env python3
"""
Pelly-Hoskins style blocking index B from coarsened ERA-40 theta (2 PVU)
files (72 longitudes x 160 Gaussian latitudes, 4 analyses per day).

For a central blocking latitude c at a given longitude,

    B = mean of theta over [c, c + 15 deg]  -  mean of theta over [c - 15 deg, c]

Each band mean is computed by linearly interpolating theta in latitude onto a
0.25 deg grid from the band's lower to its upper edge (61 points) and taking a
trapezoid-weighted average (weight 1/2 at the two ends, 1 elsewhere).  Because
that is linear in theta, each band mean is a fixed weight vector over the 160
grid latitudes; the vectors are built once and applied to every field.

Two products are written per input month, into one netCDF file:

    B    (time, lon)        c = the fixed CBL for each longitude (CBL_FIXED below)
    B2d  (time, lon, cbl)   c = 40, 41, ..., 64 deg for every longitude

Daily values are the average of the four 6-hourly values of that day.

Usage
    python blocking_index.py INPUT_DIR OUTPUT_DIR            # theta2pvu_YYYYMM.grib -> B_YYYYMM.nc
    python blocking_index.py INPUT_DIR OUTPUT_DIR --pattern 'theta2pvu_*.grib' --overwrite
    python blocking_index.py --selftest

Requirements: python-eccodes, numpy, xarray, netcdf4
"""
import argparse
import datetime as dt
import glob
import os
import re
import sys
import time

import numpy as np
import eccodes as ec

HALF_BAND = 15.0          # deg; each band is HALF_BAND wide
FINE = 0.25               # deg; fine latitude grid for the band average
CBL_2D = np.arange(40, 65, 1.0)        # 40 ... 64

# central blocking latitude for lon = 0, 5, 10, ..., 355
CBL_FIXED = np.array([
    56.5, 56.6, 57.1, 57.8, 58.2, 58.8, 58.9, 59.1, 60.0, 60.0, 60.0, 60.0,
    60.0, 59.6, 59.0, 58.9, 58.9, 58.9, 58.8, 57.8, 57.2, 55.3, 52.4, 49.1,
    46.8, 45.2, 44.4, 44.3, 44.3, 44.3, 44.4, 44.7, 45.4, 45.4, 45.4, 45.4,
    45.4, 45.5, 45.5, 45.6, 46.1, 46.4, 46.6, 47.4, 47.6, 47.6, 47.7, 47.8,
    48.5, 48.8, 48.8, 48.7, 47.8, 47.7, 46.8, 46.7, 47.4, 47.6, 47.8, 48.1,
    48.8, 49.7, 50.5, 51.3, 52.1, 52.8, 53.4, 54.2, 54.4, 55.2, 55.5, 55.5])
assert CBL_FIXED.size == 72


# ---------------------------------------------------------------------------
# weights
# ---------------------------------------------------------------------------
def band_weights(lats, lo, hi, fine=FINE):
    """Weight vector w over grid latitudes `lats` (any order) such that
    w @ theta = trapezoid-weighted mean of linearly-interpolated theta on the
    fine grid lo, lo+fine, ..., hi."""
    order = np.argsort(lats)
    asc = lats[order]
    if lo < asc[0] or hi > asc[-1]:
        raise ValueError(f"band [{lo}, {hi}] outside grid latitudes [{asc[0]:.2f}, {asc[-1]:.2f}]")
    nk = int(round((hi - lo) / fine)) + 1
    taper = np.ones(nk)
    taper[0] = taper[-1] = 0.5
    taper /= taper.sum()
    w = np.zeros(lats.size)
    for k in range(nk):
        y = lo + k * fine
        i = np.searchsorted(asc, y, side="right") - 1
        i = min(max(i, 0), asc.size - 2)
        frac = (y - asc[i]) / (asc[i + 1] - asc[i])
        w[order[i]] += taper[k] * (1.0 - frac)
        w[order[i + 1]] += taper[k] * frac
    return w


def b_weights(lats, c):
    """Weight vector giving B for central latitude c: north band minus south band."""
    return band_weights(lats, c, c + HALF_BAND) - band_weights(lats, c - HALF_BAND, c)


def build_weights(lats):
    W1 = np.stack([b_weights(lats, c) for c in CBL_FIXED])      # (72, nlat)
    W2 = np.stack([b_weights(lats, c) for c in CBL_2D])         # (25, nlat)
    return W1, W2


# ---------------------------------------------------------------------------
# GRIB reading
# ---------------------------------------------------------------------------
def read_month(path):
    """Return lats (nlat), lons (nlon), {date: [theta(nlat, nlon), ...]}."""
    lats = lons = None
    by_day = {}
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                if ec.codes_get(gid, "shortName") != "pt":
                    continue
                nj, ni = ec.codes_get(gid, "Nj"), ec.codes_get(gid, "Ni")
                if lats is None:
                    lats = ec.codes_get_array(gid, "distinctLatitudes")
                    lons = ec.codes_get_array(gid, "distinctLongitudes")
                    if ni != 72:
                        raise ValueError(f"expected 72 longitudes (coarsened file), got {ni}")
                d = ec.codes_get(gid, "validityDate")
                day = dt.date(d // 10000, d // 100 % 100, d % 100)
                by_day.setdefault(day, []).append(ec.codes_get_values(gid).reshape(nj, ni))
            finally:
                ec.codes_release(gid)
    if not by_day:
        raise ValueError("no theta (pt) messages found")
    return lats, lons, by_day


# ---------------------------------------------------------------------------
def compute_month(path):
    lats, lons, by_day = read_month(path)
    W1, W2 = build_weights(lats)
    days = sorted(by_day)
    B = np.empty((len(days), lons.size))
    B2 = np.empty((len(days), lons.size, CBL_2D.size))
    nsteps = []
    for t, day in enumerate(days):
        fields = by_day[day]
        nsteps.append(len(fields))
        theta = np.mean(fields, axis=0)                 # daily mean: linear, so B(mean) = mean(B)
        B[t] = np.einsum("jl,lj->j", W1, theta)         # per-longitude fixed CBL
        B2[t] = theta.T @ W2.T                          # (lon, cbl)
    return days, lons, B, B2, nsteps


def write_month(dst, days, lons, B, B2, nsteps, src):
    import xarray as xr
    ds = xr.Dataset(
        {"B": (("time", "lon"), B.astype("f4"),
               {"long_name": "blocking index at the fixed central blocking latitude",
                "units": "K", "description": "mean theta over [CBL, CBL+15] minus mean over [CBL-15, CBL]"}),
         "B2d": (("time", "lon", "cbl"), B2.astype("f4"),
                 {"long_name": "blocking index as a function of central latitude", "units": "K"}),
         "cbl_fixed": (("lon",), CBL_FIXED.astype("f4"),
                       {"long_name": "central blocking latitude used for B", "units": "degrees_north"}),
         "nsteps": (("time",), np.array(nsteps, "i2"), {"long_name": "6-hourly analyses averaged per day"})},
        coords={"time": np.array([np.datetime64(d) for d in days], "datetime64[ns]"),
                "lon": ("lon", lons, {"units": "degrees_east"}),
                "cbl": ("cbl", CBL_2D, {"units": "degrees_north", "long_name": "central blocking latitude"})},
        attrs={"source": os.path.basename(src),
               "method": f"band mean = trapezoid-weighted average of theta linearly interpolated to a "
                         f"{FINE} deg latitude grid over a {HALF_BAND} deg band; daily = mean of 6-hourly values",
               "created": dt.datetime.now(dt.UTC).strftime("%Y-%m-%d %H:%M UTC")})
    ds.to_netcdf(dst, encoding={k: {"zlib": True, "complevel": 4} for k in ("B", "B2d")})


# ---------------------------------------------------------------------------
def selftest():
    # synthetic Gaussian-like latitudes (N80) and a theta field linear in latitude:
    # theta = a + b*lat  ->  B = b * 15 exactly (north mean - south mean = b*HALF_BAND)
    x, _ = np.polynomial.legendre.leggauss(160)
    lats = np.degrees(np.arcsin(x))[::-1]                    # north -> south, like the files
    b = 0.8
    theta = 300 + b * lats[:, None] * np.ones((1, 72))
    W1, W2 = build_weights(lats)
    B1 = np.einsum("jl,lj->j", W1, theta)
    B2 = theta.T @ W2.T
    assert np.allclose(B1, b * HALF_BAND, atol=1e-9) and np.allclose(B2, b * HALF_BAND, atol=1e-9)
    assert np.allclose(W1.sum(axis=1), 0) and np.allclose(W2.sum(axis=1), 0)   # constant theta -> B = 0
    # brute-force check of one band mean against explicit interpolation
    c = 47.3
    fine = np.arange(c, c + HALF_BAND + 1e-9, FINE)
    prof = 300 + 3 * np.sin(np.deg2rad(4 * lats))
    interp = np.interp(fine, lats[::-1], prof[::-1])
    tap = np.ones(fine.size); tap[0] = tap[-1] = 0.5
    assert np.isclose(band_weights(lats, c, c + HALF_BAND) @ prof, (tap * interp).sum() / tap.sum())
    print("selftest ok: linear profile gives B = slope*15; constant gives 0; band mean matches explicit interpolation")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", nargs="?")
    p.add_argument("output_dir", nargs="?")
    p.add_argument("--pattern", default="*.grib")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        selftest()
        return
    if not (a.input_dir and a.output_dir):
        p.error("input_dir and output_dir are required")
    files = sorted(glob.glob(os.path.join(a.input_dir, a.pattern)))
    if not files:
        sys.exit(f"no files matching {a.pattern} in {a.input_dir}")
    os.makedirs(a.output_dir, exist_ok=True)
    t0 = time.time()
    for src in files:
        m = re.search(r"(\d{6})", os.path.basename(src))
        tag = m.group(1) if m else os.path.splitext(os.path.basename(src))[0]
        dst = os.path.join(a.output_dir, f"B_{tag}.nc")
        if os.path.exists(dst) and not a.overwrite:
            print(f"{os.path.basename(dst)} exists, skipped")
            continue
        try:
            days, lons, B, B2, nsteps = compute_month(src)
        except Exception as e:
            print(f"FAILED {os.path.basename(src)}: {e}")
            continue
        write_month(dst, days, lons, B, B2, nsteps, src)
        odd = [d for d, n in zip(days, nsteps) if n != 4]
        note = f"  (WARNING: {len(odd)} day(s) without 4 analyses)" if odd else ""
        print(f"{os.path.basename(src)}: {len(days)} days -> {dst}{note}   [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()