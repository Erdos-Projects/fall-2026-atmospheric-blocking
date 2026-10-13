#!/usr/bin/env python3
"""
Local instantaneous blocking (LIB), sector blocking (SB) and sector blocking
episodes (SBE) from the monthly B_YYYYMM.nc files written by blocking_index.py.

Definitions (per day, per 5-deg longitude, periodic in longitude):
  LIB  = 1 where B > 0.
  SB   = 1 at every longitude within +-20 deg of a longitude x for which LIB = 1
         at x-5, x and x+5 (a run of three), 0 elsewhere.
  SBE  = 1 where SB = 1 on that day and on each of the previous 4 days.

Months are grouped into winters: Dec of year Y with Jan and Feb of Y+1, and
written as LIB_Y.nc, SB_Y.nc, SBE_Y.nc (Y = the December year).  Each file
holds one int8 variable (named LIB, SB or SBE) with dimensions (time, lon).
The first four days of a winter cannot be SBE because no earlier days are in
the record; days missing from the input are treated as SB = 0.

Usage
    python blocking_episodes.py B_DIR OUT_DIR
    python blocking_episodes.py B_DIR OUT_DIR --overwrite
    python blocking_episodes.py --selftest

Requirements: numpy, xarray, netcdf4
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import xarray as xr

SB_HALFWIDTH = 4      # cells of 5 deg either side of the run centre (+-20 deg)
SBE_DAYS = 5          # today + previous 4 days


# ---------------------------------------------------------------------------
def lib_sb_sbe(B):
    """B: (ndays, nlon) array (days consecutive). Returns LIB, SB, SBE as int8."""
    lib = (B > 0)
    # centre of a run of three consecutive LIB longitudes (periodic)
    run3 = lib & np.roll(lib, 1, axis=1) & np.roll(lib, -1, axis=1)
    sb = np.zeros_like(lib)
    for k in range(-SB_HALFWIDTH, SB_HALFWIDTH + 1):
        sb |= np.roll(run3, k, axis=1)
    sbe = sb.copy()
    for d in range(1, SBE_DAYS):
        prev = np.zeros_like(sb)
        prev[d:] = sb[:-d]                    # sb shifted forward by d days; first d days -> 0
        sbe &= prev
    return lib.astype("i1"), sb.astype("i1"), sbe.astype("i1")


def lib_sb_sbe_reference(B):
    """Slow, literal implementation of the description, for the self-test."""
    nd, nl = B.shape
    lib = (B > 0).astype(int)
    sb = np.zeros_like(lib)
    for t in range(nd):
        for x in range(nl):
            if lib[t, (x - 1) % nl] and lib[t, x] and lib[t, (x + 1) % nl]:
                for k in range(-4, 5):
                    sb[t, (x + k) % nl] = 1
    sbe = np.zeros_like(lib)
    for t in range(4, nd):
        for x in range(nl):
            sbe[t, x] = int(all(sb[t - d, x] for d in range(5)))
    return lib, sb, sbe


# ---------------------------------------------------------------------------
def winters(files):
    """Group B_YYYYMM.nc paths into {Y: [paths]} for Dec Y, Jan Y+1, Feb Y+1."""
    groups = {}
    for f in files:
        m = re.search(r"B_(\d{4})(\d{2})\.nc$", os.path.basename(f))
        if not m:
            continue
        y, mo = int(m.group(1)), int(m.group(2))
        if mo == 12:
            groups.setdefault(y, []).append(f)
        elif mo in (1, 2):
            groups.setdefault(y - 1, []).append(f)
    return dict(sorted(groups.items()))


def load_winter(paths):
    parts = [xr.open_dataset(p)[["B"]].load() for p in sorted(paths)]
    ds = xr.concat(parts, dim="time").sortby("time")
    t = ds.time.values
    full = np.arange(t[0], t[-1] + np.timedelta64(1, "D"), np.timedelta64(1, "D"))
    if full.size != t.size:
        print(f"  WARNING: {full.size - t.size} missing day(s) in {np.datetime_as_string(t[0], 'D')}.."
              f"{np.datetime_as_string(t[-1], 'D')}; treated as B <= 0")
        ds = ds.reindex(time=full, fill_value=np.float32(-1.0))
    return ds


def write(path, name, data, ds, desc):
    out = xr.Dataset(
        {name: (("time", "lon"), data, {"long_name": desc, "flag_values": "0, 1"})},
        coords={"time": ds.time, "lon": ds.lon},
        attrs={"definition": desc, "source": "B from blocking_index.py (fixed CBL)",
               "sbe_window_days": SBE_DAYS, "sb_halfwidth_deg": 5 * SB_HALFWIDTH})
    out.to_netcdf(path, encoding={name: {"zlib": True, "complevel": 4}})


DESC = {"LIB": "local instantaneous blocking: B > 0",
        "SB": "sector blocking: within 20 deg of a run of 3 consecutive LIB longitudes",
        "SBE": "sector blocking episode: SB on this day and the previous 4 days"}


def selftest():
    rng = np.random.default_rng(1)
    for trial in range(20):
        B = rng.standard_normal((40, 72)) + rng.uniform(-1, 1)
        if trial % 3 == 0:                      # make blocks persist: smooth in time
            B = np.cumsum(B, axis=0) / np.sqrt(np.arange(1, 41))[:, None]
        a = lib_sb_sbe(B)
        b = lib_sb_sbe_reference(B)
        for x, y, n in zip(a, b, ("LIB", "SB", "SBE")):
            assert np.array_equal(x.astype(int), y), f"{n} mismatch in trial {trial}"
    # hand-made case: a 3-run at longitudes 71, 0, 1 (wrap-around) for 6 days
    B = -np.ones((8, 72)); B[1:7, [71, 0, 1]] = 1
    lib, sb, sbe = lib_sb_sbe(B)
    assert sb[1].sum() == 9 and set(np.flatnonzero(sb[1])) == {68, 69, 70, 71, 0, 1, 2, 3, 4}
    assert sbe[:5].sum() == 0 and sbe[5].sum() == 9 and sbe[6].sum() == 9 and sbe[7].sum() == 0
    print("selftest ok: vectorised result matches the literal definition; wrap-around and 5-day window correct")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", nargs="?", help="directory with B_YYYYMM.nc")
    p.add_argument("output_dir", nargs="?")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        selftest()
        return
    if not (a.input_dir and a.output_dir):
        p.error("input_dir and output_dir are required")
    groups = winters(glob.glob(os.path.join(a.input_dir, "B_*.nc")))
    if not groups:
        sys.exit(f"no B_YYYYMM.nc files in {a.input_dir}")
    os.makedirs(a.output_dir, exist_ok=True)
    for y, paths in groups.items():
        outs = {n: os.path.join(a.output_dir, f"{n}_{y}.nc") for n in DESC}
        if all(os.path.exists(f) for f in outs.values()) and not a.overwrite:
            print(f"winter {y}-{y + 1}: outputs exist, skipped")
            continue
        months = sorted(re.search(r"(\d{6})", os.path.basename(f)).group(1) for f in paths)
        note = "" if len(paths) == 3 else f"  (WARNING: only {len(paths)} month(s): {' '.join(months)})"
        ds = load_winter(paths)
        lib, sb, sbe = lib_sb_sbe(ds.B.values)
        for n, data in zip(("LIB", "SB", "SBE"), (lib, sb, sbe)):
            write(outs[n], n, data, ds, DESC[n])
        print(f"winter {y}-{y + 1}: {ds.time.size} days, LIB {lib.mean():.1%}, SB {sb.mean():.1%}, "
              f"SBE {sbe.mean():.1%} of day-longitude cells -> {a.output_dir}/{{LIB,SB,SBE}}_{y}.nc{note}")


if __name__ == "__main__":
    main()