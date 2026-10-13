#!/usr/bin/env python3
"""
Coarsen ERA-40 theta files from 320 longitudes (1.125 deg, regular Gaussian
N80) to 72 longitudes (0, 5, 10, ..., 355 deg), keeping the 160 Gaussian
latitudes unchanged.

Method (per latitude row, per message):
  1. linearly interpolate the 320 values onto a 0.25 deg grid (periodic in
     longitude),
  2. for each target longitude L, average the fine-grid values at
     L-2.5, L-2.25, ..., L+2.5 (21 points) with trapezoidal weights:
     1/2 at the two ends, 1 elsewhere, divided by the total weight 20.

Both steps are linear, so they are combined once into a 72 x 320 weight
matrix W and every row is coarsened as  row @ W.T.  A direct (explicit
fine-grid) implementation is kept for --selftest.

Every .grib file in the input directory is processed; every message in each
file is coarsened and written, with the same file name, to the output
directory.  GRIB metadata (parameter, level, date/time) is copied from the
source message; the grid becomes regular Gaussian with Ni=72, 5 deg spacing.

Usage
    python coarsen_theta_lon.py INPUT_DIR OUTPUT_DIR
    python coarsen_theta_lon.py INPUT_DIR OUTPUT_DIR --pattern 'theta2pvu_*.grib' --overwrite
    python coarsen_theta_lon.py --selftest

Requirements: python-eccodes, numpy
"""
import argparse
import glob
import os
import sys
import time

import numpy as np
import eccodes as ec

NLON_IN = 320
NLON_OUT = 72
DLON_OUT = 360.0 / NLON_OUT        # 5 deg
FINE = 0.25                        # fine grid spacing, deg
HALF_WIDTH = 2.5                   # averaging window = L +- HALF_WIDTH


# ---------------------------------------------------------------------------
# weights
# ---------------------------------------------------------------------------
def interp_weights(lon, nlon_in):
    """Periodic linear interpolation at longitude `lon` from a regular grid of
    nlon_in points starting at 0: returns (i0, i1, w1) with value =
    (1-w1)*v[i0] + w1*v[i1]."""
    dlon = 360.0 / nlon_in
    x = (lon % 360.0) / dlon
    i0 = int(np.floor(x)) % nlon_in
    i1 = (i0 + 1) % nlon_in
    return i0, i1, x - np.floor(x)


def build_matrix(nlon_in=NLON_IN, nlon_out=NLON_OUT, fine=FINE, half=HALF_WIDTH):
    """W (nlon_out x nlon_in): coarse = W @ fine_row_of_source."""
    nk = int(round(2 * half / fine)) + 1            # 21 fine points per window
    taper = np.ones(nk)
    taper[0] = taper[-1] = 0.5
    taper /= taper.sum()                            # total weight 20 -> normalised
    W = np.zeros((nlon_out, nlon_in))
    for j in range(nlon_out):
        L = j * 360.0 / nlon_out
        for k in range(nk):
            f = L - half + k * fine
            i0, i1, w1 = interp_weights(f, nlon_in)
            W[j, i0] += taper[k] * (1.0 - w1)
            W[j, i1] += taper[k] * w1
    assert np.allclose(W.sum(axis=1), 1.0)
    return W


def coarsen_explicit(row, nlon_out=NLON_OUT, fine=FINE, half=HALF_WIDTH):
    """Reference implementation: interpolate to the fine grid, then average."""
    n = row.size
    src = np.arange(n) * 360.0 / n
    fine_lon = np.arange(0, 360, fine)
    fine_val = np.interp(fine_lon, np.r_[src, 360.0], np.r_[row, row[0]])   # periodic
    nk = int(round(2 * half / fine)) + 1
    taper = np.ones(nk)
    taper[0] = taper[-1] = 0.5
    out = np.empty(nlon_out)
    for j in range(nlon_out):
        L = j * 360.0 / nlon_out
        idx = (np.round((L - half) / fine).astype(int) + np.arange(nk)) % fine_lon.size
        out[j] = np.sum(taper * fine_val[idx]) / taper.sum()
    return out


# ---------------------------------------------------------------------------
# GRIB processing
# ---------------------------------------------------------------------------
def coarsen_file(src, dst, W):
    nmsg = 0
    with open(src, "rb") as fi, open(dst + ".part", "wb") as fo:
        while True:
            gid = ec.codes_grib_new_from_file(fi)
            if gid is None:
                break
            try:
                gt, ni, nj = ec.codes_get(gid, "gridType"), ec.codes_get(gid, "Ni"), ec.codes_get(gid, "Nj")
                if gt != "regular_gg" or ni != W.shape[1]:
                    raise ValueError(f"{os.path.basename(src)}: expected regular_gg with Ni={W.shape[1]}, "
                                     f"got {gt} Ni={ni}")
                field = ec.codes_get_values(gid).reshape(nj, ni)
                coarse = field @ W.T                                   # (nj, 72)
                ec.codes_set(gid, "Ni", W.shape[0])
                ec.codes_set(gid, "iDirectionIncrementInDegrees", 360.0 / W.shape[0])
                ec.codes_set(gid, "longitudeOfFirstGridPointInDegrees", 0.0)
                ec.codes_set(gid, "longitudeOfLastGridPointInDegrees", 360.0 - 360.0 / W.shape[0])
                ec.codes_set_values(gid, coarse.ravel())
                ec.codes_write(gid, fo)
                nmsg += 1
            finally:
                ec.codes_release(gid)
    os.replace(dst + ".part", dst)
    return nmsg


def selftest():
    W = build_matrix()
    rng = np.random.default_rng(0)
    lon = np.arange(NLON_IN) * 360.0 / NLON_IN
    rows = [np.full(NLON_IN, 300.0),                                  # constant
            300 + 20 * np.cos(np.deg2rad(lon)),                       # wavenumber 1
            300 + 5 * np.sin(np.deg2rad(13 * lon + 40)),              # short wave
            rng.standard_normal(NLON_IN)]                             # noise
    for r in rows:
        a, b = r @ W.T, coarsen_explicit(r)
        assert np.allclose(a, b, atol=1e-10), "matrix and explicit versions differ"
    assert np.allclose(rows[0] @ W.T, 300.0)                          # constant preserved
    # 0 deg and 180 deg windows wrap/centre correctly for the wavenumber-1 field
    c = rows[1] @ W.T
    assert c[0] > c[36] and abs(c[18] - 300) < 1e-6 and abs(c[54] - 300) < 1e-6
    print("selftest ok: matrix == explicit fine-grid average; rows sum to 1; wrap-around correct")
    print(f"each output point draws on {int((W[0] > 0).sum())} source points, "
          f"weights at 0 deg: {np.round(W[0][W[0] > 0], 4)}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", nargs="?")
    p.add_argument("output_dir", nargs="?")
    p.add_argument("--pattern", default="*.grib", help="glob for input files (default *.grib)")
    p.add_argument("--overwrite", action="store_true", help="redo files that already exist in output_dir")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        selftest()
        return
    if not (a.input_dir and a.output_dir):
        p.error("input_dir and output_dir are required")
    if os.path.abspath(a.input_dir) == os.path.abspath(a.output_dir):
        p.error("output_dir must differ from input_dir")

    files = sorted(glob.glob(os.path.join(a.input_dir, a.pattern)))
    if not files:
        sys.exit(f"no files matching {a.pattern} in {a.input_dir}")
    os.makedirs(a.output_dir, exist_ok=True)
    W = build_matrix()
    t0, done, skipped = time.time(), 0, 0
    for src in files:
        dst = os.path.join(a.output_dir, os.path.basename(src))
        if os.path.exists(dst) and not a.overwrite:
            skipped += 1
            continue
        try:
            n = coarsen_file(src, dst, W)
        except Exception as e:
            print(f"FAILED {os.path.basename(src)}: {e}")
            if os.path.exists(dst + ".part"):
                os.remove(dst + ".part")
            continue
        done += 1
        print(f"{os.path.basename(src)}: {n} messages -> {dst}   [{time.time() - t0:.0f}s]", flush=True)
    print(f"done: {done} files written, {skipped} skipped (already present), {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()