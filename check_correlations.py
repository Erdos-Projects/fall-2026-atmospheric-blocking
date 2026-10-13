#!/usr/bin/env python3
"""
Correlation between two of the length-72 longitude vectors stored by
plot_sbe_frequency.py in EPISODES_DIR/longitude_vectors.nc.

Vectors:  LIB  SB  SBE  MHW-LIB  MHW-SB  MHW-SBE  CBL

Usage
    python correlate_vectors.py EPISODES_DIR SBE MHW-SBE
    python correlate_vectors.py EPISODES_DIR SBE CBL --plot          # scatter + fit
    python correlate_vectors.py EPISODES_DIR --matrix                # all pairs

Reports the Pearson correlation (with the least-squares slope and intercept)
and the Spearman rank correlation.  The p-values assume 72 independent
samples; adjacent 5-deg longitudes are strongly autocorrelated, so treat them
as optimistic.

Requirements: numpy, xarray, netcdf4; scipy and matplotlib optional
"""
import argparse
import os
import sys

import numpy as np
import xarray as xr

NAMES = ["LIB", "SB", "SBE", "MHW-LIB", "MHW-SB", "MHW-SBE", "CBL"]


def load(episodes_dir):
    path = os.path.join(episodes_dir, "longitude_vectors.nc")
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run plot_sbe_frequency.py on {episodes_dir} first")
    ds = xr.open_dataset(path)
    vecs = {}
    for n in NAMES:
        key = n.replace("-", "_")
        if key in ds:
            vecs[n] = ds[key].values.astype(float)
    return ds.lon.values, vecs


def pearson(x, y):
    r = np.corrcoef(x, y)[0, 1]
    slope, intercept = np.polyfit(x, y, 1)
    n = x.size
    try:
        from scipy import stats
        p = 2 * stats.t.sf(abs(r) * np.sqrt((n - 2) / (1 - r * r)), n - 2) if abs(r) < 1 else 0.0
        rho, p_s = stats.spearmanr(x, y)
    except ImportError:
        p = rho = p_s = np.nan
    return r, p, slope, intercept, rho, p_s


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("episodes_dir")
    p.add_argument("a", nargs="?", help="first vector, one of " + " ".join(NAMES))
    p.add_argument("b", nargs="?", help="second vector")
    p.add_argument("--matrix", action="store_true", help="Pearson correlation matrix of all available vectors")
    p.add_argument("--plot", action="store_true", help="scatter plot of b against a with the fitted line")
    p.add_argument("--out", help="save the scatter plot to this file")
    args = p.parse_args()

    lon, vecs = load(args.episodes_dir)

    if args.matrix:
        names = [n for n in NAMES if n in vecs]
        M = np.corrcoef(np.stack([vecs[n] for n in names]))
        w = max(len(n) for n in names) + 1
        print(" " * w + "".join(f"{n:>9}" for n in names))
        for n, row in zip(names, M):
            print(f"{n:<{w}}" + "".join(f"{v:9.3f}" for v in row))
        return

    if not (args.a and args.b):
        p.error("give two vector names, or --matrix")
    for n in (args.a, args.b):
        if n not in NAMES:
            p.error(f"unknown vector {n!r}; choose from {' '.join(NAMES)}")
        if n not in vecs:
            sys.exit(f"{n} is not in longitude_vectors.nc (CBL is stored only when it could be found)")
    x, y = vecs[args.a], vecs[args.b]
    r, pv, slope, intercept, rho, p_s = pearson(x, y)
    print(f"{args.a} vs {args.b}  (n = {x.size} longitudes)")
    print(f"  Pearson  r = {r:+.3f}   p = {pv:.2g}   fit: {args.b} = {slope:.4f} * {args.a} + {intercept:.3f}")
    print(f"  Spearman rho = {rho:+.3f}   p = {p_s:.2g}")

    if args.plot or args.out:
        import matplotlib
        if args.out:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.5, 5))
        sc = ax.scatter(x, y, c=lon, cmap="twilight", s=28, edgecolor="white", linewidth=0.5)
        xx = np.array([x.min(), x.max()])
        ax.plot(xx, slope * xx + intercept, color="0.3", lw=1.2)
        cb = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label("longitude [deg E]")
        ax.set_xlabel(args.a)
        ax.set_ylabel(args.b)
        ax.set_title(f"r = {r:+.3f}, rho = {rho:+.3f}", loc="left", fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        if args.out:
            fig.savefig(args.out, dpi=150)
            print("wrote", args.out)
        else:
            plt.show()


if __name__ == "__main__":
    main()