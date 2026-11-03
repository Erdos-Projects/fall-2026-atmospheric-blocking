#!/usr/bin/env python3
"""
E1 (eigenbasis AR(1)) model of the 6-hourly theta-on-2-PVU field, 25.2-75.7 N x
all 72 longitudes, and its evaluation through the blocking indices.

Model (the G1/E1 procedure of PNW_Daily/exceedance/step21_ge.py, adapted):
  * climatology  clim[d, cell] = mean of theta over the TRAINING winters and the
    four 6-hourly slots of calendar day d (d = 0 .. 89, Dec 1 .. Feb 28; Feb 29
    is pooled with Feb 28);
  * anomalies    X_s = theta_s - clim[day(s)]  at every 6-hourly step s;
  * basis        C = covariance of the within-winter 6-hourly changes
                 X_{s+1} - X_s over the training winters (mean removed);
                 U = its eigenvectors, descending eigenvalue; the first K that
                 explain VAR_KEEP of the change variance are kept;
  * dynamics     g_k(s+1) = alpha_k g_k(s) + sigma_k z,  z ~ N(0, 1), independently
                 per direction k, with alpha_k the least-squares slope over the
                 training transitions and sigma_k^2 the residual variance.

Evaluation: NSAMP start days D drawn without replacement from the training
winters with D..D+5 inside the winter.  Start state = the 18 UTC field of day D;
its component outside the K kept directions is held fixed; 20 six-hourly steps
are simulated; the daily means of the simulated theta give B (Pelly-Hoskins,
fixed annual CBL, same weights as blocking_index.py) and then LIB / SB / SBE on
the full periodic 72-longitude circle (SBE = SB on all 5 days).  The data side
uses the real theta on days D+1 .. D+5 through the same B weights (checked
against the B_vals files).

Outputs (blocking/e1_theta/): e1_fit.npz, frequencies.csv, start_days.txt,
{LIB,SB,SBE}_data_e1.png, run_log.txt (stdout).

Usage:  python e1_theta_model.py [--arrays blocking/theta_arrays] [--out blocking/e1_theta]
                                 [--ntrain 30] [--var-keep 0.99] [--seed 0]
Environment: erdos_ds_era40_env (numpy, xarray, matplotlib).
"""
import argparse
import glob
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ar1_sector_model import lib_sb_sbe, MHW, ALL_LONS, SECTOR            # noqa: E402

NDAYS = 5
NSAMP = 1000
SEED = 0
VAR_KEEP = 0.99
NTRAIN = 30
SPD = 4                      # steps per day
HALF_BAND = 15.0
FINE = 0.25
# fixed (annual) central blocking latitude for lon = 0, 5, ..., 355 (blocking_index.py)
CBL_FIXED = np.array([
    56.5, 56.6, 57.1, 57.8, 58.2, 58.8, 58.9, 59.1, 60.0, 60.0, 60.0, 60.0,
    60.0, 59.6, 59.0, 58.9, 58.9, 58.9, 58.8, 57.8, 57.2, 55.3, 52.4, 49.1,
    46.8, 45.2, 44.4, 44.3, 44.3, 44.3, 44.4, 44.7, 45.4, 45.4, 45.4, 45.4,
    45.4, 45.5, 45.5, 45.6, 46.1, 46.4, 46.6, 47.4, 47.6, 47.6, 47.7, 47.8,
    48.5, 48.8, 48.8, 48.7, 47.8, 47.7, 46.8, 46.7, 47.4, 47.6, 47.8, 48.1,
    48.8, 49.7, 50.5, 51.3, 52.1, 52.8, 53.4, 54.2, 54.4, 55.2, 55.5, 55.5])


# ---------------------------------------------------------------------------
# B weights (copied from blocking_index.py so that this script has no eccodes dependency)
# ---------------------------------------------------------------------------
def band_weights(lats, lo, hi, fine=FINE):
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


def b_weight_matrix(lats):
    """(72, nlat): B[j] = W[j] @ theta[:, j]."""
    return np.stack([band_weights(lats, c, c + HALF_BAND) - band_weights(lats, c - HALF_BAND, c) for c in CBL_FIXED])


def b_from_daily_theta(theta_day, W):
    """theta_day (..., nlat, nlon) -> B (..., nlon)."""
    return np.einsum("jl,...lj->...j", W, theta_day)


# ---------------------------------------------------------------------------
def load_arrays(d):
    files = sorted(glob.glob(os.path.join(d, "theta_*.npz")))
    out = []
    for f in files:
        y = int(re.search(r"theta_(\d{4})\.npz$", f).group(1))
        z = np.load(f)
        out.append((y, z["time"], z["theta"]))
    z = np.load(files[0])
    return out, z["lat"], z["lon"]


def calendar_day(times, year):
    """0 .. 89 for Dec 1 .. Feb 28 of winter `year`; Feb 29 -> 89."""
    d = (times.astype("datetime64[D]") - np.datetime64(f"{year}-12-01")).astype(int)
    return np.minimum(d, 89)


def fit_e1(winters, clim, var_keep):
    """Change-covariance eigenbasis and per-direction AR(1) parameters over the training winters."""
    X = [theta.reshape(theta.shape[0], -1).astype(np.float64) - clim[calendar_day(t, y)] for y, t, theta in winters]
    dX = np.concatenate([x[1:] - x[:-1] for x in X])
    dX -= dX.mean(0)
    C = dX.T @ dX / dX.shape[0]
    lam, U = np.linalg.eigh(C)
    lam, U = lam[::-1], U[:, ::-1]
    share = np.cumsum(lam) / lam.sum()
    K = int(np.searchsorted(share, var_keep) + 1)
    UK = U[:, :K]
    G = [x @ UK for x in X]
    x0 = np.concatenate([g[:-1] for g in G])
    x1 = np.concatenate([g[1:] for g in G])
    alpha = (x0 * x1).sum(0) / (x0 * x0).sum(0)
    resid = ((x1 - alpha * x0) ** 2).mean(0)
    return UK, alpha, resid, lam, K, dX.shape[0]


def simulate(X0, clim, days, UK, alpha, resid, W, nlat, nlon, rng):
    """X0 (N, M) anomalies at 18 UTC of day D; days (N, NDAYS) calendar-day indices of D+1..D+5.
    Returns B (N, NDAYS, 72) from the daily means of the simulated theta."""
    N = X0.shape[0]
    g = X0 @ UK                                  # (N, K)
    R = X0 - g @ UK.T                            # held fixed
    sig = np.sqrt(resid)
    B = np.empty((N, NDAYS, nlon))
    for d in range(NDAYS):
        acc = np.zeros_like(X0)
        for _ in range(SPD):
            g = alpha * g + sig * rng.standard_normal(g.shape)
            acc += g @ UK.T + R
        theta_day = acc / SPD + clim[days[:, d]]
        B[:, d] = b_from_daily_theta(theta_day.reshape(N, nlat, nlon), W)
    return B


def frequencies(B):
    lib, sb, sbe = lib_sb_sbe(B, periodic=True)
    return {"LIB": lib.mean(axis=(0, 1)), "SB": sb.mean(axis=(0, 1)), "SBE": sbe.mean(axis=0)}


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arrays", default=os.path.join(HERE, "theta_arrays"))
    p.add_argument("--b-vals", default=os.path.join(HERE, "B_vals"))
    p.add_argument("--out", default=os.path.join(HERE, "e1_theta"))
    p.add_argument("--ntrain", type=int, default=NTRAIN)
    p.add_argument("--var-keep", type=float, default=VAR_KEEP)
    p.add_argument("--nsamp", type=int, default=NSAMP)
    p.add_argument("--seed", type=int, default=SEED)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    t0 = time.time()

    winters, lats, lons = load_arrays(a.arrays)
    nlat, nlon = lats.size, lons.size
    M = nlat * nlon
    train = winters[:a.ntrain]
    print(f"{len(winters)} winters loaded; training on the first {len(train)} ({train[0][0]}-{train[-1][0] + 1}); "
          f"grid {nlat} lat x {nlon} lon = {M} cells, 6-hourly")
    assert np.array_equal(lons, ALL_LONS)
    W = b_weight_matrix(lats)

    # ---- check: B from the daily-mean theta arrays reproduces the B_vals files
    import xarray as xr
    worst = 0.0
    for y, t, theta in train[:3] + train[-1:]:
        days = theta.reshape(-1, SPD, nlat, nlon).mean(1)
        Bth = b_from_daily_theta(days.astype(np.float64), W)
        parts = [xr.open_dataset(f)["B"].load() for f in sorted(glob.glob(os.path.join(a.b_vals, f"B_{y}12.nc")))
                 + sorted(glob.glob(os.path.join(a.b_vals, f"B_{y + 1}0[12].nc")))]
        Bfile = xr.concat(parts, dim="time").sortby("time").values
        worst = max(worst, float(np.abs(Bth - Bfile).max()))
    print(f"check: |B(theta arrays) - B_vals| max = {worst:.2e} K over 4 winters")

    # ---- climatology (training winters, 4 slots pooled)
    clim = np.zeros((90, M))
    cnt = np.zeros(90)
    for y, t, theta in train:
        cd = calendar_day(t, y)
        np.add.at(clim, cd, theta.reshape(theta.shape[0], -1))
        np.add.at(cnt, cd, 1)
    clim /= cnt[:, None]

    # ---- fit
    UK, alpha, resid, lam, K, ntrans = fit_e1(train, clim, a.var_keep)
    share = np.cumsum(lam) / lam.sum()
    print(f"fit: {ntrans} six-hourly transitions; K = {K} of {M} directions keep {share[K - 1]:.4f} of the change variance "
          f"[{time.time() - t0:.0f}s]")
    print(f"     alpha: first 5 = {np.round(alpha[:5], 3)}, min {alpha.min():.3f}, max {alpha.max():.3f}; "
          f"6-h e-folding of direction 1 = {-1 / np.log(alpha[0]):.1f} steps")
    print(f"     anomaly sd (training, all cells) {np.sqrt(np.mean([((th.reshape(th.shape[0], -1) - clim[calendar_day(t, y)]) ** 2).mean() for y, t, th in train])):.2f} K; "
          f"6-h change rms per cell {np.sqrt(lam.sum() / M):.2f} K, kept part {np.sqrt(lam[:K].sum() / M):.2f} K")
    np.savez(os.path.join(a.out, "e1_fit.npz"), U=UK.astype(np.float32), alpha=alpha, resid=resid, eig=lam, K=K,
             clim=clim.astype(np.float32), lat=lats, lon=lons, train_years=[y for y, _, _ in train])

    # ---- start days from the training winters
    starts = [(w, d) for w, (y, t, th) in enumerate(train) for d in range(th.shape[0] // SPD - NDAYS)]
    pick = rng.choice(len(starts), size=a.nsamp, replace=False)
    chosen = [starts[k] for k in pick]
    print(f"{len(starts)} admissible start days in the training winters, {a.nsamp} drawn without replacement")
    with open(os.path.join(a.out, "start_days.txt"), "w") as fh:
        for w, d in chosen:
            fh.write(f"{train[w][0]} {d}\n")

    X0 = np.empty((a.nsamp, M))
    days = np.empty((a.nsamp, NDAYS), int)
    B_data = np.empty((a.nsamp, NDAYS, nlon))
    for i, (w, d) in enumerate(chosen):
        y, t, theta = train[w]
        s = SPD * d + SPD - 1                                        # 18 UTC of day D
        cd = calendar_day(t, y)
        X0[i] = theta[s].reshape(-1) - clim[cd[s]]
        days[i] = cd[SPD * (d + 1): SPD * (d + 1 + NDAYS): SPD]
        th_days = theta[SPD * (d + 1): SPD * (d + 1 + NDAYS)].reshape(NDAYS, SPD, nlat, nlon).mean(1)
        B_data[i] = b_from_daily_theta(th_days.astype(np.float64), W)
    B_model = simulate(X0, clim, days, UK, alpha, resid, W, nlat, nlon, rng)
    print(f"simulated {a.nsamp} x {NDAYS} days [{time.time() - t0:.0f}s]")

    f_data, f_model = frequencies(B_data), frequencies(B_model)
    sec = np.isin(lons, SECTOR)
    print("\nday-by-day means over the 72 longitudes (B in K; LIB = fraction with B > 0)")
    print("  day   B data   B model   sd B data  sd B model   LIB data  LIB model")
    for d in range(NDAYS):
        print(f"  {d + 1}    {B_data[:, d].mean():7.2f}   {B_model[:, d].mean():7.2f}     {B_data[:, d].std():6.2f}     {B_model[:, d].std():6.2f}      "
              f"{(B_data[:, d] > 0).mean():.3f}     {(B_model[:, d] > 0).mean():.3f}")
    dB = B_data[:, 1:] - B_data[:, :-1]
    dM = B_model[:, 1:] - B_model[:, :-1]
    print(f"  day-to-day change of B: rms data {np.sqrt((dB ** 2).mean()):.2f} K, model {np.sqrt((dM ** 2).mean()):.2f} K; "
          f"lag-1 corr of B across days: data {np.corrcoef(B_data[:, :-1].ravel(), B_data[:, 1:].ravel())[0, 1]:.3f}, "
          f"model {np.corrcoef(B_model[:, :-1].ravel(), B_model[:, 1:].ravel())[0, 1]:.3f}")
    print("\nmean frequencies            all 72 lons          sector 330-90")
    for k in ("LIB", "SB", "SBE"):
        print(f"  {k:3s}   data {f_data[k].mean():.3f}  E1 {f_model[k].mean():.3f}     "
              f"data {f_data[k][sec].mean():.3f}  E1 {f_model[k][sec].mean():.3f}   (MHW {MHW[k].mean():.3f} / {MHW[k][sec].mean():.3f})")

    import pandas as pd
    df = pd.DataFrame({"lon": lons})
    for k in ("LIB", "SB", "SBE"):
        df[f"{k}_data"], df[f"{k}_E1"], df[f"{k}_MHW"] = f_data[k], f_model[k], MHW[k]
    df.to_csv(os.path.join(a.out, "frequencies.csv"), index=False, float_format="%.4f")
    print(df.round(3).to_string(index=False))

    # ---- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    long = {"LIB": "local instantaneous blocking (LIB)", "SB": "sector blocking (SB)",
            "SBE": "sector blocking episode (SBE)"}
    for k in ("LIB", "SB", "SBE"):
        fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="#fcfcfb")
        ax.set_facecolor("#fcfcfb")
        for lo, hi in ((0, 90), (330, 360)):
            ax.axvspan(lo, hi, color="#e6e5e0", alpha=0.6, lw=0, label="sector 330-90" if lo == 0 else None)
        ax.plot(lons, f_data[k], color="#2a78d6", marker="o", ms=4, lw=2, label="data (days D+1..D+5)")
        ax.plot(lons, f_model[k], color="#eb6834", marker="s", ms=4, lw=2, ls="--", label=f"E1 theta model (K={K})")
        ax.set_xticks(np.arange(0, 361, 30))
        ax.set_xlim(-3, 358)
        ax.set_ylim(0, None)
        ax.set_xlabel("longitude (deg E)")
        ax.set_ylabel("frequency")
        ax.set_title(f"{long[k]} frequency, DJF, {a.nsamp} five-day windows, training winters")
        ax.grid(True, color="#e6e5e0", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        png = os.path.join(a.out, f"{k}_data_e1.png")
        fig.savefig(png, dpi=150)
        plt.close(fig)
        print("wrote", png)


if __name__ == "__main__":
    main()
