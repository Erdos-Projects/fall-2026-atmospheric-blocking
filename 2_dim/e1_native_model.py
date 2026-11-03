#!/usr/bin/env python3
"""
Daily-step E1 (eigenbasis AR(1)) model of daily-mean theta on 2 PVU on the native
ERA-40 grid (320 lon x 1.125 deg, Gaussian latitudes cropped to 25.2-75.7 N:
46 x 320 = 14720 cells), evaluated through LIB / SB / SBE.

  clim[d, cell]  mean over the training winters of calendar day d (Feb 29 -> Feb 28)
  X_t            daily-mean theta - clim                     (anomaly)
  C              covariance of X_{t+1} - X_t over the training transitions (mean removed),
                 obtained through the n x n Gram matrix; ALL non-null directions kept
  g_k            independent AR(1) per direction: alpha_k least squares, sigma_k^2 residual variance

Evaluation: NSAMP start days D from the training winters with D..D+5 in the winter;
start state = daily-mean anomaly of day D, its component outside the basis held fixed;
5 daily steps; simulated theta -> 5-deg longitudes with the coarsen_theta_long weights
-> B (fixed annual CBL) -> LIB / SB / SBE on the periodic 72-longitude circle.
Data side: the real daily theta of D+1..D+5 through the same coarsening and weights
(checked against B_vals).

Outputs in blocking/e1_native/: frequencies.csv, start_days.txt, e1_native_fit.npz,
{LIB,SB,SBE}_data_e1native.png, run_log.txt.
Usage: python e1_native_model.py [--ntrain 30] [--nsamp 1000] [--seed 0] [--lat 25 76]
Environment: erdos_ds_era40_env.
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from ar1_sector_model import lib_sb_sbe, MHW, ALL_LONS, SECTOR          # noqa: E402
from e1_theta_model import b_weight_matrix, b_from_daily_theta             # noqa: E402
from coarsen_theta_long import build_matrix                                # noqa: E402

NDAYS, NSAMP, SEED, NTRAIN = 5, 1000, 0, 30


def calday(times, year):
    return np.minimum((times.astype("datetime64[D]") - np.datetime64(f"{year}-12-01")).astype(int), 89)


def load(d, lat_range):
    out = None
    for f in sorted(glob.glob(os.path.join(d, "theta_daily_*.npz"))):
        z = np.load(f); y = int(os.path.basename(f)[-8:-4])
        if out is None:
            sel = (z["lat"] >= lat_range[0]) & (z["lat"] <= lat_range[1])
            lats, lons = z["lat"][sel], z["lon"]; out = []
        out.append((y, z["time"], z["theta"][:, sel]))
    return out, lats, lons


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arrays", default=os.path.join(HERE, "theta_native_daily"))
    p.add_argument("--b-vals", default=os.path.join(HERE, "B_vals"))
    p.add_argument("--out", default=os.path.join(HERE, "e1_native"))
    p.add_argument("--lat", type=float, nargs=2, default=(25.0, 76.0))
    p.add_argument("--ntrain", type=int, default=NTRAIN)
    p.add_argument("--nsamp", type=int, default=NSAMP)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--cbl-mask", action="store_true",
                   help="model only the cells that carry non-zero weight in some B value (the CBL +-15 deg bands)")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    t0 = time.time()

    winters, lats, lons = load(a.arrays, a.lat)
    nlat, nlon = lats.size, lons.size
    M = nlat * nlon
    train = winters[:a.ntrain]
    print(f"{len(winters)} winters; training on the first {len(train)} ({train[0][0]}-{train[-1][0] + 1}); "
          f"grid {nlat} lat ({lats.min():.2f}-{lats.max():.2f} N) x {nlon} lon = {M} cells, daily means")
    Wlon = build_matrix(nlon_in=nlon, nlon_out=72)              # (72, 320)
    Wb = b_weight_matrix(lats)                                   # (72, nlat)
    Wtot = (Wb[:, :, None] * Wlon[:, None, :]).reshape(72, -1)   # (72, nlat*nlon): B = theta_flat @ Wtot.T
    if a.cbl_mask:
        cells = np.flatnonzero((Wtot != 0).any(0))
        mask2d = np.zeros((nlat, nlon), bool); mask2d.flat[cells] = True
        rows = mask2d.sum(0)
        print(f"CBL mask: {cells.size} of {M} cells ({cells.size / M:.1%}); rows kept per longitude {rows.min()}-{rows.max()} "
              f"(latitude spans {lats[mask2d.any(1)].min():.1f}-{lats[mask2d.any(1)].max():.1f} N overall)")
        np.save(os.path.join(a.out, "cbl_mask.npy"), mask2d)
    else:
        cells = np.arange(M)
    Wc = Wtot[:, cells]
    winters = [(y, t, th.reshape(th.shape[0], -1)[:, cells]) for y, t, th in winters]
    train = winters[:a.ntrain]
    M = cells.size

    def B_of(theta_flat):
        """theta (..., M cells) -> B (..., 72)."""
        return theta_flat @ Wc.T

    # ---- check against B_vals
    import xarray as xr
    worst = 0.0
    for y, t, th in train[:2] + train[-1:]:
        Bth = B_of(th.astype(np.float64))
        fs = sorted(glob.glob(os.path.join(a.b_vals, f"B_{y}12.nc"))) + sorted(glob.glob(os.path.join(a.b_vals, f"B_{y + 1}0[12].nc")))
        Bf = xr.concat([xr.open_dataset(f)["B"].load() for f in fs], dim="time").sortby("time").values
        worst = max(worst, float(np.abs(Bth - Bf).max()))
    print(f"check: |B(native daily theta -> 5 deg) - B_vals| max = {worst:.2e} K over 3 winters")

    # ---- climatology and anomalies
    clim = np.zeros((90, M)); cnt = np.zeros(90)
    for y, t, th in train:
        cd = calday(t, y); np.add.at(clim, cd, th.astype(np.float64)); np.add.at(cnt, cd, 1)
    clim /= cnt[:, None]
    X = [th.astype(np.float64) - clim[calday(t, y)] for y, t, th in train]
    D = np.concatenate([x[1:] - x[:-1] for x in X]); mu = D.mean(0); D -= mu
    n = D.shape[0]
    lam, V = np.linalg.eigh(D @ D.T / n); lam, V = lam[::-1], V[:, ::-1]
    pos = lam > lam[0] * 1e-10; lam, V = lam[pos], V[:, pos]
    U = (D.T @ V) / np.sqrt(n * lam); del D, V
    r = lam.size
    G = [x @ U for x in X]
    x0 = np.concatenate([g[:-1] for g in G]); x1 = np.concatenate([g[1:] for g in G])
    alpha = (x0 * x1).sum(0) / (x0 * x0).sum(0); resid = ((x1 - alpha * x0) ** 2).mean(0); del G, x0, x1
    print(f"fit: {n} daily transitions, rank {r} (all kept)  [{time.time() - t0:.0f}s]")
    print(f"     alpha: first 5 = {np.round(alpha[:5], 3)}, min {alpha.min():.3f}, max {alpha.max():.3f}; "
          f"1-day e-folding of direction 1 = {-1 / np.log(alpha[0]):.1f} days")
    print(f"     anomaly sd (training, modelled cells) {np.sqrt(np.mean([(x ** 2).mean() for x in X])):.2f} K; "
          f"daily change rms per cell {np.sqrt(lam.sum() / M):.2f} K")
    np.savez(os.path.join(a.out, "e1_native_fit.npz"), U=U.astype(np.float32), alpha=alpha, resid=resid, eig=lam,
             clim=clim.astype(np.float32), lat=lats, lon=lons, train_years=[y for y, _, _ in train])

    # ---- start days
    starts = [(w, d) for w, (y, t, th) in enumerate(train) for d in range(th.shape[0] - NDAYS)]
    pick = rng.choice(len(starts), size=a.nsamp, replace=False); chosen = [starts[k] for k in pick]
    print(f"{len(starts)} admissible start days in the training winters, {a.nsamp} drawn without replacement")
    with open(os.path.join(a.out, "start_days.txt"), "w") as fh:
        for w, d in chosen:
            fh.write(f"{train[w][0]} {d}\n")
    X0 = np.stack([X[w][d] for w, d in chosen])                                           # (N, M)
    days = np.stack([calday(train[w][1], train[w][0])[d + 1: d + 1 + NDAYS] for w, d in chosen])
    B_data = np.stack([B_of(train[w][2][d + 1: d + 1 + NDAYS].astype(np.float64)) for w, d in chosen])
    del X

    # ---- simulate
    g = X0 @ U; R = X0 - g @ U.T; sig = np.sqrt(resid)
    B_model = np.empty((a.nsamp, NDAYS, 72))
    for d in range(NDAYS):
        g = alpha * g + sig * rng.standard_normal(g.shape)
        B_model[:, d] = B_of(g @ U.T + R + clim[days[:, d]])
    print(f"simulated {a.nsamp} x {NDAYS} days  [{time.time() - t0:.0f}s]")

    def freqs(B):
        lib, sb, sbe = lib_sb_sbe(B, periodic=True)
        return {"LIB": lib.mean(axis=(0, 1)), "SB": sb.mean(axis=(0, 1)), "SBE": sbe.mean(axis=0)}
    f_data, f_model = freqs(B_data), freqs(B_model)
    sec = np.isin(lons[:0].tolist() + list(ALL_LONS), SECTOR) if False else np.isin(ALL_LONS, SECTOR)
    print("\nday-by-day means over the 72 longitudes (B in K; LIB = fraction with B > 0)")
    print("  day   B data   B model   sd B data  sd B model   LIB data  LIB model")
    for d in range(NDAYS):
        print(f"  {d + 1}    {B_data[:, d].mean():7.2f}   {B_model[:, d].mean():7.2f}     {B_data[:, d].std():6.2f}     {B_model[:, d].std():6.2f}      "
              f"{(B_data[:, d] > 0).mean():.3f}     {(B_model[:, d] > 0).mean():.3f}")
    dB, dM = B_data[:, 1:] - B_data[:, :-1], B_model[:, 1:] - B_model[:, :-1]
    print(f"  day-to-day change of B: rms data {np.sqrt((dB ** 2).mean()):.2f} K, model {np.sqrt((dM ** 2).mean()):.2f} K; "
          f"lag-1 corr: data {np.corrcoef(B_data[:, :-1].ravel(), B_data[:, 1:].ravel())[0, 1]:.3f}, "
          f"model {np.corrcoef(B_model[:, :-1].ravel(), B_model[:, 1:].ravel())[0, 1]:.3f}")
    print("\nmean frequencies            all 72 lons          sector 330-90")
    for k in ("LIB", "SB", "SBE"):
        print(f"  {k:3s}   data {f_data[k].mean():.3f}  E1 {f_model[k].mean():.3f}     "
              f"data {f_data[k][sec].mean():.3f}  E1 {f_model[k][sec].mean():.3f}   (MHW {MHW[k].mean():.3f} / {MHW[k][sec].mean():.3f})")
    import pandas as pd
    df = pd.DataFrame({"lon": ALL_LONS})
    for k in ("LIB", "SB", "SBE"):
        df[f"{k}_data"], df[f"{k}_E1native"], df[f"{k}_MHW"] = f_data[k], f_model[k], MHW[k]
    df.to_csv(os.path.join(a.out, "frequencies.csv"), index=False, float_format="%.4f")
    print(df.round(3).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    long = {"LIB": "local instantaneous blocking (LIB)", "SB": "sector blocking (SB)", "SBE": "sector blocking episode (SBE)"}
    for k in ("LIB", "SB", "SBE"):
        fig, ax = plt.subplots(figsize=(9, 4.5), facecolor="#fcfcfb"); ax.set_facecolor("#fcfcfb")
        for lo, hi in ((0, 90), (330, 360)):
            ax.axvspan(lo, hi, color="#e6e5e0", alpha=0.6, lw=0, label="sector 330-90" if lo == 0 else None)
        ax.plot(ALL_LONS, f_data[k], color="#2a78d6", marker="o", ms=4, lw=2, label="data (days D+1..D+5)")
        ax.plot(ALL_LONS, f_model[k], color="#eb6834", marker="s", ms=4, lw=2, ls="--", label=f"E1 native grid, daily steps (rank {r}, {M} cells)")
        ax.set_xticks(np.arange(0, 361, 30)); ax.set_xlim(-3, 358); ax.set_ylim(0, None)
        ax.set_xlabel("longitude (deg E)"); ax.set_ylabel("frequency")
        ax.set_title(f"{long[k]} frequency, DJF, {a.nsamp} five-day windows, training winters")
        ax.grid(True, color="#e6e5e0", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False); fig.tight_layout()
        png = os.path.join(a.out, f"{k}_data_e1native.png"); fig.savefig(png, dpi=150); plt.close(fig); print("wrote", png)


if __name__ == "__main__":
    main()
