#!/usr/bin/env python3
"""
AR(1) model of the blocking index B over the Atlantic/European sector
(lon 330, 335, ..., 355, 0, 5, ..., 90; 25 longitudes) and a comparison of the
LIB / SB / SBE frequencies it produces against the data and the MHW values.

Model:  B_{t+1} - mu = PHI * (B_t - mu) + eps,   eps ~ N(0, M)
  mu   = per-longitude mean of B over every day in the dataset;
  M_xy = cov(r_x, r_y) over all within-winter consecutive-day pairs,
         r = (B_{t+1} - mu) - PHI * (B_t - mu),  for |x - y| <= 15 deg (circular distance);
  M_xy = 0 otherwise.
The tapered M is not positive semidefinite, so the nearest PSD matrix in the
Frobenius norm (negative eigenvalues clipped to 0) is used for sampling.

Experiment: 1000 start days D (uniform, without replacement) with D..D+5 inside
one winter.  B_0 = data on day D; simulate B_1..B_5.  LIB/SB are averaged over
the 5 simulated days and the 1000 samples; SBE = SB on all 5 days, averaged over
the 1000 samples.  The same statistics are computed from the data on D+1..D+5.

LIB/SB/SBE follow blocking_episodes.py, but restricted to the 25-longitude
sector: longitudes outside the sector contribute no LIB (no wrap-around), so
SB/SBE within 20 deg of the sector edges (330-350 and 70-90) are slightly low.
For reference the data SB/SBE are also computed on the full 72-longitude
periodic grid ("data_full" columns in the CSV).

Usage:  python ar1_sector_model.py [--b-dir blocking/B_vals] [--out blocking/ar1_output]
"""
import argparse
import glob
import os
import re

import numpy as np
import xarray as xr

PHI = 0.62
TAPER_DEG = 15
NSAMP = 1000
NDAYS = 5
SB_HALFWIDTH = 4
SEED = 0

ALL_LONS = np.arange(0, 360, 5)
SECTOR = np.r_[np.arange(330, 360, 5), np.arange(0, 95, 5)]          # 25 longitudes
SECTOR_IDX = np.array([int(np.flatnonzero(ALL_LONS == l)[0]) for l in SECTOR])

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


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_winters(b_dir):
    """Return list of (december_year, B[ndays, 72]) with days consecutive within each winter."""
    groups = {}
    for f in glob.glob(os.path.join(b_dir, "B_*.nc")):
        m = re.search(r"B_(\d{4})(\d{2})\.nc$", os.path.basename(f))
        if not m:
            continue
        y, mo = int(m.group(1)), int(m.group(2))
        if mo == 12:
            groups.setdefault(y, []).append(f)
        elif mo in (1, 2):
            groups.setdefault(y - 1, []).append(f)
    out = []
    for y, paths in sorted(groups.items()):
        ds = xr.concat([xr.open_dataset(p)[["B"]].load() for p in sorted(paths)], dim="time").sortby("time")
        t = ds.time.values
        if np.any(np.diff(t) != np.timedelta64(1, "D")):
            raise ValueError(f"winter {y}: days are not consecutive")
        assert np.array_equal(ds.lon.values, ALL_LONS)
        out.append((y, ds.B.values.astype(float)))
    return out


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------
def fit_mean(winters):
    """Per-longitude mean of B over every day in the dataset (25,)."""
    return np.concatenate([B[:, SECTOR_IDX] for _, B in winters]).mean(axis=0)


def fit_noise_cov(winters, mu):
    """Tapered residual covariance M (25 x 25) and its PSD projection."""
    res = np.concatenate([(B[1:, SECTOR_IDX] - mu) - PHI * (B[:-1, SECTOR_IDX] - mu) for _, B in winters])
    C = np.cov(res.T)                                       # ddof = 1
    d = np.abs(SECTOR[:, None] - SECTOR[None, :]) % 360
    d = np.minimum(d, 360 - d)
    M = np.where(d <= TAPER_DEG, C, 0.0)
    w, V = np.linalg.eigh(M)
    M_psd = (V * np.clip(w, 0, None)) @ V.T
    L = V * np.sqrt(np.clip(w, 0, None))                   # eps = L @ z, z ~ N(0, I)
    info = dict(n_pairs=res.shape[0], resid_mean=res.mean(0), eig=w,
                frob_change=np.linalg.norm(M_psd - M) / np.linalg.norm(M))
    return M, M_psd, L, info


def simulate(B0, mu, L, rng):
    """B0: (nsamp, 25).  Returns (nsamp, NDAYS, 25) with days 1..NDAYS."""
    nsamp = B0.shape[0]
    z = rng.standard_normal((nsamp, NDAYS, L.shape[1]))
    eps = z @ L.T
    out = np.empty((nsamp, NDAYS, B0.shape[1]))
    prev = B0
    for t in range(NDAYS):
        prev = mu + PHI * (prev - mu) + eps[:, t]
        out[:, t] = prev
    return out


# ---------------------------------------------------------------------------
# blocking indices on a (nsamp, ndays, nlon) array
# ---------------------------------------------------------------------------
def _shift(a, k, periodic):
    """Shift along the last axis so that result[..., x] = a[..., x - k] (zero padded unless periodic)."""
    if periodic:
        return np.roll(a, k, axis=-1)
    out = np.zeros_like(a)
    if k > 0:
        out[..., k:] = a[..., :-k]
    elif k < 0:
        out[..., :k] = a[..., -k:]
    else:
        out[...] = a
    return out


def lib_sb_sbe(B, periodic):
    """LIB, SB per (sample, day, lon); SBE per (sample, lon) = SB on every day of the window."""
    lib = B > 0
    run3 = lib & _shift(lib, 1, periodic) & _shift(lib, -1, periodic)
    sb = np.zeros_like(lib)
    for k in range(-SB_HALFWIDTH, SB_HALFWIDTH + 1):
        sb |= _shift(run3, k, periodic)
    sbe = sb.all(axis=1)
    return lib, sb, sbe


def frequencies(B, periodic=False):
    lib, sb, sbe = lib_sb_sbe(B, periodic)
    return {"LIB": lib.mean(axis=(0, 1)), "SB": sb.mean(axis=(0, 1)), "SBE": sbe.mean(axis=0)}


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--b-dir", default=os.path.join(os.path.dirname(__file__), "B_vals"))
    p.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "ar1_output"))
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--no-mhw", action="store_true",
                   help="also write figures with only the data and model curves (suffix _no_mhw)")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    winters = load_winters(a.b_dir)
    print(f"{len(winters)} winters, {sum(B.shape[0] for _, B in winters)} days")

    mu = fit_mean(winters)
    print("mu by longitude:", np.round(mu, 2))
    M, M_psd, L, info = fit_noise_cov(winters, mu)
    print(f"residual pairs: {info['n_pairs']}")
    print(f"tapered M eigenvalues: min {info['eig'].min():.2f}, max {info['eig'].max():.2f}, "
          f"{(info['eig'] < 0).sum()} negative; PSD projection changes M by "
          f"{100 * info['frob_change']:.1f}% (Frobenius)")
    print("residual mean by longitude (model uses 0):", np.round(info["resid_mean"], 2))

    # start days: (winter index, day index) with D..D+5 inside the winter
    starts = [(w, i) for w, (_, B) in enumerate(winters) for i in range(B.shape[0] - NDAYS)]
    pick = rng.choice(len(starts), size=NSAMP, replace=False)
    chosen = [starts[k] for k in pick]
    print(f"{len(starts)} admissible start days, {NSAMP} drawn without replacement")

    B0 = np.stack([winters[w][1][i, SECTOR_IDX] for w, i in chosen])                       # (N, 25)
    B_data = np.stack([winters[w][1][i + 1:i + 1 + NDAYS][:, SECTOR_IDX] for w, i in chosen])  # (N, 5, 25)
    B_full = np.stack([winters[w][1][i + 1:i + 1 + NDAYS] for w, i in chosen])             # (N, 5, 72)
    B_model = simulate(B0, mu, L, rng)

    f_model = frequencies(B_model)
    f_data = frequencies(B_data)
    f_full = {k: v[SECTOR_IDX] for k, v in frequencies(B_full, periodic=True).items()}
    f_mhw = {k: v[SECTOR_IDX] for k, v in MHW.items()}

    # ---- table
    import pandas as pd
    df = pd.DataFrame({"lon": SECTOR})
    for k in ("LIB", "SB", "SBE"):
        df[f"{k}_data"] = f_data[k]
        df[f"{k}_model"] = f_model[k]
        df[f"{k}_MHW"] = f_mhw[k]
        df[f"{k}_data_full72"] = f_full[k]
    csv = os.path.join(a.out, "frequencies.csv")
    df.to_csv(csv, index=False, float_format="%.4f")
    np.savetxt(os.path.join(a.out, "M_tapered.csv"), M, delimiter=",", fmt="%.4f")
    np.savetxt(os.path.join(a.out, "M_psd.csv"), M_psd, delimiter=",", fmt="%.4f")
    np.savetxt(os.path.join(a.out, "mu.csv"), np.c_[SECTOR, mu], delimiter=",", fmt="%.4f", header="lon,mu")
    with open(os.path.join(a.out, "start_days.txt"), "w") as fh:
        for w, i in chosen:
            fh.write(f"{winters[w][0]} {i}\n")
    print("\nday-by-day sector means (B in K; LIB = fraction of longitudes with B > 0)")
    print("  day   B data   B model   LIB data  LIB model")
    for t in range(NDAYS):
        print(f"  {t + 1}    {B_data[:, t].mean():7.2f}   {B_model[:, t].mean():7.2f}     "
              f"{(B_data[:, t] > 0).mean():.3f}     {(B_model[:, t] > 0).mean():.3f}")
    print("\nsector-mean frequencies")
    for k in ("LIB", "SB", "SBE"):
        print(f"  {k:3s}  data {f_data[k].mean():.3f}  model {f_model[k].mean():.3f}  "
              f"MHW {f_mhw[k].mean():.3f}  data(full72) {f_full[k].mean():.3f}")
    print(df.round(3).to_string(index=False))

    # ---- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.where(SECTOR >= 180, SECTOR - 360, SECTOR).astype(float)       # -30 ... 90
    ticks = np.arange(-30, 91, 15)
    labels = [f"{int(t % 360)}" for t in ticks]
    series = [("data (days D+1..D+5)", f_data, "#2a78d6", "o", "-"),
              ("AR(1) model about mu", f_model, "#eb6834", "s", "--"),
              ("MHW", f_mhw, "#1baf7a", "^", "-.")]
    long = {"LIB": "local instantaneous blocking (LIB)", "SB": "sector blocking (SB)",
            "SBE": "sector blocking episode (SBE)"}
    variants = [("", series)] + ([("_no_mhw", series[:2])] if a.no_mhw else [])
    for k, (suffix, ser) in [(k, v) for k in ("LIB", "SB", "SBE") for v in variants]:
        fig, ax = plt.subplots(figsize=(8, 4.5), facecolor="#fcfcfb")
        ax.set_facecolor("#fcfcfb")
        for name, f, col, mk, ls in ser:
            ax.plot(x, f[k], color=col, marker=mk, ms=5, lw=2, ls=ls, label=name)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_xlim(-32, 92)
        ax.set_ylim(0, None)
        ax.set_xlabel("longitude (deg E)")
        ax.set_ylabel("frequency")
        ax.set_title(f"{long[k]} frequency, DJF, {NSAMP} five-day windows")
        ax.grid(True, color="#e6e5e0", lw=0.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        png = os.path.join(a.out, f"{k}_data_model_mhw{suffix}.png")
        fig.savefig(png, dpi=150)
        plt.close(fig)
        print("wrote", png)
    print("wrote", csv)


if __name__ == "__main__":
    main()
