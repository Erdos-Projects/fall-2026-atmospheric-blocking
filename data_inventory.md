# Data inventory

## Source data

All source data is the ECMWF ERA-40 reanalysis, distributed by the NCAR Geoscience Data Exchange (GDEX). It covers September 1957 to August 2002, every 6 hours.

| Data | Source | Contents | Where it is | Size |
|---|---|---|---|---|
| Full 2 PVU fields | [GDEX d117004](https://gdex.ucar.edu/datasets/d117004/), DOI 10.5065/KHN2-AC11 | 540 monthly files `e4oper.an.pv.YYYYMM`. Each holds theta (`pt`), pressure, geopotential, u, v, humidity and ozone on the 2 PVU surface; reduced Gaussian N80 grid | `data/` (not in git); team Google Drive | 31 GB |
| Excerpt: theta on 2 PVU | Extracted from d117004 | `data_excerpt/theta/theta2pvu_YYYYMM.grib`: theta only, regular Gaussian grid, 320 x 160 | Git | 142 MB |
| Excerpt: wind at 300 hPa | [GDEX d117001](https://gdex.ucar.edu/datasets/d117001/), DOI 10.5065/WQZB-A142, via `fetch_era40_uv300_djf.py` | `data_excerpt/uv/uv300_YYYYMM.grib`: u and v, regular Gaussian grid, 320 x 160. `vod300_YYYYMM.grib`: the vorticity and divergence they were computed from. `template.json`: download byte offsets | Git | 426 MB |

The excerpt covers 4 winters: December 1998 to February 2002, December-February only (12 months per variable).

## Access

- **Full 2 PVU files:** download from GDEX d117004, or copy from the team's Google Drive folder. Place them directly in `data/`.
- **Wind at 300 hPa:** `python fetch_era40_uv300_djf.py --year 1957 --nwinters 45 --out-dir era40_uv300`. Needs `cdo`, which is in `environment.yml`.
- **Licence:** both GDEX data sets are CC BY 4.0. Downloading them also implies accepting ECMWF's Terms of Use.

## Derived data

| Data | Made by | Where it is |
|---|---|---|
| Daily B, Tyrlis & Hoskins method | `notebooks/B Value Computation.ipynb` | `results/ERA40_B_index_daily_1957-2002_TH08.parquet` |
| Daily B, seasonal central latitude | `notebooks/Seasonal B Value Computation.ipynb` | `results/ERA40_B_index_daily_1957-2002_seasonal-CBL.parquet` |
| Central latitude by season, 1979-1993 | `notebooks/Seasonal B Value Computation.ipynb` | `results/ERA40_CBL_by_season_1979-1993.csv` |
| Drift and noise tables | `notebooks/Drift and Noise Analysis.ipynb` | `results/ERA40_winter_B_analysis.sqlite` |
| LIB, SB, SBE per winter, 45 winters | `blocking_episodes.py`, from `B_YYYYMM.nc` | `episodes/` |
| Monthly B, `B_YYYYMM.nc` | `blocking_index.py` | `B_vals/`: not in git |
| Daily-mean theta, native grid | `theta_native_daily.py`: not in repo yet | `theta_native_daily/`, 794 MB: not in git |
| Model fits, simulations, test frequencies | `2_dim/2_dim_pipeline.ipynb` | `2_dim/notebook_output/` |

## Limitations

- **One reanalysis, few winters.** There are about 45 winters, and consecutive days are strongly correlated, so the effective sample is far smaller than the number of days.
- **Changing observations.** ERA-40 combines an observing system that changed over the period, so the record is not perfectly uniform in time.
- **Superseded product.** GDEX marks these data sets as historical. We use ERA-40 to compare directly with Tyrlis & Hoskins and MHW.
- **Wind level.** Tyrlis & Hoskins define the central latitude from 300 hPa wind. We have 300 hPa wind only for the excerpt, so `repro_MHW.ipynb` and the seasonal notebook use wind on the 2 PVU surface.
- **Digitized reference values.** The published central latitude and MHW's frequency curves were read from figures, so they are approximate.

## Still to document

- How `data_excerpt/theta/` was extracted from d117004 (no script in the repo).
- Download dates for the full files.
