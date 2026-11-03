# Modeling Atmospheric Blocking

Matt King, Jun Park, Richard Carini. Erdős Institute data science project, Fall 2026.

## The problem

Atmospheric blocking is a sustained, mostly stationary high-pressure system that disrupts the usual west-to-east flow. Blocks are linked to prolonged or extreme regional weather, such as winter cold spells and summer heat waves.

Blocking is measured with an index B, the difference in potential temperature (theta) between the poleward and equatorward sides of a central latitude, on the 2 PVU surface near the tropopause (Tyrlis & Hoskins 2008). B > 0 means the usual north-south pattern has reversed.

**Question:** how well can a simple statistical model reproduce how often blocking occurs at each longitude, especially persistent episodes? Masato, Hoskins & Woollings (2009, "MHW") found that a red-noise AR(1) model of B matches short-lived blocking but underestimates persistent episodes, and attributed the difference to blocks sustaining themselves. We test how much of that gap a better statistical model can close.

**Why it matters:** if red noise can't explain persistent blocking, the extra persistence needs a physical explanation and there is room for forecasts to beat simple statistics. If a statistical model can explain it, that claim needs revisiting. The audience is atmospheric scientists and anyone evaluating blocking in forecast or climate models.

**Unit of analysis:** one winter day at one longitude (72 longitudes, 0, 5, ..., 355°E). Whole winters are the independent units for splitting and uncertainty, because days within a winter are correlated.

**Scope:** Northern Hemisphere, December-February, ERA-40 reanalysis, 1957-2002.

**Not addressed:** forecasting individual blocking events; causes of or trends in blocking under climate change; seasons other than winter in the final comparison; the Southern Hemisphere; full physics-based simulation.

**Blocking measures** (the modeling targets; see [kpis.md](kpis.md)):

| Measure | Meaning |
|---|---|
| LIB (local instantaneous blocking) | B > 0 at a longitude on a day |
| SB (sector blocking) | LIB at 3 consecutive longitudes inside the 45° sector centered on a longitude |
| SBE (sector blocking episode) | SB at a longitude on that day and the previous 4 days |

## Project documents

- [data_inventory.md](data_inventory.md): data sources, access, licensing, limitations
- [kpis.md](kpis.md): metrics, formulas, baselines
- [data_audit.md](data_audit.md): data audit, feature decisions, what was tried and discarded
- [modeling_plan.md](modeling_plan.md): models tried, test results, provisional final choice, what didn't work
- [schema.json](schema.json): inputs, derived features and targets
- [evaluation_plan.md](evaluation_plan.md): train/test split, leakage checks, stress tests, open items

## Setup

```bash
conda env create -f environment.yml
conda activate erdos_ds_era40_env
```

In Jupyter or VS Code, select this environment as the notebook's kernel.

The full data set (31 GB) is not in git. Put the monthly `e4oper.an.pv.YYYYMM` files in `data/`; [data_inventory.md](data_inventory.md) says where to get them. `data_excerpt/` holds 4 winters for quick tests.

## Repository layout

| Path | Contents |
|---|---|
| `data_excerpt/` | 4 winters of theta and wind, tracked for testing |
| `data/` | Full ERA-40 files (ignored by git) |
| `notebooks/` | Analysis notebooks, listed below |
| `2_dim/` | Eigenbasis AR(1) and MHW-style pipeline with its helper modules and outputs |
| `results/` | Saved B values, seasonal central latitudes, analysis database |
| `results/eda/` | EDA animations |
| `episodes/` | LIB, SB and SBE per winter from `src/blocking_episodes.py` |
| `src/` | Scripts, listed below |

### Notebooks

Unless noted, they read the full data or the saved B files in `results/`.

| Notebook | Purpose |
|---|---|
| `eda.ipynb` | Theta and wind exploration, PCA on theta, blocking animations (uses `data_excerpt/`) |
| `data_audit.ipynb` | Completeness of the full files; distributions, outliers, correlations and lag structure of B (training winters) |
| `feature_selection.ipynb` | Tests what each candidate input adds to predicting tomorrow's B; writes `results/feature_selection.csv` |
| `B Value Computation.ipynb` | B from theta following Tyrlis & Hoskins; writes `results/ERA40_B_index_daily_1957-2002_TH08.parquet`. Run this first |
| `Blocking Frequencies - Our B vs MHW.ipynb` | Checks LIB, SB and SBE from that B against values digitized from MHW Fig. 2 |
| `repro_MHW.ipynb` | Reproduces the central blocking latitude and tests it against the published curve |
| `Seasonal B Value Computation.ipynb` | Central latitude per season and B computed with it |
| `Drift and Noise Analysis.ipynb` | How the day-to-day change in B depends on B; writes `results/ERA40_winter_B_analysis.sqlite` |
| `Model Simulations and Analysis.ipynb` | 1-D models of B: MHW-style AR(1), AR(2), regional models, bootstrap, state-dependent models |
| `modeling_baselines.ipynb` | Scores both models against a trivial baseline, with bootstrap intervals; writes `results/model_comparison.csv` |
| `2_dim/2_dim_pipeline.ipynb` | Trains on 30 winters and tests on 15: eigenbasis AR(1) vs MHW-style AR(1). Needs `B_vals/` and daily theta files that are not yet in the repo |

`pipeline_demo.ipynb` is a placeholder; the modeling pipeline is in `2_dim/2_dim_pipeline.ipynb`.

### Scripts in `src/`

| Script | Purpose |
|---|---|
| `fetch_era40_uv300_djf.py` | Downloads ERA-40 wind at 300 hPa from GDEX |
| `coarsen_theta_long.py` | Averages theta from 320 to 72 longitudes |
| `blocking_index.py` | B at one central latitude, per month (`B_YYYYMM.nc`) |
| `blocking_episodes.py` | LIB, SB and SBE per winter from those files |
| `plot_era40.py` | Maps of theta and wind; also the GRIB reader used by the notebooks |
| `plot_blocking.py` | Plots of B |
| `plot_sbe_frequency.py` | Frequencies by longitude; holds the digitized MHW curves |
| `check_correlations.py` | Correlations between frequency curves |

Notebooks import these with `sys.path.append("../src")`.

## Command-line tools

Run from the top-level folder.

Plot theta, with or without wind:

```bash
python src/plot_era40.py --theta data_excerpt/theta/theta2pvu_199812.grib --time 1998-12-15T12
python src/plot_era40.py --theta data_excerpt/theta/theta2pvu_199812.grib --uv data_excerpt/uv/uv300_199812.grib --time 1998-12-15T12
```

Plot B from `B_YYYYMM.nc` files made by `src/blocking_index.py`:

```bash
python src/plot_blocking.py B_199812.nc                                        # time x longitude heatmap
python src/plot_blocking.py B_199812.nc B_199901.nc B_199902.nc --out w9899.png  # a whole winter
python src/plot_blocking.py B_199812.nc --day 1998-12-15                       # longitude x central-latitude map for one day
```

Plot and compare blocking frequencies:

```bash
python src/plot_sbe_frequency.py episodes --all               # LIB, SB and SBE by longitude
python src/plot_sbe_frequency.py episodes --all --mhw         # compare with MHW
python src/plot_sbe_frequency.py episodes --all --mhw --cbl   # add the central latitude
python src/check_correlations.py episodes SBE MHW-SBE
python src/check_correlations.py episodes SBE CBL
```

## Status

Final stage, with results under review. The B definition and some evaluation details are still being reconciled between the pipelines; see "Open items" in [evaluation_plan.md](evaluation_plan.md). Treat numbers as provisional until those are closed.

## References

- Masato, G., B. J. Hoskins, and T. J. Woollings, 2009: Can the frequency of blocking be described by a red noise process? *J. Atmos. Sci.*, **66**, 2143-2149.
- Tyrlis, E., and B. J. Hoskins, 2008: Aspects of a Northern Hemisphere atmospheric blocking climatology. *J. Atmos. Sci.*, **65**, 1638-1652.
- Pelly, J. L., and B. J. Hoskins, 2003: A new perspective on blocking. *J. Atmos. Sci.*, **60**, 743-755.
- Woollings, T., et al., 2018: Blocking and its response to climate change. *Curr. Clim. Change Rep.*, **4**, 287-300.
