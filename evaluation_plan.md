# Evaluation plan

Metrics are in [kpis.md](kpis.md).

## Definitions

All models and scores use the definitions in `2_dim/2_dim_pipeline.ipynb`:

| Quantity | Definition used for modeling | Published definition (Tyrlis & Hoskins, MHW) |
|---|---|---|
| B | One fixed central latitude per longitude | Maximum over the central latitude and 4° either side |
| SB | A run of 3 LIB longitudes centered within 20° of a longitude | The whole run lies inside the 45° sector centered on a longitude |
| SBE | SB on all 5 days of a window | SB on that day and the previous 4 days |

`notebooks/B Value Computation.ipynb`, `Blocking Frequencies - Our B vs MHW.ipynb`, `data_audit.ipynb` and `feature_selection.ipynb` use the published B. `Blocking Frequencies - Our B vs MHW.ipynb` and `data_audit.ipynb` also use the published SB and SBE.

The choice changes frequencies substantially: SBE at 30°E on the test winters is 0.053 with one central latitude and 0.137 with the maximum. Model scores can therefore be compared with each other, and can't be compared with MHW's published curves.

## Unit of analysis

A row is one winter day at one longitude. Days within a winter are strongly correlated, and longitudes up to about 25° apart move together. The independent unit is the **winter**: splits and bootstrap resampling keep whole winters together.

## Split

Models learn from past winters and are checked on later ones. Winters are labeled by the year of their December.

| Set | Winters | Count |
|---|---|---|
| Training | Dec 1957 to Feb 1987 | 30 |
| Test | Dec 1987 to Feb 2002 | 15 |

Fitted on training winters only:
- **Eigenbasis model:** climatology, eigenbasis and AR(1) coefficients.
- **MHW-style model:** mean and noise covariance.
- **Regional model:** the same quantities as MHW-style, plus memory by region.

**Test procedure.** From each possible start day D in a test winter, a model starts from the observed state on day D and runs 5 days, 10 times. SBE is the fraction of windows with SB on all 5 days. MHW's model runs freely, so its curves measure long-run behavior and ours measure a 5-day outcome.

## Leakage review

| Risk | Status |
|---|---|
| Days from one winter on both sides of the split | Avoided: the split is by whole winters |
| Future data used in fitting | Avoided: training winters all come first |
| Fitted preprocessing | Avoided: climatology and means come from training winters |
| Start states taken from test data | Intended: models simulate from the observed state on day D |
| Central latitude overlaps test years | Present and small: the fixed curve comes from ERA-15 (1979-1993), which overlaps winters from 1987. Only its smoothed shape is used |
| Metrics chosen after modeling | Present: the KPIs and data audit were written after most of the modeling |
| Same data for choosing and scoring models | To confirm: see open item 2 |

## Stress tests (planned)

- **Shuffled time:** refit after shuffling days within training winters. Skill at SBE should collapse; if it doesn't, something is leaking.
- **Distribution shift:** train on the later winters and test on the earlier ones, since blocking varies a lot between winters.
- **Outlier winters:** drop one test winter at a time, such as 1995-96, the winter with the highest SBE.
- **Sensitivity:** vary the number of simulations per start day and the random seed.

Bootstrap intervals on RMSE differences are in place; see `notebooks/modeling_baselines.ipynb`.

## Open items

1. **Cross-validation.** The repo has one train/test split. Any cross-validation should use training winters only.
2. **Use of test winters.** Confirm whether model settings were chosen without looking at test results. The final comparison should be run once.
3. **SBE against the trivial baseline.** The regional and eigenbasis models beat the MHW-style model on SBE but aren't clearly better than the trivial baseline; see [modeling_plan.md](modeling_plan.md). Presentation claims should reflect this.
4. **Reproducibility.** The daily theta files and `B_vals/` live outside git, and `2_dim/2_dim_pipeline.ipynb` reads both. The code that makes the daily theta files, and the step that puts ERA-40 theta on a regular grid, aren't in the repo. See [data_inventory.md](data_inventory.md).
