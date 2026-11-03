# Evaluation plan

Metrics are defined in [kpis.md](kpis.md): LIB, SB and SBE frequency by longitude, plus the number of 3-, 4- and 5-day blocks.

## Unit of analysis

A row is one winter day at one longitude. Days within a winter are strongly correlated and longitudes up to about 25° apart move together, so neither is an independent sample. The independent unit is the **winter**: splits and bootstrap resampling keep whole winters together.

## Intended use and split

The goal is to learn a model from past winters and check how well it reproduces blocking statistics in later ones.

`2_dim/2_dim_pipeline.ipynb` uses a chronological split by winter, labeled by the year of its December:

| Set | Winters | Count |
|---|---|---|
| Training | Dec 1957 to Feb 1987 | 30 |
| Test | Dec 1987 to Feb 2002 | 15 |

Fitted on training winters only: the climatology, eigenbasis and AR(1) coefficients of the eigenbasis model, and the mean and noise covariance of the MHW-style model.

**Test procedure.** From each possible start day D in a test winter, each model starts from the observed state on day D and runs 5 days, 10 times. SBE is the fraction of windows with SB on all 5 days. These are not free-running simulations like MHW's, so frequencies can't be compared directly with MHW's model curves.

## Leakage review

| Risk | Status |
|---|---|
| Days from one winter on both sides of the split | Avoided: the split is by whole winters |
| Future data used in fitting | Avoided in the pipeline: training winters all come first |
| Fitted preprocessing | Avoided: climatology and mean are fitted on training winters |
| Start states taken from test data | Intended: models start from the observed state on day D, then simulate. So the test measures a 5-day outcome given a real start, not long-run behavior |
| Central latitude overlaps test years | Present and small: the fixed curve comes from ERA-15 (1979-1993), which overlaps winters from 1987. Only its smoothed shape is used |
| Same data for choosing and scoring models | To confirm: see open items |
| Fitting and scoring on the same winters | Present in `notebooks/Model Simulations and Analysis.ipynb`: the 1-D models fit and score on all 44 winters, so they are exploratory and not comparable with test results |

## Stress tests (planned)

- **Shuffled time:** refit after shuffling days within training winters. Skill at SBE should collapse; if it doesn't, something is leaking.
- **Distribution shift:** swap the split (train on the later winters, test on the earlier ones), since blocking varies a lot between winters.
- **Outlier winters:** drop one test winter at a time and see whether results hinge on it, such as 1995-96, the winter with the highest SBE.
- **Sampling noise:** bootstrap test winters for intervals on RMSE differences.
- **Sensitivity:** vary the number of simulations per start day and the random seed.

## Open items

1. **B definition.** The pipeline computes B at one central latitude. Tyrlis & Hoskins and MHW take the maximum over the central latitude and 4° either side, as `notebooks/B Value Computation.ipynb` does. Frequencies differ substantially (SBE at 30°E on the test winters: 0.053 vs 0.137). The project needs one definition, or a clear statement of why there are two.
2. **Cross-validation.** The repo has one train/test split and no cross-validation. Any cross-validation should use training winters only.
3. **Use of test winters.** Confirm whether model settings were chosen without looking at test results. The final comparison should be run once.
4. **Trivial baseline.** Not yet computed; see [kpis.md](kpis.md).
5. **Significance.** The train/test comparison has no bootstrap interval yet.
6. **Reproducibility.** `B_vals/`, `theta_native_daily/` and four scripts named in the pipeline notebook aren't in the repo: `theta_native_daily.py`, `e1_native_windows.py`, `mhw_ar1_windows.py`, `plot_test_both_models.py`. `THETA_DIR` in the notebook and a path in `e1_native_model.py` are machine-specific.
