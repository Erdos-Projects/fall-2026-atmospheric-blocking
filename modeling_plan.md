# Modeling plan

**Status: draft.** The final choice depends on the open items in [evaluation_plan.md](evaluation_plan.md), especially which definition of B the project uses. This document was written at the end of the project, looking back at what was done.

## Goal

Reproduce how often blocking occurs at each longitude, especially persistent episodes (SBE), with a simple statistical model. KPIs are in [kpis.md](kpis.md).

We aren't predicting a label for each row, so standard classifiers and tree ensembles don't fit the question. The output has to be simulated B over time and longitude, because SB and SBE are defined by runs across longitudes and days. A model also has to be simple enough to say something about MHW's claim that blocks sustain themselves. This is why we stayed with autoregressive models.

## Models tried

### Train/test comparison

Trained on 30 winters and tested on 15, in `2_dim/2_dim_pipeline.ipynb`. Each run starts from the observed state on a test day and simulates 5 days. Scores are computed in `notebooks/modeling_baselines.ipynb`.

| Model | Description | Fitted quantities | Why tried |
|---|---|---|---|
| Trivial | Each longitude's frequency over the training winters | 72 per measure | Reference point, using no model and no test data |
| MHW-style AR(1) | `B[t+1] - mu = phi (B[t] - mu) + eps`, with spatially correlated noise | 72 means and a noise covariance limited to ±15°; `phi` fixed at 0.62 | Simplified version of the published red-noise model |
| Eigenbasis AR(1) | An independent AR(1) for each eigenvector of the daily change in the theta field, then B computed from simulated theta | 2,589 directions, each with its own coefficient and noise size, plus the eigenvectors and a climatology | Models theta rather than B, so spatial structure comes from the data |

### Exploratory 1-D models

Models of B in `notebooks/Model Simulations and Analysis.ipynb`, using B with the Tyrlis & Hoskins maximum and a winter central latitude. These fit and score on all 44 winters and run freely, so they **can't be compared with the table above**. They were used to learn what matters.

| Model | SBE RMSE, same winters | Finding |
|---|---|---|
| AR(1), one set of parameters for the hemisphere | 0.0183 | Over Eurasia it predicts too much LIB, because it gives every longitude the same spread |
| AR(1) per Tyrlis & Hoskins region (7) | 0.0101 | Fitting each region separately captures most of the geography |
| AR(1) per data-driven region (7) | 0.0089 | Cuts chosen from memory and spread alone. A small but clear improvement over the Tyrlis & Hoskins regions (bootstrap) |

## Results on the test winters

RMSE against the 15 test winters, lower is better. From `results/model_comparison.csv`; these use the pipeline's definitions of B and SB (see open item 1).

| KPI | Trivial | MHW-style AR(1) | Eigenbasis AR(1) |
|---|---|---|---|
| **SBE RMSE (primary)** | 0.0123 | 0.0142 | **0.0064** |
| LIB RMSE | 0.0208 | **0.0128** | 0.0132 |
| SB RMSE | 0.0321 | **0.0145** | 0.0242 |
| SBE, % of longitudes inside data band | 79.2 | 73.6 | **86.1** |
| LIB, % of longitudes inside data band | 73.6 | **87.5** | 84.7 |
| SB, % of longitudes inside data band | 72.2 | **100.0** | 95.8 |

95% bootstrap intervals for differences in RMSE, resampling test winters (`results/model_comparison_rmse_differences.csv`):

| Comparison | SBE | LIB | SB |
|---|---|---|---|
| Eigenbasis minus MHW-style | **-0.0119 to -0.0006** | -0.0021 to +0.0015 | -0.0053 to +0.0157 |
| Eigenbasis minus trivial | -0.0127 to +0.0032 | **-0.0126 to -0.0031** | -0.0195 to +0.0002 |
| MHW-style minus trivial | -0.0116 to +0.0142 | **-0.0121 to -0.0027** | **-0.0247 to -0.0042** |

Bold marks intervals that exclude 0.

**What this supports:**
- The eigenbasis model beats the MHW-style model on the primary KPI.
- Both models clearly beat the trivial baseline on LIB, and the MHW-style model does on SB.
- The two models aren't clearly different on LIB or SB. The higher SB RMSE for the eigenbasis model is within sampling noise.

**What it doesn't support:**
- The eigenbasis model isn't clearly better than the trivial baseline on SBE, though its RMSE is about half. With 15 winters the interval is wide.
- The MHW-style model does no better than the trivial baseline on SBE.

The trivial baseline's limits matter too. It gives only frequency curves: it can't simulate B or say anything about a particular winter. It also can't test MHW's claim.

## Final model (provisional)

**Eigenbasis AR(1)**, because it has the best primary KPI and clearly improves on the MHW-style model there, without clearly losing on any secondary KPI.

Reasons for caution:
- **Complexity.** It has far more fitted quantities than the alternatives. The guide asks for clear cross-validated improvement to justify this, and the repo has a single split so far.
- **Stationarity.** Some fitted coefficients reach 1 (largest 1.003). That does no harm over 5-day runs, but the model can't be run freely for long, so it can't be directly compared with MHW's free-running results.
- **Direction.** Tomorrow's B depends on longitudes to the west far more than to the east (see [data_audit.md](data_audit.md)). Neither AR(1) appears able to use this.
- **Different experiment.** Both models start from observed states, so a good SBE score partly reflects those starts. This differs from MHW's question of whether red noise alone produces the observed episodes.

## What didn't work

| Approach | Why it was dropped |
|---|---|
| AR(2) from lags 1 and 2 | Second coefficient came out negative, giving less memory than AR(1) |
| Two-timescale AR(1) | Fit autocorrelation well, without clear SBE improvement, at the cost of complexity. Compared before the SBE counting was corrected |
| Three data-driven regions | Too coarse, lumping Europe with Asia |
| Per-longitude parameters | Rejected as too close to fitting the data directly |
| Quantile mapping to observed distributions | Rejected for the same reason: it forces LIB to match by construction |
| Counting every day of a 5-day spell as SBE | Doesn't match MHW's curves. Counting from day 5 onward does (RMSE 0.003) |
| Whole grid rows in the B latitude bands | Caused a 0.6 K bias over the Pacific. Interpolated bands fixed it (`notebooks/B Value Computation.ipynb`) |

## Built but not yet run

`notebooks/Model Simulations and Analysis.ipynb` includes two models with no saved results, so no conclusions can be drawn:
- **Model A:** kick size depends on B, aimed at the overshoot near 210-260°E, where B is strongly skewed.
- **Model B:** memory is higher when B > 0, testing "self-sustaining" blocks directly.

## Next steps

1. Settle the open items in [evaluation_plan.md](evaluation_plan.md), then rerun the comparison once.
2. Cross-validate within the training winters.
3. Test a model that lets information travel eastward.
4. Run Models A and B, if they are to be reported.
