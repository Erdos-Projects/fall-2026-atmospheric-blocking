# Modeling plan

**Status: draft.** The final choice depends on the open items in [evaluation_plan.md](evaluation_plan.md#open-items).

## Goal

Reproduce how often blocking occurs at each longitude, especially persistent episodes (SBE), with a simple statistical model. KPIs are in [kpis.md](kpis.md).

SB and SBE are defined by runs across longitudes and days, so a model has to produce B over time and longitude. Autoregressive models do this with few assumptions and stay simple enough to say something about MHW's claim that blocks sustain themselves.

## Models

We tried two models, a regional AR(1) of B and an eigenbasis AR(1) of the theta field, and compared them with two baselines. All train on 30 winters and are tested on 15. Each run starts from the observed state on a test day and simulates 5 days.

| Model | Description | Fitted quantities | Where |
|---|---|---|---|
| Trivial (baseline) | Each longitude's frequency over the training winters | 72 per measure | `notebooks/modeling_baselines.ipynb` |
| MHW-style AR(1) (baseline) | `B[t+1] - mu = phi (B[t] - mu) + eps`, with spatially correlated noise, a simplified version of the published red-noise model | 72 means and a noise covariance limited to ±15°; `phi` fixed at 0.62 | `2_dim/2_dim_pipeline.ipynb` |
| Regional AR(1) | MHW-style, with one fitted `phi` for each of 7 Tyrlis & Hoskins regions | As MHW-style, plus 7 values of `phi` | `notebooks/Model Simulations and Analysis.ipynb` |
| Eigenbasis AR(1) | An independent AR(1) for each eigenvector of the daily change in the theta field, then B computed from simulated theta | 2,589 directions, each with its own coefficient and noise size, plus the eigenvectors and a climatology | `2_dim/2_dim_pipeline.ipynb` |

## Results on the test winters

RMSE against the 15 test winters, lower is better. From `results/model_comparison.csv` and `results/regional_ar1/`.

| KPI | Trivial | MHW-style AR(1) | Regional AR(1) | Eigenbasis AR(1) |
|---|---|---|---|---|
| **SBE RMSE (primary)** | 0.0123 | 0.0142 | 0.0092 | **0.0064** |
| LIB RMSE | 0.0208 | **0.0128** | 0.0186 | 0.0132 |
| SB RMSE | 0.0321 | **0.0145** | 0.0192 | 0.0242 |
| SBE, % of longitudes inside data band | 79.2 | 73.6 | 83.3 | **86.1** |
| LIB, % of longitudes inside data band | 73.6 | **87.5** | 63.9 | 84.7 |
| SB, % of longitudes inside data band | 72.2 | **100.0** | 97.2 | 95.8 |

95% bootstrap intervals for differences in RMSE, resampling test winters. Bold marks intervals that exclude 0.

| Comparison | SBE | LIB | SB |
|---|---|---|---|
| Regional minus MHW-style | **-0.0065 to -0.0011** | **+0.0021 to +0.0060** | -0.0050 to +0.0088 |
| Regional minus eigenbasis | -0.0010 to +0.0057 | **+0.0031 to +0.0057** | -0.0084 to +0.0013 |
| Regional minus trivial | -0.0128 to +0.0082 | -0.0082 to +0.0014 | **-0.0210 to -0.0050** |
| Eigenbasis minus MHW-style | **-0.0119 to -0.0006** | -0.0021 to +0.0015 | -0.0053 to +0.0157 |
| Eigenbasis minus trivial | -0.0127 to +0.0032 | **-0.0126 to -0.0031** | -0.0195 to +0.0002 |
| MHW-style minus trivial | -0.0116 to +0.0142 | **-0.0121 to -0.0027** | **-0.0247 to -0.0042** |

- Both models beat the MHW-style baseline on SBE, the primary KPI. Neither is clearly better than the trivial baseline there.
- The regional and eigenbasis models aren't clearly different on SBE. The eigenbasis model is clearly better on LIB.
- The regional model is clearly worse than MHW-style on LIB.
- Regional `phi` ranges from 0.65 to 0.73, and every value is above MHW-style's fixed 0.62. The SBE gain can't be credited to differences between regions alone.

The trivial baseline gives only frequency curves. It can't simulate B, describe a particular winter, or test MHW's claim.

## Final model (provisional)

**Eigenbasis AR(1)**, because it has the best primary KPI and clearly improves on the MHW-style baseline there, without clearly losing on any secondary KPI.

Limits of this choice:
- **Complexity.** It has far more fitted quantities than the regional model and isn't clearly better on SBE. The repo has a single split so far, and the project guide asks for clear cross-validated improvement to justify this complexity.
- **Stationarity.** Seven of the 2,589 fitted coefficients are at or above 1 (largest 1.003). This does no harm over 5-day runs, but the model can't run freely for long, which rules out a direct comparison with MHW's free-running results.
- **Direction.** Tomorrow's B depends on longitudes to the west far more than to the east (see [data_audit.md](data_audit.md)). None of these models appears able to use this.
- **Experiment.** All models start from observed states, so a good SBE score partly reflects those starts. MHW asked whether red noise alone produces the observed episodes.

## Next steps

1. Settle the open items in [evaluation_plan.md](evaluation_plan.md#open-items), then rerun the comparison once.
2. Cross-validate within the training winters.
3. Test a model that lets information travel eastward.
