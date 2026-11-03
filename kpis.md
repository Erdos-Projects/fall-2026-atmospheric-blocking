# Key performance indicators

These were written down after most of the modeling, based on the targets the team used during the project. They were not fixed in advance.

## What is scored

For each model, simulated B gives a frequency curve for each measure: the fraction of winter days with LIB, SB or SBE at each of the 72 longitudes. Each curve is compared with the same curve from ERA-40 on winters held out from fitting.

Write `f_model(x)` and `f_data(x)` for the model and data frequency at longitude `x`.

## Primary KPI

**SBE RMSE. Lower is better.**

    RMSE = sqrt( (1/72) * sum over x of (f_model(x) - f_data(x))^2 )

SBE is primary because persistent episodes are what MHW's red-noise model fails to capture, and they are the most consequential kind of blocking.

## Secondary KPIs

| KPI | Definition | Better | Why it matters |
|---|---|---|---|
| LIB RMSE | As above, for LIB | Lower | A model that improves SBE by distorting single-day behavior isn't really better |
| SB RMSE | As above, for SB | Lower | Checks the structure across longitudes |
| Share of longitudes inside the data band | % of longitudes with `abs(f_model - f_data) <= 2 * se`, where `se` = standard deviation across test winters / sqrt(number of test winters) | Higher | Measures misses relative to sampling noise |
| Number of 3-, 4- and 5-day blocks | Count of SB spells of each length | Closer to data | Checks duration directly (planned, not yet computed) |

RMSE weights all longitudes equally, so a small miss where blocking is rare counts as much as the same miss over Europe. The tables should also report the 330-90°E sector, as the pipeline already does.

## Baselines

| Baseline | Definition | Status |
|---|---|---|
| Trivial | Predict each longitude's frequency from the training winters, with no model | Planned, not yet computed |
| MHW-style AR(1) | `B[t+1] - mu = phi * (B[t] - mu) + eps`, with `eps` drawn from a Gaussian with tapered spatial covariance | In `2_dim/2_dim_pipeline.ipynb`. A simplified version of MHW's model, not an exact reproduction |

## Uncertainty

Model differences should be compared by resampling whole winters with replacement (a block bootstrap), and improvements called significant only if the 95% interval for the RMSE difference excludes 0. `notebooks/Model Simulations and Analysis.ipynb` does this for the 1-D models; the train/test comparison doesn't have it yet.

## Current values (provisional)

From `2_dim/notebook_output/frequencies_test.csv`, 15 test winters. These use the single-central-latitude B, and may change; see "Open items" in [evaluation_plan.md](evaluation_plan.md).

| KPI | Eigenbasis AR(1) | MHW-style AR(1) |
|---|---|---|
| SBE RMSE | 0.0064 | 0.0142 |
| LIB RMSE | 0.0132 | 0.0128 |
| SB RMSE | 0.0242 | 0.0145 |
| SBE, % of longitudes inside band | 86.1 | 73.6 |
| LIB, % of longitudes inside band | 84.7 | 87.5 |
| SB, % of longitudes inside band | 95.8 | 100.0 |
