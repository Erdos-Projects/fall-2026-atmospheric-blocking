# Modeling Atmospheric Blocking: executive summary

Matt King, Jun Park, Richard Carini. Erdős Institute, Fall 2026.

**Status: provisional.** Cross-validation is still to come, and the numbers below may change.

## Problem

Atmospheric blocking is a persistent high-pressure system that diverts the usual west-to-east flow, and it is linked to cold spells and heat waves. A 2009 study (Masato, Hoskins & Woollings, "MHW") found that a simple red-noise model reproduces how often short-lived blocking occurs but underestimates persistent episodes, and proposed that blocks sustain themselves.

We measured how much of that gap a better statistical model can close. The answer bears on whether persistent blocking needs a physical explanation beyond chance, which matters to anyone interpreting blocking in forecasts or climate models.

## Data

ERA-40 reanalysis from ECMWF, distributed by NCAR GDEX: 45 winters (December 1957 to February 2002), every 6 hours, on the 2 PVU surface near the tropopause. All 540 monthly files are complete, with no missing values. We model December-February in the Northern Hemisphere, at 72 longitudes.

The blocking index B compares temperature north and south of a central latitude; B > 0 means the normal pattern has reversed. Three measures build on it:
- **LIB:** B > 0 at a longitude on a day.
- **SB:** a run of 3 neighboring LIB longitudes centered within 20° of a longitude.
- **SBE:** SB on all 5 days of a window.

## Approach and KPIs

Models train on the first 30 winters and are tested on the last 15, split by whole winters. From each test day, a model starts at the observed state and simulates 5 days.

The primary KPI is the RMSE of SBE frequency across longitudes, with LIB and SB RMSE as secondary KPIs. Bootstrap intervals come from resampling the test winters. We tried two models and compared them with two baselines:
- **Trivial (baseline):** each longitude's frequency from the training winters, with no model.
- **MHW-style AR(1) (baseline):** red noise on B with correlated noise across longitudes, a simplified version of MHW's model.
- **Regional AR(1):** the MHW-style model with one fitted memory value for each of 7 published regions.
- **Eigenbasis AR(1):** an independent AR(1) for each eigenvector of daily changes in the temperature field, with B computed from the simulated field.

## Results

RMSE on the 15 test winters (lower is better):

| Measure | Trivial | MHW-style AR(1) | Regional AR(1) | Eigenbasis AR(1) |
|---|---|---|---|---|
| **SBE (primary)** | 0.0123 | 0.0142 | 0.0092 | **0.0064** |
| LIB | 0.0208 | **0.0128** | 0.0186 | 0.0132 |
| SB | 0.0321 | **0.0145** | 0.0192 | 0.0242 |

- Both models beat the MHW-style baseline on SBE: the 95% intervals for the differences are -0.0065 to -0.0011 (regional) and -0.0119 to -0.0006 (eigenbasis).
- The regional and eigenbasis models aren't clearly different on SBE (-0.0010 to +0.0057). The eigenbasis model is clearly better on LIB.
- Regional memory ranges from 0.65 to 0.73, and every value is above MHW-style's fixed 0.62, so the SBE gain can't be credited to differences between regions alone.
- Neither model is clearly better than the trivial baseline on SBE. Fifteen winters leave wide intervals, because SBE varies greatly from winter to winter.
- Tomorrow's B depends far more on longitudes to the west than to the east, and none of these models appears able to use this.

## Final model

**Eigenbasis AR(1), provisionally.** It has the best primary KPI and clearly improves on the MHW-style model, without clearly losing on a secondary KPI. It is also far more complex than the regional model, with about 2,600 fitted directions, and isn't clearly better on SBE.

## Limitations

- **B definition.** We compute B at one central latitude. The published index takes the maximum over three, which gives much higher frequencies, so these numbers can be compared with each other and can't be compared with MHW's.
- **Experiment.** Our models run for 5 days from observed states and MHW's runs freely, so our results don't directly test their claim that blocks sustain themselves.
- **One split, 15 test winters.** Cross-validation is still to do.
- **Short runs.** Seven of the eigenbasis model's coefficients are at or above 1, so it can't run freely for long.
- **Data.** ERA-40 is one superseded reanalysis, built from an observing system that changed over time.

## Next steps

1. Cross-validate within the training winters.
2. Test a model that lets information travel eastward.

## More detail

[README](../README.md), [data inventory](../data_inventory.md), [KPIs](../kpis.md), [evaluation plan](../evaluation_plan.md), [data audit](../data_audit.md), [modeling plan](../modeling_plan.md). Results: `results/model_comparison.csv`, from `notebooks/modeling_baselines.ipynb`.
