# Modeling Atmospheric Blocking: executive summary

Matt King, Jun Park, Richard Carini. Erdős Institute, Fall 2026.

**Status: provisional.** The team is settling how blocking is defined in the final comparison (see Limitations), and the numbers below may change.

## Problem

Atmospheric blocking is a persistent high-pressure system that diverts the usual west-to-east flow, and it is linked to cold spells and heat waves. A 2009 study (Masato, Hoskins & Woollings, "MHW") found that a simple red-noise model reproduces how often short-lived blocking occurs, but underestimates persistent episodes, and proposed that blocks sustain themselves.

We asked how much of that gap a better statistical model can close. This matters to anyone interpreting blocking in forecasts or climate models: it bears on whether persistent blocking needs a physical explanation beyond chance.

## Data

ERA-40 reanalysis from ECMWF, distributed by NCAR GDEX: 45 winters (December 1957 to February 2002), every 6 hours, on the 2 PVU surface near the tropopause. All 540 monthly files are complete, with no missing values. We model December-February in the Northern Hemisphere, at 72 longitudes.

The blocking index B compares temperature north and south of a central latitude; B > 0 means the normal pattern has reversed. Three measures build on it:
- **LIB:** B > 0 at a longitude on a day.
- **SB:** LIB across 3 neighboring longitudes within a 45° sector.
- **SBE:** SB lasting at least 5 days.

## Approach and KPIs

Models train on the first 30 winters and are tested on the last 15, split by whole winters. From each test day, a model starts at the observed state and simulates 5 days.

The primary KPI is the RMSE of SBE frequency across longitudes, with LIB and SB RMSE as secondary KPIs. Bootstrap intervals come from resampling the test winters. Three models are compared:
- **Trivial:** each longitude's frequency from the training winters, with no model.
- **MHW-style AR(1):** red noise on B with correlated noise across longitudes, a simplified version of MHW's model.
- **Eigenbasis AR(1), our model:** an independent AR(1) for each eigenvector of daily changes in the temperature field, with B computed from the simulated field.

## Results

RMSE on the 15 test winters (lower is better):

| Measure | Trivial | MHW-style AR(1) | Eigenbasis AR(1) |
|---|---|---|---|
| **SBE (primary)** | 0.0123 | 0.0142 | **0.0064** |
| LIB | 0.0208 | **0.0128** | 0.0132 |
| SB | 0.0321 | **0.0145** | 0.0242 |

- The eigenbasis model cuts SBE error by more than half compared with the MHW-style model, and the 95% interval for the difference excludes zero (-0.0119 to -0.0006).
- The two models aren't clearly different on LIB or SB.
- The eigenbasis model isn't yet clearly better than the trivial baseline on SBE (-0.0127 to +0.0032). Fifteen winters leave wide intervals, because SBE varies greatly from winter to winter.

Supporting analyses found that:
- fitting red noise separately by region captures most of the geographic pattern of blocking;
- data-driven regions land close to the published ones;
- tomorrow's B depends far more on longitudes to the west than to the east, a signal neither model appears able to use.

## Final model

**Eigenbasis AR(1), provisionally.** It has the best primary KPI and clearly improves on the MHW-style model, without clearly losing on a secondary KPI. Its cost is complexity: about 2,600 fitted directions.

## Limitations

- **B definition.** The test comparison computes B at one central latitude. The published index takes the maximum over three, which gives much higher frequencies, so these numbers can't be compared directly with MHW's.
- **Different experiment.** Our models run for 5 days from observed states; MHW's runs freely. Our results don't directly test their claim that blocks sustain themselves.
- **One split, 15 test winters.** There's no cross-validation yet.
- **Limited to 5-day runs.** Seven of the eigenbasis model's coefficients are at or above 1, so it can't run freely for long.
- **Data.** ERA-40 is one superseded reanalysis, built from an observing system that changed over time.

## Next steps

1. Use the published B and SB definitions throughout, then rerun the comparison once.
2. Cross-validate within the training winters.
3. Test a model that lets information travel eastward.
4. Test state-dependent noise and memory, which directly address the "self-sustaining" idea.

## More detail

[README](../README.md), [data inventory](../data_inventory.md), [KPIs](../kpis.md), [evaluation plan](../evaluation_plan.md), [data audit](../data_audit.md), [modeling plan](../modeling_plan.md). Results: `results/model_comparison.csv`, from `notebooks/modeling_baselines.ipynb`.
