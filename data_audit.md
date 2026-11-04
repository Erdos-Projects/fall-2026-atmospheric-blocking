# Data audit and feature decisions

Numbers come from `notebooks/data_audit.ipynb` and `notebooks/feature_selection.ipynb`. File checks cover the whole record. Everything else uses the 30 training winters (December 1957 to February 1987, 2,707 days), so no decisions rest on test winters. B is the Tyrlis & Hoskins version; see open item 1 in [evaluation_plan.md](evaluation_plan.md). These notebooks were written at the end of the project, using hindsight about which questions mattered.

`notebooks/eda.ipynb` holds the earlier exploration on the 4-winter excerpt: theta distributions, PCA and blocking animations.

## Data audit

### Completeness

- All 540 monthly files are present, September 1957 to August 2002, each with all 7 variables at every 6-hourly time.
- No missing values: no message flags missing data and no theta value is NaN.
- Theta ranges from 225 to 448 K over the whole record, which is physically plausible at the 2 PVU surface.
- Daily B has 16,436 days at 72 longitudes, with no gaps or NaNs.

Nothing needed imputation.

### Distribution of B

| Across the 72 longitudes | Lowest | Median | Highest |
|---|---|---|---|
| Mean (K) | -15.3 | -9.6 | -3.0 |
| Standard deviation (K) | 7.4 | 9.9 | 11.7 |
| Skewness | -0.64 | -0.04 | +0.44 |
| Excess kurtosis | -0.24 | +0.16 | +1.21 |
| Share of days with B > 0 | 0.07 | 0.18 | 0.36 |

Mean B varies far more by longitude than its spread does, which is why blocking concentrates at certain longitudes. Spread and skewness also vary, so a single symmetric noise distribution can't fit every longitude. This is visible in `results/eda/ERA40_B_train_moments_by_longitude.png`.

### Outliers

About 0.03% of values lie more than 4 standard deviations from their longitude's mean, compared with 0.006% for a normal distribution. They split evenly between high (26) and low (30). Since B ranges only from -66 to +36 K and theta passes the physical range check, we treat them as real extreme flows and removed none.

### Targets

| Measure | Mean over longitudes | Peak | Longitudes below 1% |
|---|---|---|---|
| LIB | 0.192 | 0.356 at 10°E | 0 |
| SB | 0.273 | 0.454 at 20°E | 0 |
| SBE | 0.060 | 0.145 at 25°E | 9 |

- SBE is rare at some longitudes, so small absolute errors there are large relative errors.
- SBE varies greatly between winters: at its peak longitude the winter-to-winter standard deviation is 0.084, more than half the mean. Fifteen test winters therefore carry real sampling uncertainty.

### Representativeness

Mean B in 1979-86 differs from 1957-78 by -0.10 K averaged over longitudes, and by -2.8 to +2.1 K at individual longitudes. The second group is short, so this doesn't separate changes in the observing system from natural variability. There is no large overall shift.

### Correlations

| Separation between longitudes | 5° | 10° | 15° | 20° | 25° | 30° |
|---|---|---|---|---|---|---|
| Mean correlation of B | 0.90 | 0.69 | 0.48 | 0.30 | 0.18 | 0.10 |

Longitudes within about 25° carry overlapping information, so the 72 longitudes are not independent samples.

### Lag structure

| Lag (days) | 1 | 2 | 3 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|
| Autocorrelation | 0.643 | 0.376 | 0.284 | 0.183 | 0.067 | 0.027 |
| AR(1) prediction | 0.643 | 0.413 | 0.266 | 0.110 | 0.012 | 0.000 |

- Partial autocorrelation is 0.643 at lag 1, -0.063 at lag 2 and +0.115 at lag 3, then small.
- Autocorrelation falls faster than AR(1) at lag 2, but much more slowly beyond about 5 days. Long-lasting behavior is not captured by a single AR(1) rate.
- Lag-1 autocorrelation varies from 0.59 to 0.73 across longitudes, highest at 135°E.

## Feature decisions

The project asks how far theta on 2 PVU alone can explain blocking, so scope decides first. As supporting evidence, `notebooks/feature_selection.ipynb` predicts tomorrow's B at each longitude from today's B plus one candidate, by linear regression, scored on held-out winters (5-fold by winter).

Today's B alone gives a median out-of-sample R² of 0.43. All candidates and reasons are in `results/feature_selection.csv`.

### Kept

| Feature | R² gain | Why |
|---|---|---|
| Theta on 2 PVU | n/a | The quantity the project models; B is built from it |

### Not used

| Feature | R² gain | Why |
|---|---|---|
| B today, 5-15° west | 0.11-0.14 | The strongest signal found, but current models can't use its direction (see below) |
| Northward wind, north minus south | 0.045 | Out of scope: useful, but the model would have to simulate wind too. A possible extension |
| Eastward wind, band average | 0.029 | Out of scope, and partly redundant with B (correlation 0.76). Used only to locate the central latitude |
| Other wind indices | 0.003-0.006 | Little gain; out of scope |
| B today, 5-15° east | 0.000-0.022 | Far weaker than the same distance west |
| Geopotential, north minus south | 0.010 | Redundant with B (correlation 0.88) |
| Pressure, north minus south | 0.009 | Redundant with B (correlation 0.86) |
| B yesterday (lag 2) | 0.003 | Agrees with AR(2) fitting a negative second coefficient |
| Humidity, ozone, other band averages | 0.000-0.004 | Little or no gain |

**Limits of this test:** it is linear, predicts one day ahead rather than a 5-day episode, and adds one candidate at a time.

**The west-east difference matters for modeling.** Information about tomorrow comes from the west (gain 0.14 at 10° west, against 0.007 at 10° east), consistent with patterns moving eastward. Both current models update with a symmetric rule (a single φ, or eigenvectors with independent coefficients), so as far as the code shows, neither can represent that direction. This is the largest unused signal found here, and a natural next step for the models.

### Engineered features

| Feature | Rationale |
|---|---|
| B | Published index that summarizes reversal of the north-south theta pattern in one number per longitude |
| LIB, SB, SBE | Published targets that add spatial and then time constraints |
| Central latitude, fixed and seasonal | Places B where storms are most active; seasonal version takes every ingredient from the season's own data |
| 2-6 day band-passed wind variance | Isolates storm-length fluctuations, for the central latitude |
| Anomalies from climatology | Removes each longitude's or calendar day's average, which models don't need to learn |
| Eigenbasis of daily theta changes | Decorrelates the field, so independent AR(1)s can be fitted to each direction |
| Latitude restriction, 25-76°N and B-weighted cells | Limits the state to the cells that determine B |
| Regions (Tyrlis & Hoskins; data-driven) | Lets memory and spread vary by region with few parameters (1-D models) |

### Tried and discarded

| Idea | Outcome |
|---|---|
| AR(2) from lags 1 and 2 | Second coefficient negative, so less memory than AR(1), not more |
| Two-timescale AR(1) | Fits autocorrelation well, but no clear SBE improvement; dropped as too complex |
| Three data-driven regions | Too coarse: lumps Europe with Asia, worse than seven |
| Per-longitude parameters | Rejected: too close to fitting the data directly |
| EOF-filtered theta for detecting episodes | Episodes didn't match those from unfiltered theta (`eda.ipynb`) |

The schema of inputs, features and targets is in [schema.json](schema.json).
