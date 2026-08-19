# Martingale R-learner simulation code

This repository contains code for simulation studies for the Martingale R-learner, a method for estimating time-varying heterogeneous treatment effects for right-censored time-to-event outcomes.

The code includes:

- simulation data generation under multiple settings
- Martingale R-learner with estimated nuisance functions
- Martingale R-learner with oracle nuisance functions for comparison
- benchmark meta-learners including S-learner and T-learner

## Overview

For each simulation replicate, the code:

1. generates one training dataset under a specified data-generating mechanism
2. generates or reuses a fixed independent test dataset
3. fits the Martingale R-learner and benchmark learners
4. evaluates performance by mean squared error on the test dataset

The Martingale R-learner implementation uses:

- cross-fitted nuisance estimation
- pseudo-data construction over ordered event times
- spline basis expansion for the target treatment effect function
- weighted ridge regression for target estimation

## Requirements

The code is written in R.

Required R packages:

- `tidyverse`
- `survival`
- `gam`
- `polspline`
- `caret`
- `glmnet`

You can install them in R with:

```r
install.packages(c(
  "tidyverse",
  "survival",
  "gam",
  "polspline",
  "caret",
  "glmnet"
))
```

## Simple example

A copy of a simulated dataset from simulation Setup A is provided in
`example/simulated_data.rds`.

The dataset contains the subject ID, effect modifier `z2`, treatment
indicator, observed follow-up time, and event indicator.

### 1. Load the simulated data

```r
dat <- readRDS("example/simulated_data.rds")

str(dat)
head(dat)
```

### 2. Run the Martingale R-learner

The functions used for this simple example are provided in
`R code/MRL_example_functions.R`.

```r
source("R code/MRL_example_functions.R")

mrl_fit <- fit_mrl(
  data = dat,
  K = 5,
  Q = 5,
  degree_of_freedom_grid = 5:9,
  batch_size = 50
)
```

The selected spline degree of freedom and ridge penalty parameter can be
viewed using:

```r
mrl_fit$degree_of_freedom
mrl_fit$lambda
```

The nuisance estimation includes a treatment propensity model and a
conditional survival model. The marginal survival function is obtained
by averaging the treatment-specific conditional survival functions using
the estimated propensity score.

### 3. Predict the estimated HTE

For example, consider three combinations of follow-up time and the effect
modifier `z2`:

```r
new_data <- data.frame(
  time = c(0.5, 1.0, 1.5),
  z2 = c(-1, 0, 1)
)

new_data$tau_hat <- predict_mrl(
  mrl_fit,
  newdata = new_data
)

new_data
```

`tau_hat` gives the estimated instantaneous heterogeneous treatment effect
at each specified combination of follow-up time and covariates.

## Full simulation scripts

The complete simulation code for the Martingale R-learner is available
below:

- [Setup A](R%20code/setup_A_revised.R)
- [Setup B](R%20code/setup_B_revised.R)
- [Setup C](R%20code/setup_C_revised.R)
- [Setup D](R%20code/setup_D_revised.R)
- [Setup E](R%20code/setup_E_revised.R)

### S-learner and T-learner

The simulation code for the benchmark S-learner and T-learner methods is
also provided:

- [Setup A](R%20code/S%20and%20T%20learner%20for%20setup%20A.rmd)
- [Setup B](R%20code/S%20and%20T%20learner%20for%20setup%20B.rmd)
- [Setup C](R%20code/S%20and%20T%20learner%20for%20setup%20C.rmd)
- [Setup D](R%20code/S%20and%20T%20learner%20for%20setup%20D.rmd)
- [Setup E](R%20code/S%20and%20T%20learner%20for%20setup%20E.rmd)
