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
