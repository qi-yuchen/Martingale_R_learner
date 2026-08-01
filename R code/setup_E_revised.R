#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args)) # data seed

library(tidyverse)
library(survival)
library(gam)
library(polspline)
library(caret)
library(glmnet)

# Data
t85 = 0.4618303

result.all = rep(0, 2)
names(result.all) = c("estimated mse", "oracle mse")

b0 = 3
n = 500

m.trt0 = tibble(trt = 0)
m.trt1 = tibble(trt = 1)

# f for uniroot
f.uniroot = function(t, a, b, c) {
  a * t + b * cos(t) + c
}

df = data.frame(
  id = 1:n,
  z1 = runif(n, -2, 2),
  z2 = runif(n, -2, 2),
  z3 = runif(n, -2, 2),
  z4 = runif(n, -2, 2),
  z5 = runif(n, -2, 2),
  u = runif(n, 0, 1),
  censortime = runif(n, 0, 2)
) %>%
  mutate(
    trt.pi = 1 / (1 + exp(z1 + z2 + z3))
  ) %>%
  rowwise() %>%
  mutate(
    trt = rbinom(1, 1, trt.pi),
    tau_no_t = log(6 + z3 + z4) + sin(z5),
    a = b0 + z1 + trt * tau_no_t,
    b = -trt,
    c = trt + log(u),
    eventtime = if_else(
      trt == 0,
      -log(u) / (b0 + z1),
      uniroot(
        f.uniroot,
        c(0, 100),
        a = a,
        b = b,
        c = c,
        tol = 1e-10
      )$root
    )
  ) %>%
  mutate(
    observedtime = min(eventtime, censortime),
    observedevent = as.integer(eventtime <= censortime)
  ) %>%
  ungroup() %>%
  arrange(observedtime)

# CV folds, used for lambda and degree-of-freedom selection
K = 5
cvSplits = createFolds(df$observedevent, k = K, list = FALSE)
df = df %>% cbind(cvSplits)

# Cross-fitting folds, used for nuisance models
Q = 5
cfSplits = createFolds(df$observedevent, k = Q, list = FALSE)
df = df %>% cbind(cfSplits)

# clean df
df.raw = df

df = df %>%
  dplyr::select(
    id,
    z1,
    z2,
    z3,
    z4,
    z5,
    trt,
    trt.pi,
    observedtime,
    observedevent,
    cvSplits,
    cfSplits
  )

# time grid
t.vector = df %>%
  dplyr::filter(observedevent == 1) %>%
  pull(observedtime) %>%
  unique() %>%
  sort()

delta_t.vector = t.vector - lag(t.vector, default = 0)


# Test data
set.seed(123)

test.data = data.frame(
  id = 1:n,
  z1 = runif(n, -2, 2),
  z2 = runif(n, -2, 2),
  z3 = runif(n, -2, 2),
  z4 = runif(n, -2, 2),
  z5 = runif(n, -2, 2),
  u = runif(n, 0, 1),
  censortime = runif(n, 0, 2)
) %>%
  mutate(
    trt.pi = 1 / (1 + exp(z1 + z2 + z3))
  ) %>%
  rowwise() %>%
  mutate(
    trt = rbinom(1, 1, trt.pi),
    tau_no_t = log(6 + z3 + z4) + sin(z5),
    a = b0 + z1 + trt * tau_no_t,
    b = -trt,
    c = trt + log(u),
    eventtime = if_else(
      trt == 0,
      -log(u) / (b0 + z1),
      uniroot(
        f.uniroot,
        c(0, 100),
        a = a,
        b = b,
        c = c,
        tol = 1e-10
      )$root
    )
  ) %>%
  mutate(
    observedtime = min(eventtime, censortime),
    observedevent = as.integer(eventtime <= censortime)
  ) %>%
  ungroup() %>%
  arrange(observedtime)

set.seed(NULL)


# estimated nuisance models
covariates = c("z1", "z2", "z3", "z4", "z5")

pi_hat.list = vector("list", Q)
surv_c.list = vector("list", Q)

for (cf_fold in 1:Q) {
  df.m = df %>%
    dplyr::filter(cfSplits != cf_fold) %>%
    as.matrix()
  
  pi_hat.list[[cf_fold]] = gam(
    trt ~ s(z1, 3) + s(z2, 3) + s(z3, 3),
    data = df %>% dplyr::filter(cfSplits != cf_fold),
    family = "binomial"
  )
  
  surv_c.list[[cf_fold]] = hare(
    df.m[, "observedtime"],
    df.m[, "observedevent"],
    df.m[, c(covariates, "trt")]
  )
}


# interaction basis helper
cproduct = function(m1, m2) {
  m1 = as.matrix(m1)
  m2 = as.matrix(m2)
  
  out = do.call(
    cbind,
    lapply(seq_len(ncol(m1)), function(j) {
      m1[, j] * m2
    })
  )
  
  return(out)
}


# make design matrix
make_design_E = function(t.basis, z1.basis, z2.basis, z3.basis, z4.basis, z5.basis) {
  x.basis = cbind(t.basis, z1.basis, z2.basis, z3.basis, z4.basis, z5.basis)
  
  X_ns = cbind(
    x.basis,
    cproduct(t.basis, z1.basis),
    cproduct(t.basis, z2.basis),
    cproduct(t.basis, z3.basis),
    cproduct(z1.basis, z2.basis),
    cproduct(z2.basis, z3.basis),
    cproduct(z1.basis, z3.basis),
    cproduct(t.basis, z4.basis),
    cproduct(t.basis, z5.basis),
    cproduct(z1.basis, z4.basis),
    cproduct(z1.basis, z5.basis),
    cproduct(z2.basis, z4.basis),
    cproduct(z2.basis, z5.basis),
    cproduct(z3.basis, z4.basis),
    cproduct(z3.basis, z5.basis),
    cproduct(z4.basis, z5.basis)
  )
  
  basis = cbind(1, X_ns)
  
  return(basis)
}


# split subjects by cv folds, then batches
batch_size = 50

rows_by_cv = split(
  seq_len(nrow(df)),
  factor(df$cvSplits, levels = 1:K)
)

rows_by_cv_batch = lapply(rows_by_cv, function(rows) {
  split(rows, ceiling(seq_along(rows) / batch_size))
})


# candidate degree of freedom values
degree_of_freedom_grid = 5:9

cv_result = tibble(
  degree_of_freedom = integer(),
  p_tau = integer(),
  lambda_min = numeric(),
  cv_min = numeric(),
  lambda_index = integer()
)

fit_by_df = vector("list", length(degree_of_freedom_grid))
names(fit_by_df) = paste0("df", degree_of_freedom_grid)


# outer loop over basis degree of freedom
for (df_index in seq_along(degree_of_freedom_grid)) {
  
  degree_of_freedom = degree_of_freedom_grid[df_index]
  
  cat("Degree of freedom:", degree_of_freedom, "\n")
  
  # define basis for this degree of freedom
  t.basis.obj = splines::ns(t.vector, df = degree_of_freedom)
  z1.basis.obj = splines::ns(df$z1, df = degree_of_freedom)
  z2.basis.obj = splines::ns(df$z2, df = degree_of_freedom)
  z3.basis.obj = splines::ns(df$z3, df = degree_of_freedom)
  z4.basis.obj = splines::ns(df$z4, df = degree_of_freedom)
  z5.basis.obj = splines::ns(df$z5, df = degree_of_freedom)
  
  t.basis.test = predict(t.basis.obj, t.vector[1])
  z1.basis.test = predict(z1.basis.obj, df$z1[1])
  z2.basis.test = predict(z2.basis.obj, df$z2[1])
  z3.basis.test = predict(z3.basis.obj, df$z3[1])
  z4.basis.test = predict(z4.basis.obj, df$z4[1])
  z5.basis.test = predict(z5.basis.obj, df$z5[1])
  
  basis.test = make_design_E(
    t.basis = t.basis.test,
    z1.basis = z1.basis.test,
    z2.basis = z2.basis.test,
    z3.basis = z3.basis.test,
    z4.basis = z4.basis.test,
    z5.basis = z5.basis.test
  )
  
  p_tau = ncol(basis.test)
  
  # penalty matrix
  penalty_factor_tau = c(0, rep(1, p_tau - 1))
  P_tau = diag(penalty_factor_tau, nrow = p_tau, ncol = p_tau)
  
  # compute G and h for one computational batch
  compute_Gh_batch = function(batch_rows) {
    G_batch = matrix(0, nrow = p_tau, ncol = p_tau)
    h_batch = rep(0, p_tau)
    
    eps = 1e-8
    
    for (rr in batch_rows) {
      one = df[rr, ]
      
      df.id = tibble(
        id = one$id,
        z1 = one$z1,
        z2 = one$z2,
        z3 = one$z3,
        z4 = one$z4,
        z5 = one$z5,
        trt = one$trt,
        trt.pi = one$trt.pi,
        observedtime = one$observedtime,
        observedevent = one$observedevent,
        cvSplits = one$cvSplits,
        cfSplits = one$cfSplits,
        k = seq_along(t.vector),
        t_k = t.vector,
        delta_t = delta_t.vector
      ) %>%
        dplyr::filter(t_k <= observedtime) %>%
        mutate(
          counting = if_else(
            observedevent == 1 & abs(t_k - observedtime) < 1e-10,
            1,
            0
          )
        )
      
      if (nrow(df.id) == 0) next
      
      cf_fold = df.id$cfSplits[1]
      
      survival0.id = 1 - phare(
        df.id$t_k,
        c(df.id[1, covariates], m.trt0),
        surv_c.list[[cf_fold]]
      )
      
      survival1.id = 1 - phare(
        df.id$t_k,
        c(df.id[1, covariates], m.trt1),
        surv_c.list[[cf_fold]]
      )
      
      pi.id = predict.glm(
        pi_hat.list[[cf_fold]],
        newdata = tibble(
          z1 = df.id$z1[1],
          z2 = df.id$z2[1],
          z3 = df.id$z3[1]
        ),
        type = "response"
      )
      
            survival.id =
        (1 - as.numeric(pi.id)) * survival0.id +
        as.numeric(pi.id) * survival1.id

df.id = df.id %>%
        mutate(
          survival_tk = survival.id,
          survival0 = survival0.id,
          survival1 = survival1.id,
          pi = as.numeric(pi.id)
        ) %>%
        mutate(
          survival_tk = pmin(pmax(survival_tk, eps), 1 - eps),
          survival0 = pmin(pmax(survival0, eps), 1 - eps),
          survival1 = pmin(pmax(survival1, eps), 1 - eps),
          rho = pi * survival1 / survival0 /
            (1 + (survival1 / survival0 - 1) * pi),
          Lambda_tk = -log(survival_tk),
          Lambda_tkpre = lag(Lambda_tk, default = 0),
          Lambda_delta = Lambda_tk - Lambda_tkpre
        )
      
      y_i = (df.id$counting - df.id$Lambda_delta) / df.id$delta_t
      w_i = df.id$delta_t
      
      t.basis_i = predict(t.basis.obj, df.id$t_k)
      z1.basis_i = predict(z1.basis.obj, df.id$z1)
      z2.basis_i = predict(z2.basis.obj, df.id$z2)
      z3.basis_i = predict(z3.basis.obj, df.id$z3)
      z4.basis_i = predict(z4.basis.obj, df.id$z4)
      z5.basis_i = predict(z5.basis.obj, df.id$z5)
      
      basis_i = make_design_E(
        t.basis = t.basis_i,
        z1.basis = z1.basis_i,
        z2.basis = z2.basis_i,
        z3.basis = z3.basis_i,
        z4.basis = z4.basis_i,
        z5.basis = z5.basis_i
      )
      
      x_i = basis_i * as.numeric(df.id$trt - df.id$rho)
      
      G_batch = G_batch + crossprod(x_i * sqrt(w_i))
      h_batch = h_batch + as.vector(crossprod(x_i, w_i * y_i))
    }
    
    return(list(
      G = G_batch,
      h = h_batch
    ))
  }
  
  
  # compute G and h for each cv fold
  G_by_cv = vector("list", K)
  h_by_cv = vector("list", K)
  
  for (cv_fold in 1:K) {
    cat("  CV fold:", cv_fold, "\n")
    
    G_cv = matrix(0, nrow = p_tau, ncol = p_tau)
    h_cv = rep(0, p_tau)
    
    for (batch_id in seq_along(rows_by_cv_batch[[cv_fold]])) {
      cat("    batch:", batch_id, "\n")
      
      batch_rows = rows_by_cv_batch[[cv_fold]][[batch_id]]
      
      Gh_batch = compute_Gh_batch(batch_rows)
      
      G_cv = G_cv + Gh_batch$G
      h_cv = h_cv + Gh_batch$h
    }
    
    G_by_cv[[cv_fold]] = G_cv
    h_by_cv[[cv_fold]] = h_cv
  }
  
  
  # combine all cv folds
  G_all = Reduce("+", G_by_cv)
  h_all = Reduce("+", h_by_cv)
  
  
  # lambda grid for this degree of freedom
  lambda_base = mean(diag(G_all)[penalty_factor_tau == 1])
  
  if (!is.finite(lambda_base) || lambda_base <= 0) {
    lambda_base = 1
  }
  
  lambda.grid = lambda_base * 10^seq(-6, 2, length.out = 80)
  
  
  # cv over lambda for this degree of freedom
  cv.score = rep(NA, length(lambda.grid))
  
  for (lambda_index in seq_along(lambda.grid)) {
    lambda_value = lambda.grid[lambda_index]
    
    score_by_cv = rep(NA, K)
    
    for (cv_fold in 1:K) {
      G_valid = G_by_cv[[cv_fold]]
      h_valid = h_by_cv[[cv_fold]]
      
      G_train = G_all - G_valid
      h_train = h_all - h_valid
      
      theta_v = qr.solve(
        G_train + lambda_value * P_tau,
        h_train
      )
      
      score_v = as.numeric(
        crossprod(theta_v, G_valid %*% theta_v) -
          2 * crossprod(h_valid, theta_v)
      )
      
      score_by_cv[cv_fold] = score_v
    }
    
    cv.score[lambda_index] = sum(score_by_cv)
  }
  
  lambda.min = lambda.grid[which.min(cv.score)]
  cv.min = min(cv.score)
  lambda.index = which.min(cv.score)
  
  
  # final fit for this degree of freedom
  theta_final = qr.solve(
    G_all + lambda.min * P_tau,
    h_all
  )
  
  
  # save result for this degree of freedom
  cv_result = bind_rows(
    cv_result,
    tibble(
      degree_of_freedom = degree_of_freedom,
      p_tau = p_tau,
      lambda_min = lambda.min,
      cv_min = cv.min,
      lambda_index = lambda.index
    )
  )
  
  fit_by_df[[paste0("df", degree_of_freedom)]] = list(
    degree_of_freedom = degree_of_freedom,
    p_tau = p_tau,
    penalty_factor_tau = penalty_factor_tau,
    P_tau = P_tau,
    G_by_cv = G_by_cv,
    h_by_cv = h_by_cv,
    G_all = G_all,
    h_all = h_all,
    lambda.grid = lambda.grid,
    cv.score = cv.score,
    lambda.min = lambda.min,
    theta_final = theta_final,
    t.basis.obj = t.basis.obj,
    z1.basis.obj = z1.basis.obj,
    z2.basis.obj = z2.basis.obj,
    z3.basis.obj = z3.basis.obj,
    z4.basis.obj = z4.basis.obj,
    z5.basis.obj = z5.basis.obj
  )
}


# select best degree of freedom and lambda
best_row = cv_result %>%
  dplyr::slice_min(cv_min, n = 1, with_ties = FALSE)

best_degree_of_freedom = best_row$degree_of_freedom
best_lambda = best_row$lambda_min

best_fit = fit_by_df[[paste0("df", best_degree_of_freedom)]]

theta_final = best_fit$theta_final
t.basis.obj = best_fit$t.basis.obj
z1.basis.obj = best_fit$z1.basis.obj
z2.basis.obj = best_fit$z2.basis.obj
z3.basis.obj = best_fit$z3.basis.obj
z4.basis.obj = best_fit$z4.basis.obj
z5.basis.obj = best_fit$z5.basis.obj


# mse on test data
df.sample = test.data %>%
  dplyr::filter(observedtime <= t85)

z1.sample = df.sample$z1
z2.sample = df.sample$z2
z3.sample = df.sample$z3
z4.sample = df.sample$z4
z5.sample = df.sample$z5
t.sample = df.sample$observedtime

tau0 = log(6 + z3.sample + z4.sample) + sin(z5.sample) + sin(t.sample)

newt.basis = predict(t.basis.obj, t.sample)
newz1.basis = predict(z1.basis.obj, z1.sample)
newz2.basis = predict(z2.basis.obj, z2.sample)
newz3.basis = predict(z3.basis.obj, z3.sample)
newz4.basis = predict(z4.basis.obj, z4.sample)
newz5.basis = predict(z5.basis.obj, z5.sample)

new_basis = make_design_E(
  t.basis = newt.basis,
  z1.basis = newz1.basis,
  z2.basis = newz2.basis,
  z3.basis = newz3.basis,
  z4.basis = newz4.basis,
  z5.basis = newz5.basis
)

tau_hat = as.vector(new_basis %*% theta_final)

MSE = mean((tau0 - tau_hat)^2)

result.all[1] = MSE


# oracle uses degree of freedom selected from estimated nuisance
oracle_degree_of_freedom = best_degree_of_freedom

t.basis.obj.oracle = splines::ns(t.vector, df = oracle_degree_of_freedom)
z1.basis.obj.oracle = splines::ns(df$z1, df = oracle_degree_of_freedom)
z2.basis.obj.oracle = splines::ns(df$z2, df = oracle_degree_of_freedom)
z3.basis.obj.oracle = splines::ns(df$z3, df = oracle_degree_of_freedom)
z4.basis.obj.oracle = splines::ns(df$z4, df = oracle_degree_of_freedom)
z5.basis.obj.oracle = splines::ns(df$z5, df = oracle_degree_of_freedom)


# get oracle basis dimension
t.basis.test = predict(t.basis.obj.oracle, t.vector[1])
z1.basis.test = predict(z1.basis.obj.oracle, df$z1[1])
z2.basis.test = predict(z2.basis.obj.oracle, df$z2[1])
z3.basis.test = predict(z3.basis.obj.oracle, df$z3[1])
z4.basis.test = predict(z4.basis.obj.oracle, df$z4[1])
z5.basis.test = predict(z5.basis.obj.oracle, df$z5[1])

basis.test = make_design_E(
  t.basis = t.basis.test,
  z1.basis = z1.basis.test,
  z2.basis = z2.basis.test,
  z3.basis = z3.basis.test,
  z4.basis = z4.basis.test,
  z5.basis = z5.basis.test
)

p_tau_oracle = ncol(basis.test)


# oracle penalty matrix
penalty_factor_tau_oracle = c(0, rep(1, p_tau_oracle - 1))

P_tau_oracle = diag(
  penalty_factor_tau_oracle,
  nrow = p_tau_oracle,
  ncol = p_tau_oracle
)


# compute oracle G and h for one computational batch
compute_Gh_batch_oracle = function(batch_rows) {
  G_batch = matrix(0, nrow = p_tau_oracle, ncol = p_tau_oracle)
  h_batch = rep(0, p_tau_oracle)
  
  for (rr in batch_rows) {
    one = df[rr, ]
    
    df.id = tibble(
      id = one$id,
      z1 = one$z1,
      z2 = one$z2,
      z3 = one$z3,
      z4 = one$z4,
      z5 = one$z5,
      trt = one$trt,
      trt.pi = one$trt.pi,
      observedtime = one$observedtime,
      observedevent = one$observedevent,
      cvSplits = one$cvSplits,
      cfSplits = one$cfSplits,
      k = seq_along(t.vector),
      t_k = t.vector,
      delta_t = delta_t.vector
    ) %>%
      dplyr::filter(t_k <= observedtime) %>%
      mutate(
        counting = if_else(
          observedevent == 1 & abs(t_k - observedtime) < 1e-10,
          1,
          0
        )
      )
    
    if (nrow(df.id) == 0) next
    
    df.id = df.id %>%
      mutate(
        tau_no_t = log(6 + z3 + z4) + sin(z5),
        expTerm = exp(-(tau_no_t * t_k + 1 - cos(t_k))),
        rho = trt.pi * expTerm / (1 + (expTerm - 1) * trt.pi),
        hazardrate_tk = b0 + z1 + (tau_no_t + sin(t_k)) * rho,
        Lambda_delta = hazardrate_tk * delta_t
      )
    
    y_i = (df.id$counting - df.id$Lambda_delta) / df.id$delta_t
    w_i = df.id$delta_t
    
    t.basis_i = predict(t.basis.obj.oracle, df.id$t_k)
    z1.basis_i = predict(z1.basis.obj.oracle, df.id$z1)
    z2.basis_i = predict(z2.basis.obj.oracle, df.id$z2)
    z3.basis_i = predict(z3.basis.obj.oracle, df.id$z3)
    z4.basis_i = predict(z4.basis.obj.oracle, df.id$z4)
    z5.basis_i = predict(z5.basis.obj.oracle, df.id$z5)
    
    basis_i = make_design_E(
      t.basis = t.basis_i,
      z1.basis = z1.basis_i,
      z2.basis = z2.basis_i,
      z3.basis = z3.basis_i,
      z4.basis = z4.basis_i,
      z5.basis = z5.basis_i
    )
    
    x_i = basis_i * as.numeric(df.id$trt - df.id$rho)
    
    G_batch = G_batch + crossprod(x_i * sqrt(w_i))
    h_batch = h_batch + as.vector(crossprod(x_i, w_i * y_i))
  }
  
  return(list(
    G = G_batch,
    h = h_batch
  ))
}


# compute oracle G and h for each cv fold
G_by_cv_oracle = vector("list", K)
h_by_cv_oracle = vector("list", K)

for (cv_fold in 1:K) {
  cat("Oracle CV fold:", cv_fold, "\n")
  
  G_cv = matrix(0, nrow = p_tau_oracle, ncol = p_tau_oracle)
  h_cv = rep(0, p_tau_oracle)
  
  for (batch_id in seq_along(rows_by_cv_batch[[cv_fold]])) {
    cat("  batch:", batch_id, "\n")
    
    batch_rows = rows_by_cv_batch[[cv_fold]][[batch_id]]
    
    Gh_batch = compute_Gh_batch_oracle(batch_rows)
    
    G_cv = G_cv + Gh_batch$G
    h_cv = h_cv + Gh_batch$h
  }
  
  G_by_cv_oracle[[cv_fold]] = G_cv
  h_by_cv_oracle[[cv_fold]] = h_cv
}


# combine oracle cv folds
G_all_oracle = Reduce("+", G_by_cv_oracle)
h_all_oracle = Reduce("+", h_by_cv_oracle)


# oracle lambda grid
lambda_base_oracle = mean(
  diag(G_all_oracle)[penalty_factor_tau_oracle == 1]
)

if (!is.finite(lambda_base_oracle) || lambda_base_oracle <= 0) {
  lambda_base_oracle = 1
}

lambda.grid.oracle = lambda_base_oracle * 10^seq(-6, 2, length.out = 80)


# oracle cv over lambda, using pseudo-loss metric
cv.score.oracle = rep(NA, length(lambda.grid.oracle))

for (lambda_index in seq_along(lambda.grid.oracle)) {
  lambda_value = lambda.grid.oracle[lambda_index]
  
  score_by_cv = rep(NA, K)
  
  for (cv_fold in 1:K) {
    G_valid = G_by_cv_oracle[[cv_fold]]
    h_valid = h_by_cv_oracle[[cv_fold]]
    
    G_train = G_all_oracle - G_valid
    h_train = h_all_oracle - h_valid
    
    theta_v = qr.solve(
      G_train + lambda_value * P_tau_oracle,
      h_train
    )
    
    score_v = as.numeric(
      crossprod(theta_v, G_valid %*% theta_v) -
        2 * crossprod(h_valid, theta_v)
    )
    
    score_by_cv[cv_fold] = score_v
  }
  
  cv.score.oracle[lambda_index] = sum(score_by_cv)
}

lambda.min.oracle = lambda.grid.oracle[which.min(cv.score.oracle)]


# final oracle fit using selected oracle lambda
theta_final_oracle = qr.solve(
  G_all_oracle + lambda.min.oracle * P_tau_oracle,
  h_all_oracle
)


# oracle mse on test data
newt.basis = predict(t.basis.obj.oracle, t.sample)
newz1.basis = predict(z1.basis.obj.oracle, z1.sample)
newz2.basis = predict(z2.basis.obj.oracle, z2.sample)
newz3.basis = predict(z3.basis.obj.oracle, z3.sample)
newz4.basis = predict(z4.basis.obj.oracle, z4.sample)
newz5.basis = predict(z5.basis.obj.oracle, z5.sample)

new_basis_oracle = make_design_E(
  t.basis = newt.basis,
  z1.basis = newz1.basis,
  z2.basis = newz2.basis,
  z3.basis = newz3.basis,
  z4.basis = newz4.basis,
  z5.basis = newz5.basis
)

tau_hat_oracle = as.vector(new_basis_oracle %*% theta_final_oracle)

MSE_oracle = mean((tau0 - tau_hat_oracle)^2)

result.all[2] = MSE_oracle


# save output
output = c(result.all, best_degree_of_freedom)
saveRDS(output, "result")