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
result.all = rep(0,2)
names(result.all) = c("estimated mse", "oracle mse")

t85 = 1.922834
b0 = 1/2
n = 500

m.trt0 = tibble(trt = 0)
m.trt1 = tibble(trt = 1)

# f for uniroot
f.uniroot = function(t, a, b, c, d) {
  a*t^3 + b*t^2 + c*t + d
}

df = data.frame(id = 1:n,  
                z2 = runif(n, -2, 2),
                trt = rbinom(n, 1, 0.5), 
                u = runif(n, 0 ,1),
                censortime = runif(n, 0, 8.5)) %>% 
  rowwise() %>%
  mutate(a = 1/3,
         b = 1/2,
         c = b0 + z2 + z2^2,
         d = log(u),
         eventtime = 
           if_else(trt == 0, -log(u)/b0,
                   uniroot(f.uniroot, c(0, 100), 
                           a=a, b=b, c=c, d=d, tol=1e-10)$root)) %>%
  mutate(observedtime = min(eventtime, censortime),
         observedevent = as.integer(eventtime<=censortime)
  ) %>% 
  ungroup() %>% 
  arrange(observedtime)
# CV folds, used in glmnet, ridge or lasso
K = 5 # K folds for CV
cvSplits = createFolds(df$observedevent, k = K, list = F)
# Add fold id for CV
df = df %>% cbind(cvSplits)
# Cross-fitting folds
Q = 5 # Q folds for cross-fitting
cfSplits = createFolds(df$observedevent, k = Q, list = F)
# Add fold id for cross-fitting
df = df %>% cbind(cfSplits)

# clean df
df.raw = df
df = df %>%
  dplyr::select(
    id,
    z2,
    trt,
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
set.seed(1234)

test.data = data.frame(id = 1:n,  
                       z2 = runif(n, -2, 2),
                       trt = rbinom(n, 1, 0.5), 
                       u = runif(n, 0 ,1),
                       censortime = runif(n, 0, 8.5)) %>% 
  rowwise() %>%
  mutate(a = 1/3,
         b = 1/2,
         c = b0 + z2 + z2^2,
         d = log(u),
         eventtime = 
           if_else(trt == 0, -log(u)/b0,
                   uniroot(f.uniroot, c(0, 100), 
                           a=a, b=b, c=c, d=d, tol=1e-10)$root)) %>%
  mutate(observedtime = min(eventtime, censortime),
         observedevent = as.integer(eventtime<=censortime)
  ) %>% 
  ungroup() %>% 
  arrange(observedtime)

set.seed(NULL)

# estimated
# Q estimated models
pi_hat.list = vector("list", Q)
surv_c.list = vector("list", Q)
for (cf_fold in 1:Q) {
  df.m = df %>% 
    dplyr::filter(cfSplits != cf_fold) %>% 
    as.matrix()
  
  pi_hat.list[[cf_fold]] = gam(
    trt ~ s(z2, 3),
    data = df %>% dplyr::filter(cfSplits != cf_fold),
    family = "binomial"
  )
  
  surv_c.list[[cf_fold]] = hare(
    df.m[, "observedtime"],
    df.m[, "observedevent"],
    df.m[, c("z2", "trt")]
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
  z2.basis.obj = splines::ns(df$z2, df = degree_of_freedom)
  
  t.basis.test = predict(t.basis.obj, t.vector[1])
  z2.basis.test = predict(z2.basis.obj, df$z2[1])
  
  x.basis.test = cbind(t.basis.test, z2.basis.test)
  X_ns.test = cbind(x.basis.test, cproduct(t.basis.test, z2.basis.test))
  
  p_tau = ncol(cbind(1, X_ns.test))
  
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
        z2 = one$z2,
        trt = one$trt,
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
        c(df.id[1, c("z2")], m.trt0),
        surv_c.list[[cf_fold]]
      )
      
      survival1.id = 1 - phare(
        df.id$t_k,
        c(df.id[1, c("z2")], m.trt1),
        surv_c.list[[cf_fold]]
      )
      
      pi.id = predict.glm(
        pi_hat.list[[cf_fold]],
        newdata = tibble(z2 = df.id$z2[1]),
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
      z2.basis_i = predict(z2.basis.obj, df.id$z2)
      
      x.basis_i = cbind(t.basis_i, z2.basis_i)
      X_ns_i = cbind(x.basis_i, cproduct(t.basis_i, z2.basis_i))
      
      basis_i = cbind(1, X_ns_i)
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
    z2.basis.obj = z2.basis.obj
  )
}


# select best degree of freedom and lambda
cv_result

best_row = cv_result %>%
  dplyr::slice_min(cv_min, n = 1, with_ties = FALSE)

best_degree_of_freedom = best_row$degree_of_freedom
best_lambda = best_row$lambda_min

best_fit = fit_by_df[[paste0("df", best_degree_of_freedom)]]

best_degree_of_freedom
best_lambda
best_row


# final selected fit
theta_final = best_fit$theta_final
t.basis.obj = best_fit$t.basis.obj
z2.basis.obj = best_fit$z2.basis.obj


# mse on test data
df.sample = test.data %>%
  dplyr::filter(observedtime <= t85)

z.sample = df.sample$z2
t.sample = df.sample$observedtime

tau0 = z.sample + t.sample +
  z.sample^2 + t.sample^2

newt.basis = predict(t.basis.obj, t.sample)
newz2.basis = predict(z2.basis.obj, z.sample)

newx.basis = cbind(newt.basis, newz2.basis)
newX_ns = cbind(newx.basis, cproduct(newt.basis, newz2.basis))

tau_hat = as.vector(cbind(1, newX_ns) %*% theta_final)

MSE = mean((tau0 - tau_hat)^2)

result.all[1] = MSE

# oracle uses degree of freedom selected from estimated nuisance
oracle_degree_of_freedom = best_degree_of_freedom

t.basis.obj.oracle = splines::ns(
  t.vector,
  df = oracle_degree_of_freedom
)

z2.basis.obj.oracle = splines::ns(
  df$z2,
  df = oracle_degree_of_freedom
)


# get oracle basis dimension
t.basis.test = predict(t.basis.obj.oracle, t.vector[1])
z2.basis.test = predict(z2.basis.obj.oracle, df$z2[1])

x.basis.test = cbind(t.basis.test, z2.basis.test)
X_ns.test = cbind(x.basis.test, cproduct(t.basis.test, z2.basis.test))

p_tau_oracle = ncol(cbind(1, X_ns.test))


# oracle penalty matrix
penalty_factor_tau_oracle = c(0, rep(1, p_tau_oracle - 1))

P_tau_oracle = diag(
  penalty_factor_tau_oracle,
  nrow = p_tau_oracle,
  ncol = p_tau_oracle
)


# split subjects by cv folds, then batches
# If rows_by_cv_batch was already created in the estimated section,
# this block is not strictly necessary, but keeping it here makes
# the oracle section self-contained.
batch_size = 50

rows_by_cv = split(
  seq_len(nrow(df)),
  factor(df$cvSplits, levels = 1:K)
)

rows_by_cv_batch = lapply(rows_by_cv, function(rows) {
  split(rows, ceiling(seq_along(rows) / batch_size))
})


# compute oracle G and h for one computational batch
compute_Gh_batch_oracle = function(batch_rows) {
  G_batch = matrix(0, nrow = p_tau_oracle, ncol = p_tau_oracle)
  h_batch = rep(0, p_tau_oracle)
  
  for (rr in batch_rows) {
    one = df[rr, ]
    
    df.id = tibble(
      id = one$id,
      z2 = one$z2,
      trt = one$trt,
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
    
    # oracle rho and Lambda_m
    df.id = df.id %>%
      mutate(
        A_tk = (z2 + z2^2) * t_k + 1 / 2 * t_k^2 + 1 / 3 * t_k^3,
        expTerm = exp(-A_tk),
        rho = expTerm / (1 + expTerm),
        Lambda_tk = b0 * t_k - log(0.5 + 0.5 * expTerm),
        Lambda_delta = Lambda_tk - lag(Lambda_tk, default = 0)
      )
    
    # oracle pseudo y, w, x
    y_i = (df.id$counting - df.id$Lambda_delta) / df.id$delta_t
    w_i = df.id$delta_t
    
    t.basis_i = predict(t.basis.obj.oracle, df.id$t_k)
    z2.basis_i = predict(z2.basis.obj.oracle, df.id$z2)
    
    x.basis_i = cbind(t.basis_i, z2.basis_i)
    X_ns_i = cbind(x.basis_i, cproduct(t.basis_i, z2.basis_i))
    
    basis_i = cbind(1, X_ns_i)
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
df.sample = test.data %>%
  dplyr::filter(observedtime <= t85)

z.sample = df.sample$z2
t.sample = df.sample$observedtime

tau0 = z.sample + t.sample +
  z.sample^2 + t.sample^2

newt.basis = predict(t.basis.obj.oracle, t.sample)
newz2.basis = predict(z2.basis.obj.oracle, z.sample)

newx.basis = cbind(newt.basis, newz2.basis)
newX_ns = cbind(newx.basis, cproduct(newt.basis, newz2.basis))

tau_hat_oracle = as.vector(cbind(1, newX_ns) %*% theta_final_oracle)

MSE_oracle = mean((tau0 - tau_hat_oracle)^2)

result.all[2] = MSE_oracle

output = c(result.all, best_degree_of_freedom)
saveRDS(output, "result")
