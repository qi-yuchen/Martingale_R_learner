# Helper functions for the simple Martingale R-learner example.
# This is the estimated-nuisance portion of setup_A_revised.R wrapped
# into reusable fit/predict functions. It omits the oracle calculation
# and simulation MSE evaluation used only for the full simulation study.

cproduct <- function(m1, m2) {
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)

  do.call(
    cbind,
    lapply(seq_len(ncol(m1)), function(j) {
      m1[, j] * m2
    })
  )
}

fit_mrl <- function(
    data,
    K = 5,
    Q = 5,
    degree_of_freedom_grid = 5:9,
    batch_size = 50,
    fold_seed = 2026) {

  required_columns <- c(
    "id", "z2", "trt", "observedtime", "observedevent"
  )

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  df <- data[, required_columns]

  # CV folds for selecting the spline degree of freedom and ridge penalty.
  set.seed(fold_seed)
  df$cvSplits <- caret::createFolds(
    df$observedevent,
    k = K,
    list = FALSE
  )

  # Cross-fitting folds for nuisance estimation.
  df$cfSplits <- caret::createFolds(
    df$observedevent,
    k = Q,
    list = FALSE
  )

  # Event-time grid.
  t.vector <- sort(unique(
    df$observedtime[df$observedevent == 1]
  ))

  delta_t.vector <- t.vector - c(0, head(t.vector, -1))

  # Keep the original gam::gam smooth-term syntax while avoiding
  # a global library(gam) call.
  s <- gam::s

  # Cross-fitted nuisance models:
  # 1) treatment propensity pi(Z)
  # 2) conditional survival S(t | Z, D)
  pi_hat.list <- vector("list", Q)
  surv_c.list <- vector("list", Q)

  for (cf_fold in seq_len(Q)) {
    train_df <- df[df$cfSplits != cf_fold, ]

    pi_hat.list[[cf_fold]] <- gam::gam(
      trt ~ s(z2, 3),
      data = train_df,
      family = "binomial"
    )

    surv_c.list[[cf_fold]] <- polspline::hare(
      train_df$observedtime,
      train_df$observedevent,
      as.matrix(train_df[, c("z2", "trt")])
    )
  }

  # Split subjects by CV fold and then into computational batches.
  rows_by_cv <- split(
    seq_len(nrow(df)),
    factor(df$cvSplits, levels = seq_len(K))
  )

  rows_by_cv_batch <- lapply(rows_by_cv, function(rows) {
    split(rows, ceiling(seq_along(rows) / batch_size))
  })

  cv_result <- data.frame(
    degree_of_freedom = integer(),
    p_tau = integer(),
    lambda_min = numeric(),
    cv_min = numeric(),
    lambda_index = integer()
  )

  fit_by_df <- vector("list", length(degree_of_freedom_grid))
  names(fit_by_df) <- paste0("df", degree_of_freedom_grid)

  for (df_index in seq_along(degree_of_freedom_grid)) {
    degree_of_freedom <- degree_of_freedom_grid[df_index]

    cat("Degree of freedom:", degree_of_freedom, "\n")

    # Natural-spline basis for t and z2, including their tensor interaction.
    t.basis.obj <- splines::ns(
      t.vector,
      df = degree_of_freedom
    )

    z2.basis.obj <- splines::ns(
      df$z2,
      df = degree_of_freedom
    )

    t.basis.test <- predict(t.basis.obj, t.vector[1])
    z2.basis.test <- predict(z2.basis.obj, df$z2[1])

    x.basis.test <- cbind(t.basis.test, z2.basis.test)
    X_ns.test <- cbind(
      x.basis.test,
      cproduct(t.basis.test, z2.basis.test)
    )

    p_tau <- ncol(cbind(1, X_ns.test))

    # Ridge penalty; do not penalize the intercept.
    penalty_factor_tau <- c(0, rep(1, p_tau - 1))
    P_tau <- diag(
      penalty_factor_tau,
      nrow = p_tau,
      ncol = p_tau
    )

    compute_Gh_batch <- function(batch_rows) {
      G_batch <- matrix(0, nrow = p_tau, ncol = p_tau)
      h_batch <- rep(0, p_tau)
      eps <- 1e-8

      for (rr in batch_rows) {
        one <- df[rr, ]

        keep <- t.vector <= one$observedtime
        if (!any(keep)) next

        t_k <- t.vector[keep]
        delta_t <- delta_t.vector[keep]

        counting <- as.integer(
          one$observedevent == 1 &
            abs(t_k - one$observedtime) < 1e-10
        )

        cf_fold <- one$cfSplits

        # Conditional survival under treatment 0 and 1.
        survival0 <- 1 - polspline::phare(
          t_k,
          c(one$z2, 0),
          surv_c.list[[cf_fold]]
        )

        survival1 <- 1 - polspline::phare(
          t_k,
          c(one$z2, 1),
          surv_c.list[[cf_fold]]
        )

        pi_i <- as.numeric(stats::predict.glm(
          pi_hat.list[[cf_fold]],
          newdata = data.frame(z2 = one$z2),
          type = "response"
        ))

        # Marginal survival is derived from the conditional survival model
        # by averaging over treatment using the estimated propensity.
        survival_m <-
          (1 - pi_i) * survival0 + pi_i * survival1

        survival_m <- pmin(pmax(survival_m, eps), 1 - eps)
        survival0 <- pmin(pmax(survival0, eps), 1 - eps)
        survival1 <- pmin(pmax(survival1, eps), 1 - eps)

        rho <- pi_i * survival1 / survival0 /
          (1 + (survival1 / survival0 - 1) * pi_i)

        Lambda_tk <- -log(survival_m)
        Lambda_delta <- Lambda_tk - c(0, head(Lambda_tk, -1))

        # Martingale R-learner pseudo-response/design.
        y_i <- (counting - Lambda_delta) / delta_t
        w_i <- delta_t

        t.basis_i <- predict(t.basis.obj, t_k)
        z2.basis_i <- predict(
          z2.basis.obj,
          rep(one$z2, length(t_k))
        )

        x.basis_i <- cbind(t.basis_i, z2.basis_i)
        X_ns_i <- cbind(
          x.basis_i,
          cproduct(t.basis_i, z2.basis_i)
        )

        basis_i <- cbind(1, X_ns_i)
        x_i <- basis_i * as.numeric(one$trt - rho)

        G_batch <- G_batch + crossprod(x_i * sqrt(w_i))
        h_batch <- h_batch + as.vector(
          crossprod(x_i, w_i * y_i)
        )
      }

      list(G = G_batch, h = h_batch)
    }

    # Fold-specific sufficient statistics.
    G_by_cv <- vector("list", K)
    h_by_cv <- vector("list", K)

    for (cv_fold in seq_len(K)) {
      cat("  CV fold:", cv_fold, "\n")

      G_cv <- matrix(0, nrow = p_tau, ncol = p_tau)
      h_cv <- rep(0, p_tau)

      for (batch_id in seq_along(rows_by_cv_batch[[cv_fold]])) {
        batch_rows <- rows_by_cv_batch[[cv_fold]][[batch_id]]
        Gh_batch <- compute_Gh_batch(batch_rows)

        G_cv <- G_cv + Gh_batch$G
        h_cv <- h_cv + Gh_batch$h
      }

      G_by_cv[[cv_fold]] <- G_cv
      h_by_cv[[cv_fold]] <- h_cv
    }

    G_all <- Reduce("+", G_by_cv)
    h_all <- Reduce("+", h_by_cv)

    # Candidate ridge penalties.
    lambda_base <- mean(
      diag(G_all)[penalty_factor_tau == 1]
    )

    if (!is.finite(lambda_base) || lambda_base <= 0) {
      lambda_base <- 1
    }

    lambda.grid <- lambda_base * 10^seq(
      -6,
      2,
      length.out = 80
    )

    # CV over lambda using the pseudo-loss.
    cv.score <- rep(NA_real_, length(lambda.grid))

    for (lambda_index in seq_along(lambda.grid)) {
      lambda_value <- lambda.grid[lambda_index]
      score_by_cv <- rep(NA_real_, K)

      for (cv_fold in seq_len(K)) {
        G_valid <- G_by_cv[[cv_fold]]
        h_valid <- h_by_cv[[cv_fold]]

        G_train <- G_all - G_valid
        h_train <- h_all - h_valid

        theta_v <- qr.solve(
          G_train + lambda_value * P_tau,
          h_train
        )

        score_by_cv[cv_fold] <- as.numeric(
          crossprod(theta_v, G_valid %*% theta_v) -
            2 * crossprod(h_valid, theta_v)
        )
      }

      cv.score[lambda_index] <- sum(score_by_cv)
    }

    lambda.index <- which.min(cv.score)
    lambda.min <- lambda.grid[lambda.index]
    cv.min <- cv.score[lambda.index]

    theta_final <- qr.solve(
      G_all + lambda.min * P_tau,
      h_all
    )

    cv_result <- rbind(
      cv_result,
      data.frame(
        degree_of_freedom = degree_of_freedom,
        p_tau = p_tau,
        lambda_min = lambda.min,
        cv_min = cv.min,
        lambda_index = lambda.index
      )
    )

    fit_by_df[[paste0("df", degree_of_freedom)]] <- list(
      degree_of_freedom = degree_of_freedom,
      lambda.min = lambda.min,
      theta_final = theta_final,
      t.basis.obj = t.basis.obj,
      z2.basis.obj = z2.basis.obj
    )
  }

  # Select the spline dimension and ridge penalty with minimum CV loss.
  best_index <- which.min(cv_result$cv_min)
  best_row <- cv_result[best_index, ]
  best_fit <- fit_by_df[[
    paste0("df", best_row$degree_of_freedom)
  ]]

  out <- list(
    theta_final = best_fit$theta_final,
    t.basis.obj = best_fit$t.basis.obj,
    z2.basis.obj = best_fit$z2.basis.obj,
    degree_of_freedom = best_fit$degree_of_freedom,
    lambda = best_fit$lambda.min,
    cv_result = cv_result
  )

  class(out) <- "mrl_fit"
  out
}

predict_mrl <- function(object, newdata) {
  if (!inherits(object, "mrl_fit")) {
    stop("object must be returned by fit_mrl().")
  }

  required_columns <- c("time", "z2")
  missing_columns <- setdiff(required_columns, names(newdata))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in newdata: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  newt.basis <- predict(
    object$t.basis.obj,
    newdata$time
  )

  newz2.basis <- predict(
    object$z2.basis.obj,
    newdata$z2
  )

  newx.basis <- cbind(newt.basis, newz2.basis)
  newX_ns <- cbind(
    newx.basis,
    cproduct(newt.basis, newz2.basis)
  )

  as.vector(
    cbind(1, newX_ns) %*% object$theta_final
  )
}
