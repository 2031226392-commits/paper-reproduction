make_design_matrix <- function(n, d, scenario = c("independent", "correlated"), rho = 0.8) {
  scenario <- match.arg(scenario)
  if (scenario == "independent") {
    return(matrix(rnorm(n * d), nrow = n, ncol = d))
  }
  idx <- 1:d
  Sigma <- rho ^ abs(outer(idx, idx, "-"))
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  Z %*% L
}

generate_fig1_data <- function(cfg, scenario) {
  X <- make_design_matrix(cfg$n, cfg$d, scenario = scenario, rho = cfg$rho)
  beta <- c(rep(cfg$beta_signal, cfg$d1), rep(0, cfg$d - cfg$d1))
  y <- as.numeric(X %*% beta + rnorm(cfg$n, sd = cfg$sigma))
  list(X = X, y = y, beta = beta)
}

compute_ols_pvalues <- function(X, y) {
  n <- nrow(X)
  X1 <- cbind(1, X)
  xtx_inv <- solve(crossprod(X1))
  beta_hat <- as.numeric(xtx_inv %*% crossprod(X1, y))
  resid <- y - as.numeric(X1 %*% beta_hat)
  df <- n - ncol(X1)
  sigma2_hat <- sum(resid ^ 2) / max(df, 1)
  se <- sqrt(diag(xtx_inv) * sigma2_hat)
  tval <- beta_hat / pmax(se, 1e-12)
  2 * pt(-abs(tval[-1]), df = max(df, 1))
}

build_lambda_grid <- function(cfg) {
  set.seed(cfg$seed)
  pilot <- generate_fig1_data(cfg, "independent")
  fit <- glmnet::glmnet(pilot$X, pilot$y, standardize = FALSE, intercept = TRUE)
  lmax <- max(fit$lambda)
  lmin <- lmax * 1e-3
  exp(seq(log(lmax), log(lmin), length.out = cfg$lambda_len))
}

compute_curve_metrics <- function(fit_full, lambda_grid, d1, pvals, c_thresh) {
  coef_mat <- as.matrix(stats::coef(fit_full, s = lambda_grid))
  M <- ncol(coef_mat)
  true_fdr <- numeric(M)
  est_fdr <- numeric(M)

  for (k in seq_len(M)) {
    sel <- which(abs(coef_mat[-1, k]) > 0)
    R <- length(sel)
    if (R == 0L) {
      true_fdr[k] <- 0
      est_fdr[k] <- 0
      next
    }
    false_n <- sum(sel > d1)
    true_fdr[k] <- false_n / R
    est_fdr[k] <- min(1, sum(pvals[sel] > c_thresh) / ((1 - c_thresh) * R))
  }

  list(true_fdr = true_fdr, est_fdr = est_fdr)
}

compute_oof_losses_lambda <- function(X, y, lambda_grid, folds = 5L, seed = NULL) {
  n <- nrow(X)
  if (!is.null(seed)) set.seed(as.integer(seed))
  fold_id <- sample(rep(seq_len(folds), length.out = n))
  M <- length(lambda_grid)
  losses <- matrix(NA_real_, nrow = n, ncol = M)

  for (v in seq_len(folds)) {
    te <- which(fold_id == v)
    tr <- which(fold_id != v)
    fit <- glmnet::glmnet(X[tr, , drop = FALSE], y[tr],
                          lambda = lambda_grid, standardize = FALSE, intercept = TRUE)
    pred <- as.matrix(stats::predict(fit, newx = X[te, , drop = FALSE], s = lambda_grid))
    losses[te, ] <- (matrix(y[te], nrow = length(te), ncol = M) - pred) ^ 2
  }

  losses
}

run_etfdr_scenario <- function(cfg, scenario, lambda_grid, reps, with_cvc = FALSE) {
  M <- length(lambda_grid)
  cvm_mat <- matrix(NA_real_, nrow = reps, ncol = M)
  true_fdr_mat <- matrix(NA_real_, nrow = reps, ncol = M)
  est_fdr_mat <- matrix(NA_real_, nrow = reps, ncol = M)

  sel_records <- list()

  for (r in seq_len(reps)) {
    set.seed(cfg$seed + r + ifelse(scenario == "correlated", 100000L, 0L))
    dat <- generate_fig1_data(cfg, scenario)

    fit_full <- glmnet::glmnet(dat$X, dat$y, lambda = lambda_grid,
                               standardize = FALSE, intercept = TRUE)
    pvals <- compute_ols_pvalues(dat$X, dat$y)
    fdr_obj <- compute_curve_metrics(fit_full, lambda_grid, cfg$d1, pvals, cfg$fdr_c)

    cvfit <- glmnet::cv.glmnet(dat$X, dat$y,
                               lambda = lambda_grid,
                               nfolds = cfg$folds,
                               standardize = FALSE,
                               intercept = TRUE)

    cvm <- as.numeric(cvfit$cvm)
    cvm_mat[r, ] <- cvm
    true_fdr_mat[r, ] <- fdr_obj$true_fdr
    est_fdr_mat[r, ] <- fdr_obj$est_fdr

    idx_min <- which.min(cvm)
    idx_1se <- min(which(cvm <= cvm[idx_min] + cvfit$cvsd[idx_min]))

    rec <- data.frame(
      scenario = scenario,
      rep = r,
      method = c("cv_min", "cv_1se"),
      idx = c(idx_min, idx_1se),
      lambda = c(lambda_grid[idx_min], lambda_grid[idx_1se]),
      cv_mse = c(cvm[idx_min], cvm[idx_1se]),
      true_fdr = c(fdr_obj$true_fdr[idx_min], fdr_obj$true_fdr[idx_1se]),
      est_fdr = c(fdr_obj$est_fdr[idx_min], fdr_obj$est_fdr[idx_1se]),
      stringsAsFactors = FALSE
    )

    if (isTRUE(with_cvc)) {
      losses <- compute_oof_losses_lambda(dat$X, dat$y, lambda_grid,
                                          folds = cfg$folds,
                                          seed = cfg$seed + 500000L + r)
      p_cvc <- compute_cvc_pvalues(losses,
                                   alpha_prime = cfg$cvc_alpha_prime,
                                   B = cfg$cvc_B,
                                   seed = cfg$seed + 600000L + r)
      Acv <- which(p_cvc >= cfg$cvc_alpha)
      if (length(Acv) == 0L) {
        Acv <- idx_min
      }
      # more regularized model = larger lambda
      idx_cvc <- Acv[which.max(lambda_grid[Acv])]
      rec <- rbind(
        rec,
        data.frame(
          scenario = scenario,
          rep = r,
          method = "cvc_max_lambda",
          idx = idx_cvc,
          lambda = lambda_grid[idx_cvc],
          cv_mse = cvm[idx_cvc],
          true_fdr = fdr_obj$true_fdr[idx_cvc],
          est_fdr = fdr_obj$est_fdr[idx_cvc],
          stringsAsFactors = FALSE
        )
      )
    }

    sel_records[[r]] <- rec

    if (r %% 10L == 0L || r == reps) {
      message(sprintf("[etfdr][%s] %d/%d", scenario, r, reps))
    }
  }

  list(
    scenario = scenario,
    lambda = lambda_grid,
    cvm_mean = colMeans(cvm_mat),
    cvm_sd = apply(cvm_mat, 2, sd),
    true_fdr_mean = colMeans(true_fdr_mat),
    true_fdr_sd = apply(true_fdr_mat, 2, sd),
    est_fdr_mean = colMeans(est_fdr_mat),
    est_fdr_sd = apply(est_fdr_mat, 2, sd),
    selected = do.call(rbind, sel_records)
  )
}

plot_etfdr_fig1 <- function(res_list, out_path) {
  png(out_path, width = 1500, height = 700, res = 120)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 4))

  for (res in res_list) {
    x <- -log(res$lambda)
    y1 <- res$cvm_mean

    plot(x, y1, type = "l", lwd = 2, col = "#1F77B4",
         xlab = "-log(lambda)", ylab = "CV MSE",
         main = ifelse(res$scenario == "independent",
                       "(a) Independent X",
                       "(b) Correlated X (AR(0.8))"))

    idx_min <- which.min(res$cvm_mean)
    idx_1se <- min(which(res$cvm_mean <= res$cvm_mean[idx_min] + res$cvm_sd[idx_min]))
    abline(v = x[idx_min], col = "#1F77B4", lwd = 2)
    abline(v = x[idx_1se], col = "#1F77B4", lwd = 2, lty = 2)

    par(new = TRUE)
    plot(x, res$true_fdr_mean, type = "l", lwd = 2, col = "#111111",
         axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    lines(x, res$est_fdr_mean, lwd = 2, col = "#D62728")
    axis(4)
    mtext("False Discovery Rate", side = 4, line = 2.5)

    legend("topleft",
           legend = c("CV MSE", "True FDR", "FDR Est.", "CV min", "1se"),
           col = c("#1F77B4", "#111111", "#D62728", "#1F77B4", "#1F77B4"),
           lty = c(1, 1, 1, 1, 2), lwd = c(2, 2, 2, 2, 2),
           bg = "white", cex = 0.85)
  }

  dev.off()
}

plot_cv_cvc_compare <- function(sel_df, out_path) {
  # compare cv_min vs cvc_max_lambda on CV-MSE and True FDR
  d <- sel_df[sel_df$method %in% c("cv_min", "cvc_max_lambda"), ]
  if (nrow(d) == 0) return(invisible(NULL))

  agg_mean <- aggregate(cbind(cv_mse, true_fdr) ~ scenario + method, data = d, FUN = mean)
  agg_sd <- aggregate(cbind(cv_mse, true_fdr) ~ scenario + method, data = d, FUN = sd)
  names(agg_mean)[3:4] <- c("cv_mse_mean", "true_fdr_mean")
  names(agg_sd)[3:4] <- c("cv_mse_sd", "true_fdr_sd")
  a <- merge(agg_mean, agg_sd, by = c("scenario", "method"))

  png(out_path, width = 1200, height = 700, res = 120)
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

  for (sc in c("independent", "correlated")) {
    cur <- a[a$scenario == sc, ]
    cur <- cur[match(c("cv_min", "cvc_max_lambda"), cur$method), ]

    mids <- barplot(cur$cv_mse_mean,
                    names.arg = c("CV", "CVC"),
                    col = c("#1F77B4", "#FF7F0E"),
                    ylim = c(0, max(cur$cv_mse_mean + cur$cv_mse_sd) * 1.25),
                    ylab = "CV MSE",
                    main = paste0(ifelse(sc == "independent", "Independent", "Correlated"), " (MSE/FDR)"))
    arrows(mids, cur$cv_mse_mean - cur$cv_mse_sd,
           mids, cur$cv_mse_mean + cur$cv_mse_sd,
           angle = 90, code = 3, length = 0.05)

    par(new = TRUE)
    plot(mids, cur$true_fdr_mean, type = "b", pch = 16, lwd = 2, col = "#D62728",
         axes = FALSE, xlab = "", ylab = "",
         ylim = c(0, max(cur$true_fdr_mean + cur$true_fdr_sd, 0.5)))
    arrows(mids, pmax(0, cur$true_fdr_mean - cur$true_fdr_sd),
           mids, cur$true_fdr_mean + cur$true_fdr_sd,
           angle = 90, code = 3, length = 0.05, col = "#D62728")
    axis(4)
    mtext("True FDR", side = 4, line = 2)
    legend("topleft", legend = c("MSE bar", "FDR line"),
           fill = c("#1F77B4", NA), border = c(NA, NA),
           lty = c(NA, 1), pch = c(NA, 16), col = c(NA, "#D62728"),
           bty = "n", cex = 0.85)
  }

  dev.off()
}

run_etfdr_fig1 <- function(cfg) {
  safe_dir_create("results")
  safe_dir_create("figures")

  lambda_grid <- build_lambda_grid(cfg)

  # Main Figure 1 reproduction
  res_ind <- run_etfdr_scenario(cfg, "independent", lambda_grid, reps = cfg$reps_main, with_cvc = FALSE)
  res_cor <- run_etfdr_scenario(cfg, "correlated", lambda_grid, reps = cfg$reps_main, with_cvc = FALSE)

  curve_df <- rbind(
    data.frame(scenario = res_ind$scenario, lambda = res_ind$lambda,
               cv_mse_mean = res_ind$cvm_mean, cv_mse_sd = res_ind$cvm_sd,
               true_fdr_mean = res_ind$true_fdr_mean, true_fdr_sd = res_ind$true_fdr_sd,
               est_fdr_mean = res_ind$est_fdr_mean, est_fdr_sd = res_ind$est_fdr_sd),
    data.frame(scenario = res_cor$scenario, lambda = res_cor$lambda,
               cv_mse_mean = res_cor$cvm_mean, cv_mse_sd = res_cor$cvm_sd,
               true_fdr_mean = res_cor$true_fdr_mean, true_fdr_sd = res_cor$true_fdr_sd,
               est_fdr_mean = res_cor$est_fdr_mean, est_fdr_sd = res_cor$est_fdr_sd)
  )

  curve_path <- file.path("results", paste0(cfg$save_prefix, "_curves.csv"))
  write.csv(curve_df, curve_path, row.names = FALSE)

  fig1_path <- file.path("figures", paste0(cfg$save_prefix, "_figure1_repro.png"))
  plot_etfdr_fig1(list(res_ind, res_cor), fig1_path)

  # CV vs CVC FDR-MSE comparison in similar settings
  cmp_ind <- run_etfdr_scenario(cfg, "independent", lambda_grid, reps = cfg$reps_cvc_compare, with_cvc = TRUE)
  cmp_cor <- run_etfdr_scenario(cfg, "correlated", lambda_grid, reps = cfg$reps_cvc_compare, with_cvc = TRUE)
  sel_df <- rbind(cmp_ind$selected, cmp_cor$selected)

  sel_path <- file.path("results", paste0(cfg$save_prefix, "_cv_cvc_selected.csv"))
  write.csv(sel_df, sel_path, row.names = FALSE)

  cmp_fig_path <- file.path("figures", paste0(cfg$save_prefix, "_cv_cvc_fdr_mse.png"))
  plot_cv_cvc_compare(sel_df, cmp_fig_path)

  list(
    curve_path = curve_path,
    sel_path = sel_path,
    fig1_path = fig1_path,
    cmp_fig_path = cmp_fig_path
  )
}
