source("R/etfdr_fig1.R")
source("R/preprocess_data.R")

compute_oof_with_fold <- function(X, y, lambda_grid, folds = 4L, seed = NULL) {
  n <- nrow(X)
  if (!is.null(seed)) set.seed(as.integer(seed))
  fold_id <- sample(rep(seq_len(folds), length.out = n))

  M <- length(lambda_grid)
  losses <- matrix(NA_real_, nrow = n, ncol = M)
  for (v in seq_len(folds)) {
    te <- which(fold_id == v)
    tr <- which(fold_id != v)
    fit <- glmnet::glmnet(X[tr, , drop = FALSE], y[tr],
                          lambda = lambda_grid,
                          standardize = FALSE,
                          intercept = TRUE)
    pred <- as.matrix(stats::predict(fit, newx = X[te, , drop = FALSE], s = lambda_grid))
    losses[te, ] <- (matrix(y[te], nrow = length(te), ncol = M) - pred) ^ 2
  }
  list(losses = losses, fold_id = fold_id)
}

select_inner_min_1se <- function(losses, fold_id) {
  folds <- sort(unique(fold_id))
  fm <- matrix(NA_real_, nrow = length(folds), ncol = ncol(losses))
  for (i in seq_along(folds)) {
    fm[i, ] <- colMeans(losses[fold_id == folds[i], , drop = FALSE])
  }
  mu <- colMeans(fm)
  se <- apply(fm, 2, sd) / sqrt(nrow(fm))
  se[is.na(se)] <- 0

  idx_min <- which.min(mu)
  thr <- mu[idx_min] + se[idx_min]
  cand <- which(mu <= thr)
  # stronger regularization = larger lambda = smaller -log(lambda)
  idx_1se <- cand[which.max(cand)]
  list(idx_min = idx_min, idx_1se = idx_1se)
}

select_ncv01 <- function(X, y, lambda_grid, outer_folds = 5L, inner_folds = 4L, seed = NULL) {
  n <- nrow(X)
  if (!is.null(seed)) set.seed(as.integer(seed))
  outer_id <- sample(rep(seq_len(outer_folds), length.out = n))

  idx_min_vec <- integer(outer_folds)
  idx_1se_vec <- integer(outer_folds)
  err_min <- numeric(outer_folds)
  err_1se <- numeric(outer_folds)

  for (o in seq_len(outer_folds)) {
    te <- which(outer_id == o)
    tr <- which(outer_id != o)

    X_tr <- X[tr, , drop = FALSE]
    y_tr <- y[tr]
    X_te <- X[te, , drop = FALSE]
    y_te <- y[te]

    inner <- compute_oof_with_fold(X_tr, y_tr, lambda_grid,
                                   folds = inner_folds,
                                   seed = seed + 1000L + o)
    sel <- select_inner_min_1se(inner$losses, inner$fold_id)

    fit_tr <- glmnet::glmnet(X_tr, y_tr, lambda = lambda_grid,
                             standardize = FALSE, intercept = TRUE)

    pred_min <- as.numeric(stats::predict(fit_tr, newx = X_te, s = lambda_grid[sel$idx_min]))
    pred_1se <- as.numeric(stats::predict(fit_tr, newx = X_te, s = lambda_grid[sel$idx_1se]))

    err_min[o] <- mean((y_te - pred_min) ^ 2)
    err_1se[o] <- mean((y_te - pred_1se) ^ 2)
    idx_min_vec[o] <- sel$idx_min
    idx_1se_vec[o] <- sel$idx_1se
  }

  se_outer_min <- stats::sd(err_min) / sqrt(length(err_min))
  if (is.na(se_outer_min)) se_outer_min <- 0
  thr_min <- min(err_min) + se_outer_min
  cand_min <- which(err_min <= thr_min)

  se_outer_1se <- stats::sd(err_1se) / sqrt(length(err_1se))
  if (is.na(se_outer_1se)) se_outer_1se <- 0
  thr_1se <- min(err_1se) + se_outer_1se
  cand_1se <- which(err_1se <= thr_1se)

  list(
    idx_ncv0 = idx_min_vec[which.min(err_min)],
    idx_ncv0_1se = idx_min_vec[cand_min[which.max(lambda_grid[idx_min_vec[cand_min]])]],
    idx_ncv1 = idx_1se_vec[which.min(err_1se)],
    idx_ncv1_1se = idx_1se_vec[cand_1se[which.max(lambda_grid[idx_1se_vec[cand_1se]])]]
  )
}

compute_ncv_curve <- function(X, y, lambda_grid, d1, fdr_c,
                              outer_folds = 5L, seed = NULL) {
  n <- nrow(X)
  M <- length(lambda_grid)
  if (!is.null(seed)) set.seed(as.integer(seed))
  outer_id <- sample(rep(seq_len(outer_folds), length.out = n))

  mse_fold <- matrix(NA_real_, nrow = outer_folds, ncol = M)
  true_fdr_fold <- matrix(NA_real_, nrow = outer_folds, ncol = M)
  est_fdr_fold <- matrix(NA_real_, nrow = outer_folds, ncol = M)

  for (o in seq_len(outer_folds)) {
    te <- which(outer_id == o)
    tr <- which(outer_id != o)
    X_tr <- X[tr, , drop = FALSE]
    y_tr <- y[tr]
    X_te <- X[te, , drop = FALSE]
    y_te <- y[te]

    fit_tr <- glmnet::glmnet(X_tr, y_tr, lambda = lambda_grid,
                             standardize = FALSE, intercept = TRUE)
    pred <- as.matrix(stats::predict(fit_tr, newx = X_te, s = lambda_grid))
    mse_fold[o, ] <- colMeans((matrix(y_te, nrow = length(te), ncol = M) - pred) ^ 2)

    coef_mat <- as.matrix(stats::coef(fit_tr, s = lambda_grid))
    pvals_tr <- compute_ols_pvalues(X_tr, y_tr)

    for (k in seq_len(M)) {
      sel <- which(abs(coef_mat[-1, k]) > 0)
      Rk <- length(sel)
      if (Rk == 0L) {
        true_fdr_fold[o, k] <- 0
        est_fdr_fold[o, k] <- 0
      } else {
        true_fdr_fold[o, k] <- sum(sel > d1) / Rk
        est_fdr_fold[o, k] <- min(1, sum(pvals_tr[sel] > fdr_c) / ((1 - fdr_c) * Rk))
      }
    }
  }

  list(
    mse = colMeans(mse_fold),
    true_fdr = colMeans(true_fdr_fold),
    est_fdr = colMeans(est_fdr_fold)
  )
}

build_overlay_curves_from_base <- function(base_curve_df,
                                           methods = c("cv", "cvc", "ncv0", "ncv1")) {
  rows <- list()
  rid <- 1L
  for (scenario in c("independent", "correlated")) {
    base <- base_curve_df[base_curve_df$scenario == scenario, ]
    base <- base[order(base$lambda, decreasing = TRUE), ]
    for (m in methods) {
      rows[[rid]] <- data.frame(
        scenario = scenario,
        method = m,
        lambda = base$lambda,
        mse_mean = base$cv_mse_mean,
        true_fdr_mean = base$true_fdr_mean,
        est_fdr_mean = base$est_fdr_mean,
        fdp_mean = base$true_fdr_mean,
        stringsAsFactors = FALSE
      )
      rid <- rid + 1L
    }
  }
  do.call(rbind, rows)
}

run_ncv1_overlay <- function(cfg) {
  safe_dir_create("results")
  safe_dir_create("figures")

  base_curve_path <- file.path("results", paste0(cfg$base_prefix, "_curves.csv"))
  base_sel_path <- file.path("results", paste0(cfg$base_prefix, "_cv_cvc_selected.csv"))

  if (!file.exists(base_curve_path) || !file.exists(base_sel_path)) {
    stop("Base ETFDR results not found. Please run scripts/run_etfdr_fig1.R first.")
  }

  curve_df <- read.csv(base_curve_path)
  sel_base <- read.csv(base_sel_path)

  lambda_grid <- sort(unique(curve_df$lambda), decreasing = TRUE)
  M <- length(lambda_grid)
  scen <- c("independent", "correlated")
  ncv_curve_store <- list(
    independent = list(
      mse = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M),
      true_fdr = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M),
      est_fdr = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M)
    ),
    correlated = list(
      mse = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M),
      true_fdr = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M),
      est_fdr = matrix(NA_real_, nrow = cfg$reps_ncv, ncol = M)
    )
  )

  rows <- list()
  rid <- 1L
  for (scenario in scen) {
    for (r in seq_len(cfg$reps_ncv)) {
      set.seed(cfg$seed + r + ifelse(scenario == "correlated", 100000L, 0L))
      dat <- generate_fig1_data(cfg, scenario)

      sel <- select_ncv01(dat$X, dat$y, lambda_grid,
                          outer_folds = cfg$outer_folds,
                          inner_folds = cfg$inner_folds,
                          seed = cfg$seed + 200000L + r + ifelse(scenario == "correlated", 100000L, 0L))
      ncv_curve <- compute_ncv_curve(dat$X, dat$y, lambda_grid,
                                     d1 = cfg$d1,
                                     fdr_c = cfg$fdr_c,
                                     outer_folds = cfg$outer_folds,
                                     seed = cfg$seed + 300000L + r + ifelse(scenario == "correlated", 100000L, 0L))
      ncv_curve_store[[scenario]]$mse[r, ] <- ncv_curve$mse
      ncv_curve_store[[scenario]]$true_fdr[r, ] <- ncv_curve$true_fdr
      ncv_curve_store[[scenario]]$est_fdr[r, ] <- ncv_curve$est_fdr

      rows[[rid]] <- data.frame(scenario = scenario, rep = r,
                                method = "ncv0_min",
                                idx = sel$idx_ncv0,
                                lambda = lambda_grid[sel$idx_ncv0],
                                stringsAsFactors = FALSE)
      rid <- rid + 1L

      rows[[rid]] <- data.frame(scenario = scenario, rep = r,
                                method = "ncv0_1se",
                                idx = sel$idx_ncv0_1se,
                                lambda = lambda_grid[sel$idx_ncv0_1se],
                                stringsAsFactors = FALSE)
      rid <- rid + 1L

      rows[[rid]] <- data.frame(scenario = scenario, rep = r,
                                method = "ncv1_min",
                                idx = sel$idx_ncv1,
                                lambda = lambda_grid[sel$idx_ncv1],
                                stringsAsFactors = FALSE)
      rid <- rid + 1L

      rows[[rid]] <- data.frame(scenario = scenario, rep = r,
                                method = "ncv1_1se",
                                idx = sel$idx_ncv1_1se,
                                lambda = lambda_grid[sel$idx_ncv1_1se],
                                stringsAsFactors = FALSE)
      rid <- rid + 1L

      if (r %% 10L == 0L || r == cfg$reps_ncv) {
        message(sprintf("[ncv1_overlay][%s] %d/%d", scenario, r, cfg$reps_ncv))
      }
    }
  }

  sel_ncv <- do.call(rbind, rows)
  out_sel <- file.path("results", paste0(cfg$save_prefix, "_ncv1_selected.csv"))
  write.csv(sel_ncv, out_sel, row.names = FALSE)

  sel_all <- rbind(
    sel_base[, c("scenario", "rep", "method", "idx", "lambda")],
    sel_ncv
  )
  out_all <- file.path("results", paste0(cfg$save_prefix, "_all_selected.csv"))
  write.csv(sel_all, out_all, row.names = FALSE)

  curve_rows <- list()
  cid <- 1L
  for (scenario in scen) {
    base <- curve_df[curve_df$scenario == scenario, ]
    base <- base[order(base$lambda, decreasing = TRUE), ]
    ncv_mse <- colMeans(ncv_curve_store[[scenario]]$mse)
    ncv_true <- colMeans(ncv_curve_store[[scenario]]$true_fdr)
    ncv_est <- colMeans(ncv_curve_store[[scenario]]$est_fdr)

    curve_rows[[cid]] <- data.frame(
      scenario = scenario, method = "cv", lambda = base$lambda,
      mse_mean = base$cv_mse_mean,
      true_fdr_mean = base$true_fdr_mean,
      est_fdr_mean = base$est_fdr_mean,
      fdp_mean = base$true_fdr_mean,
      stringsAsFactors = FALSE
    )
    cid <- cid + 1L

    curve_rows[[cid]] <- data.frame(
      scenario = scenario, method = "cvc", lambda = base$lambda,
      mse_mean = base$cv_mse_mean,
      true_fdr_mean = base$true_fdr_mean,
      est_fdr_mean = base$est_fdr_mean,
      fdp_mean = base$true_fdr_mean,
      stringsAsFactors = FALSE
    )
    cid <- cid + 1L

    curve_rows[[cid]] <- data.frame(
      scenario = scenario, method = "ncv0", lambda = base$lambda,
      mse_mean = ncv_mse,
      true_fdr_mean = ncv_true,
      est_fdr_mean = ncv_est,
      fdp_mean = ncv_true,
      stringsAsFactors = FALSE
    )
    cid <- cid + 1L

    curve_rows[[cid]] <- data.frame(
      scenario = scenario, method = "ncv1", lambda = base$lambda,
      mse_mean = ncv_mse,
      true_fdr_mean = ncv_true,
      est_fdr_mean = ncv_est,
      fdp_mean = ncv_true,
      stringsAsFactors = FALSE
    )
    cid <- cid + 1L
  }

  curve_overlay <- do.call(rbind, curve_rows)
  out_curve <- file.path("results", paste0(cfg$save_prefix, "_curve_overlay.csv"))
  write.csv(curve_overlay, out_curve, row.names = FALSE)

  plot_overlay(curve_overlay, sel_all,
               out_path = file.path("figures", paste0(cfg$save_prefix, "_figure.png")))

  list(sel_ncv = out_sel,
       sel_all = out_all,
       curve_overlay = out_curve,
       fig = file.path("figures", paste0(cfg$save_prefix, "_figure.png")))
}

plot_overlay <- function(curve_df, sel_all, out_path) {
  png(out_path, width = 1700, height = 760, res = 120)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 4))

  if (!("fdp_mean" %in% names(curve_df))) {
    curve_df$fdp_mean <- curve_df$true_fdr_mean
  }

  mse_cols <- c(cv = "#1F77B4", cvc = "#FF7F0E", ncv0 = "#7F7F7F", ncv1 = "#2CA02C")
  fdr_cols <- c(cv = "#1F77B4", cvc = "#FF7F0E", ncv0 = "#7F7F7F", ncv1 = "#2CA02C")
  methods <- c("cv", "cvc", "ncv0", "ncv1")

  for (scenario in c("independent", "correlated")) {
    cur <- curve_df[curve_df$scenario == scenario, ]
    cur <- cur[order(cur$lambda, decreasing = TRUE), ]
    x <- sort(unique(cur$lambda), decreasing = TRUE)

    y_mse_max <- max(cur$mse_mean, na.rm = TRUE)
    plot(x, rep(NA_real_, length(x)),
         type = "n",
         xlab = "lambda", ylab = "MSE",
         log = "x",
         ylim = c(0, y_mse_max * 1.05),
         main = ifelse(scenario == "independent", "Independent", "Correlated (AR(0.8))"))
    for (m in methods) {
      tmp <- cur[cur$method == m, ]
      tmp <- tmp[order(tmp$lambda, decreasing = TRUE), ]
      lines(tmp$lambda, tmp$mse_mean, lwd = 2, col = mse_cols[m], lty = 1)
    }

    sub <- sel_all[sel_all$scenario == scenario, ]
    mean_by_method <- aggregate(lambda ~ method, data = sub, FUN = mean)
    xline <- function(method, col, lty, lwd = 2) {
      val <- mean_by_method$lambda[mean_by_method$method == method]
      if (length(val) == 1L) abline(v = val, col = col, lty = lty, lwd = lwd)
    }

    par(new = TRUE)
    plot(x, rep(NA_real_, length(x)),
         type = "n",
         axes = FALSE, xlab = "", ylab = "",
         log = "x",
         ylim = c(0, 1))
    for (m in methods) {
      tmp <- cur[cur$method == m, ]
      tmp <- tmp[order(tmp$lambda, decreasing = TRUE), ]
      lines(tmp$lambda, tmp$true_fdr_mean, lwd = 1.6, col = fdr_cols[m], lty = 1)
      lines(tmp$lambda, tmp$est_fdr_mean, lwd = 1.6, col = fdr_cols[m], lty = 2)
      lines(tmp$lambda, tmp$fdp_mean, lwd = 1.6, col = fdr_cols[m], lty = 3)
    }
    axis(4)
    mtext("FDR / FDP", side = 4, line = 2)

    xline("cv_min", "#1F77B4", 1)
    xline("cvc_max_lambda", "#FF7F0E", 1)
    xline("ncv0_min", "#7F7F7F", 1)
    xline("ncv1_min", "#2CA02C", 1)
    xline("cv_1se", "#1F77B4", 2)
    xline("ncv0_1se", "#7F7F7F", 2)
    xline("ncv1_1se", "#2CA02C", 2)

    legend("topleft",
           legend = c("MSE: cv/cvc/ncv0/ncv1 (color)",
                      "FDR_true (solid)", "FDR_hat (dashed)", "FDP (dotted)",
                      "cv", "cvc", "ncv0", "ncv1",
                      "cv_1se", "ncv0_1se", "ncv1_1se"),
           col = c("#111111",
                   "#111111", "#111111", "#111111",
                   "#1F77B4", "#FF7F0E", "#7F7F7F", "#2CA02C",
                   "#1F77B4", "#7F7F7F", "#2CA02C"),
           lty = c(NA, 1, 2, 3, 1, 1, 1, 1, 2, 2, 2),
           lwd = c(NA, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
           pch = c(NA, NA, NA, NA, 15, 15, 15, 15, NA, NA, NA),
           pt.cex = c(NA, NA, NA, NA, 0.9, 0.9, 0.9, 0.9, NA, NA, NA),
           cex = 0.78, bg = "white")
  }

  dev.off()
}
