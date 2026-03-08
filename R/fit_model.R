fit_lm_subset <- function(X, y, active_idx) {
  x_sub <- X[, active_idx, drop = FALSE]
  fit <- lm.fit(x = x_sub, y = y)
  list(coefficients = fit$coefficients, active_idx = active_idx)
}

predict_lm_subset <- function(fit, X) {
  x_sub <- X[, fit$active_idx, drop = FALSE]
  as.numeric(x_sub %*% fit$coefficients)
}

compute_oof_loss_matrix <- function(X, y, models, fold_id) {
  n <- nrow(X)
  M <- length(models)
  losses <- matrix(NA_real_, nrow = n, ncol = M)

  folds <- sort(unique(fold_id))
  for (v in folds) {
    te <- which(fold_id == v)
    tr <- which(fold_id != v)

    X_tr <- X[tr, , drop = FALSE]
    y_tr <- y[tr]
    X_te <- X[te, , drop = FALSE]
    y_te <- y[te]

    for (m in seq_len(M)) {
      fit <- fit_lm_subset(X_tr, y_tr, models[[m]])
      pred <- predict_lm_subset(fit, X_te)
      losses[te, m] <- (y_te - pred) ^ 2
    }
  }

  colnames(losses) <- vapply(models, model_label, character(1))
  losses
}

select_cv_model <- function(losses, models) {
  cv_risk <- colMeans(losses)
  sizes <- vapply(models, model_size, integer(1))

  ord <- order(cv_risk, sizes, seq_along(cv_risk))
  m_cv <- ord[1]

  list(
    selected = m_cv,
    cv_risk = cv_risk
  )
}

select_bic_model <- function(X, y, models) {
  n <- nrow(X)
  bic <- numeric(length(models))
  for (m in seq_along(models)) {
    fit <- fit_lm_subset(X, y, models[[m]])
    pred <- predict_lm_subset(fit, X)
    rss <- sum((y - pred) ^ 2)
    k <- length(models[[m]])
    bic[m] <- n * log(rss / n + 1e-12) + log(n) * k
  }
  which.min(bic)
}

compute_fy_confset <- function(X, y, models, alpha = 0.05) {
  n <- nrow(X)
  M <- length(models)
  full_id <- which.max(vapply(models, length, integer(1)))
  full_fit <- fit_lm_subset(X, y, models[[full_id]])
  full_pred <- predict_lm_subset(full_fit, X)
  rss_full <- sum((y - full_pred) ^ 2)
  k_full <- length(models[[full_id]])
  df2 <- max(1, n - k_full)

  pvals <- rep(1, M)
  keep <- logical(M)
  rss_all <- numeric(M)

  for (m in seq_len(M)) {
    fit_m <- fit_lm_subset(X, y, models[[m]])
    pred_m <- predict_lm_subset(fit_m, X)
    rss_m <- sum((y - pred_m) ^ 2)
    rss_all[m] <- rss_m
    k_m <- length(models[[m]])

    if (m == full_id || k_m == k_full) {
      pvals[m] <- 1
      keep[m] <- TRUE
      next
    }

    df1 <- k_full - k_m
    num <- (rss_m - rss_full) / max(1, df1)
    den <- rss_full / df2
    f_stat <- num / max(den, 1e-12)
    pvals[m] <- pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
    keep[m] <- (pvals[m] >= alpha)
  }

  fy_set <- which(keep)
  if (length(fy_set) == 0L) {
    fy_set <- full_id
  }

  sizes <- vapply(models, model_size, integer(1))
  ord <- order(sizes[fy_set], rss_all[fy_set], fy_set)
  selected <- fy_set[ord[1]]

  list(
    selected = selected,
    set = fy_set,
    set_size = length(fy_set),
    pvals = pvals
  )
}
