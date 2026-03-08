compute_cvc_pvalues <- function(losses, alpha_prime = 0.005, B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))

  n <- nrow(losses)
  M <- ncol(losses)
  eps <- 1e-8

  pvals <- numeric(M)

  for (m in seq_len(M)) {
    others <- setdiff(seq_len(M), m)
    D <- matrix(losses[, m], nrow = n, ncol = length(others)) - losses[, others, drop = FALSE]
    D <- as.matrix(D)
    mu_hat <- colMeans(D)
    sd_hat <- apply(D, 2, sd)
    sd_hat[sd_hat < eps] <- eps

    T_obs_vec <- sqrt(n) * mu_hat / sd_hat
    T_obs <- max(T_obs_vec)

    # Inequality screening: keep competitors not clearly dominated by m
    z_thr <- qnorm(1 - alpha_prime)
    keep <- which(T_obs_vec > -z_thr)
    if (length(keep) == 0L) {
      pvals[m] <- 1
      next
    }

    D_keep <- D[, keep, drop = FALSE]
    mu_keep <- mu_hat[keep]
    sd_keep <- sd_hat[keep]

    centered <- sweep(D_keep, 2, mu_keep, FUN = "-")
    T_boot <- numeric(B)

    for (b in seq_len(B)) {
      e <- rnorm(n)
      e <- e - mean(e)
      z <- as.numeric((crossprod(centered, e) / n) * sqrt(n) / sd_keep)
      T_boot[b] <- max(z)
    }

    pvals[m] <- (1 + sum(T_boot >= T_obs)) / (B + 1)
  }

  pvals
}

select_cvc_variants <- function(losses, models, alpha, alpha_prime = 0.005, B = 200, seed = NULL) {
  cv_obj <- select_cv_model(losses, models)
  m_cv <- cv_obj$selected

  pvals <- compute_cvc_pvalues(losses, alpha_prime = alpha_prime, B = B, seed = seed)
  sizes <- vapply(models, model_size, integer(1))
  cv_risk <- cv_obj$cv_risk

  Acv <- which(pvals >= alpha)
  if (length(Acv) == 0L) {
    Acv <- m_cv
  }
  min_size <- min(sizes[Acv])
  min_size_models <- Acv[sizes[Acv] == min_size]

  # Variant 1: choose one model among minimal-size set by highest p-value.
  ord_pmax <- order(-pvals[min_size_models], cv_risk[min_size_models], min_size_models)
  selected_pmax <- min_size_models[ord_pmax[1]]

  list(
    pvals = pvals,
    m_cv = m_cv,
    acv = Acv,
    acv_size = length(Acv),
    min_size_models = min_size_models,
    selected_pmax = selected_pmax
  )
}

summarise_sim1_results <- function(raw_df) {
  # cvc variants
  cvc_df <- raw_df[raw_df$method %in% c("cvc_pmax", "cvc_frac"), , drop = FALSE]
  cvc_summary <- aggregate(
    cbind(correct, acv_size) ~ model_name + sigma2 + n + method + alpha,
    data = cvc_df,
    FUN = mean
  )
  names(cvc_summary)[names(cvc_summary) == "correct"] <- "selection_rate"
  names(cvc_summary)[names(cvc_summary) == "acv_size"] <- "mean_acv_size"

  # cv, bic, fy
  other_df <- raw_df[raw_df$method %in% c("cv", "bic", "fy"), , drop = FALSE]
  other_summary <- aggregate(
    correct ~ model_name + sigma2 + n + method,
    data = other_df,
    FUN = mean
  )
  names(other_summary)[names(other_summary) == "correct"] <- "selection_rate"

  list(cvc = cvc_summary, other = other_summary)
}
