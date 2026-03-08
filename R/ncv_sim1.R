source("R/utils.R")
source("R/preprocess_data.R")
source("R/fit_model.R")
source("R/evaluate_model.R")

compute_fold_means <- function(losses, fold_id) {
  folds <- sort(unique(fold_id))
  out <- matrix(NA_real_, nrow = length(folds), ncol = ncol(losses))
  for (i in seq_along(folds)) {
    idx <- which(fold_id == folds[i])
    out[i, ] <- colMeans(losses[idx, , drop = FALSE])
  }
  out
}

select_cv_1se_model <- function(losses, models, fold_id) {
  fm <- compute_fold_means(losses, fold_id)
  mu <- colMeans(fm)
  se <- apply(fm, 2, sd) / sqrt(nrow(fm))
  se[is.na(se)] <- 0

  m_min <- which.min(mu)
  thr <- mu[m_min] + se[m_min]
  cand <- which(mu <= thr)
  sizes <- vapply(models, model_size, integer(1))
  ord <- order(sizes[cand], mu[cand], cand)
  cand[ord[1]]
}

choose_outer_best <- function(errors, model_ids, models) {
  sizes <- vapply(models, model_size, integer(1))
  ord <- order(errors, sizes[model_ids], model_ids)
  model_ids[ord[1]]
}

choose_outer_1se <- function(errors, model_ids, models) {
  sizes <- vapply(models, model_size, integer(1))
  se_outer <- stats::sd(errors) / sqrt(length(errors))
  if (is.na(se_outer)) se_outer <- 0
  thr <- min(errors) + se_outer
  cand <- which(errors <= thr)
  ord <- order(sizes[model_ids[cand]], errors[cand], model_ids[cand])
  model_ids[cand][ord[1]]
}

run_one_ncv <- function(X, y, models, cfg, seed_base) {
  n <- nrow(X)
  outer_id <- make_folds(n, v = cfg$outer_folds, seed = seed_base + 9000L)

  e_best <- numeric(cfg$outer_folds)
  e_1se <- numeric(cfg$outer_folds)
  m_best <- integer(cfg$outer_folds)
  m_1se <- integer(cfg$outer_folds)

  for (o in seq_len(cfg$outer_folds)) {
    te <- which(outer_id == o)
    tr <- which(outer_id != o)

    X_tr <- X[tr, , drop = FALSE]
    y_tr <- y[tr]
    X_te <- X[te, , drop = FALSE]
    y_te <- y[te]

    inner_id <- make_folds(length(tr), v = cfg$inner_folds, seed = seed_base + 10000L + o)
    inner_losses <- compute_oof_loss_matrix(X_tr, y_tr, models, inner_id)

    m_ib <- select_cv_model(inner_losses, models)$selected
    m_i1 <- select_cv_1se_model(inner_losses, models, inner_id)

    fit_ib <- fit_lm_subset(X_tr, y_tr, models[[m_ib]])
    fit_i1 <- fit_lm_subset(X_tr, y_tr, models[[m_i1]])

    e_best[o] <- mean((y_te - predict_lm_subset(fit_ib, X_te)) ^ 2)
    e_1se[o] <- mean((y_te - predict_lm_subset(fit_i1, X_te)) ^ 2)
    m_best[o] <- m_ib
    m_1se[o] <- m_i1
  }

  list(
    ncv1 = choose_outer_best(e_1se, m_1se, models),
    ncv2 = choose_outer_1se(e_best, m_best, models),
    ncv3 = choose_outer_1se(e_1se, m_1se, models)
  )
}

summarise_ncv_results <- function(raw_df) {
  agg <- aggregate(correct ~ method + model_name + sigma2 + n,
                   data = raw_df,
                   FUN = function(z) c(mean = mean(z), se = stats::sd(z) / sqrt(length(z))))
  data.frame(
    method = agg$method,
    model_name = agg$model_name,
    sigma2 = agg$sigma2,
    n = agg$n,
    selection_rate = agg$correct[, "mean"],
    se = agg$correct[, "se"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

plot_ncv_selection <- function(summary_df, out_path) {
  safe_dir_create(dirname(out_path))
  png(out_path, width = 1400, height = 900, res = 120)
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

  methods <- c("cv", "cv_1se", "cvc", "ncv1", "ncv2", "ncv3")
  cols <- c("#1F77B4", "#17BECF", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD")
  pchs <- c(16, 15, 17, 18, 8, 3)
  ltys <- c(1, 2, 1, 1, 2, 3)

  for (model_name in c("model1", "model2")) {
    for (sigma2 in c(1, 4)) {
      sub <- summary_df[summary_df$model_name == model_name & summary_df$sigma2 == sigma2, ]
      x <- sort(unique(sub$n))
      plot(NA, xlim = range(x), ylim = c(0, 1),
           xlab = "n", ylab = "Correct Selection Rate",
           main = paste0(model_name, ", sigma2=", sigma2))
      grid()

      for (i in seq_along(methods)) {
        m <- methods[i]
        cur <- sub[sub$method == m, ]
        if (nrow(cur) == 0) next
        cur <- cur[order(cur$n), ]
        lines(cur$n, cur$selection_rate, type = "b", lwd = 2,
              col = cols[i], pch = pchs[i], lty = ltys[i])
      }

      legend("bottomright", legend = methods, col = cols,
             lty = ltys, pch = pchs, cex = 0.82, bty = "n")
    }
  }

  dev.off()
}

run_sim1_ncv <- function(cfg) {
  safe_dir_create("results")
  safe_dir_create("figures")

  models <- enumerate_subset_models()
  sizes <- vapply(models, model_size, integer(1))

  rows <- list()
  rid <- 1L

  total <- length(cfg$beta_list) * length(cfg$sigma2_grid) * length(cfg$n_grid) * cfg$repetitions
  cnt <- 0L

  for (model_name in names(cfg$beta_list)) {
    beta <- cfg$beta_list[[model_name]]
    true_id <- find_true_model_id(models, beta)

    for (sigma2 in cfg$sigma2_grid) {
      for (n in cfg$n_grid) {
        for (rep_id in seq_len(cfg$repetitions)) {
          cnt <- cnt + 1L
          seed_rep <- cfg$seed + cnt
          set.seed(seed_rep)

          dat <- generate_sim1_data(n = n, beta = beta, sigma2 = sigma2, noise = cfg$noise)

          # shared 5-fold for cv/cv_1se/cvc
          fold_id <- make_folds(n, v = cfg$cvc_folds, seed = seed_rep + 1000L)
          losses <- compute_oof_loss_matrix(dat$X, dat$y, models, fold_id)

          m_cv <- select_cv_model(losses, models)$selected
          m_cv1 <- select_cv_1se_model(losses, models, fold_id)
          cvc <- select_cvc_variants(losses, models,
                                     alpha = cfg$cvc_alpha,
                                     alpha_prime = cfg$cvc_alpha_prime,
                                     B = cfg$cvc_B,
                                     seed = seed_rep + 2000L)$selected_pmax

          ncv <- run_one_ncv(dat$X, dat$y, models, cfg, seed_base = seed_rep + 3000L)

          for (pair in list(
            list("cv", m_cv),
            list("cv_1se", m_cv1),
            list("cvc", cvc),
            list("ncv1", ncv$ncv1),
            list("ncv2", ncv$ncv2),
            list("ncv3", ncv$ncv3)
          )) {
            method <- pair[[1]]
            mid <- pair[[2]]
            rows[[rid]] <- data.frame(
              method = method,
              model_name = model_name,
              sigma2 = sigma2,
              n = n,
              rep = rep_id,
              selected_model_id = mid,
              selected_size = sizes[mid],
              true_model_id = true_id,
              correct = as.numeric(mid == true_id),
              stringsAsFactors = FALSE
            )
            rid <- rid + 1L
          }

          if (cnt %% 20L == 0L || cnt == total) {
            message(sprintf("[sim1_ncv] %d/%d", cnt, total))
          }
        }
      }
    }
  }

  raw <- do.call(rbind, rows)
  raw_path <- file.path("results", paste0(cfg$save_prefix, "_raw.csv"))
  write.csv(raw, raw_path, row.names = FALSE)

  sm <- summarise_ncv_results(raw)
  sum_path <- file.path("results", paste0(cfg$save_prefix, "_summary.csv"))
  write.csv(sm, sum_path, row.names = FALSE)

  fig_path <- file.path("figures", paste0(cfg$save_prefix, "_selection_rate.png"))
  plot_ncv_selection(sm, fig_path)

  list(raw_path = raw_path, summary_path = sum_path, fig_path = fig_path)
}
