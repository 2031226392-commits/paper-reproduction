source("R/utils.R")
source("R/preprocess_data.R")
source("R/fit_model.R")
source("R/evaluate_model.R")
source("R/make_figures.R")
source("configs/sim1_config.R")

run_simulation1 <- function(cfg) {
  safe_dir_create("results")
  safe_dir_create("figures")

  models <- enumerate_subset_models()
  model_sizes <- vapply(models, model_size, integer(1))

  rows <- list()
  row_id <- 1L

  t0 <- proc.time()[3]
  combo_total <- length(cfg$beta_list) * length(cfg$sigma2_grid) * length(cfg$n_grid) * cfg$repetitions
  combo_count <- 0L

  for (model_name in names(cfg$beta_list)) {
    beta <- cfg$beta_list[[model_name]]
    true_model_id <- find_true_model_id(models, beta)

    for (sigma2 in cfg$sigma2_grid) {
      for (n in cfg$n_grid) {
        for (rep_id in seq_len(cfg$repetitions)) {
          combo_count <- combo_count + 1L
          seed_rep <- cfg$seed + combo_count
          set.seed(seed_rep)

          dat <- generate_sim1_data(n = n, beta = beta, sigma2 = sigma2, noise = cfg$noise)
          fold_id <- make_folds(n, v = cfg$folds, seed = seed_rep + 10000L)

          losses <- compute_oof_loss_matrix(dat$X, dat$y, models, fold_id)
          cv_obj <- select_cv_model(losses, models)
          bic_id <- select_bic_model(dat$X, dat$y, models)
          fy_obj <- compute_fy_confset(dat$X, dat$y, models, alpha = cfg$alpha)
          cvc_obj <- select_cvc_variants(
            losses = losses,
            models = models,
            alpha = cfg$alpha,
            alpha_prime = cfg$alpha_prime,
            B = cfg$B,
            seed = seed_rep + 20000L
          )

          rows[[row_id]] <- data.frame(
            method = "cv",
            alpha = NA_real_,
            model_name = model_name,
            sigma2 = sigma2,
            n = n,
            rep = rep_id,
            selected_model_id = cv_obj$selected,
            selected_size = model_sizes[cv_obj$selected],
            true_model_id = true_model_id,
            correct = as.numeric(cv_obj$selected == true_model_id),
            acv_size = NA_real_,
            tie_size_min = NA_real_,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L

          rows[[row_id]] <- data.frame(
            method = "bic",
            alpha = NA_real_,
            model_name = model_name,
            sigma2 = sigma2,
            n = n,
            rep = rep_id,
            selected_model_id = bic_id,
            selected_size = model_sizes[bic_id],
            true_model_id = true_model_id,
            correct = as.numeric(bic_id == true_model_id),
            acv_size = NA_real_,
            tie_size_min = NA_real_,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L

          rows[[row_id]] <- data.frame(
            method = "fy",
            alpha = cfg$alpha,
            model_name = model_name,
            sigma2 = sigma2,
            n = n,
            rep = rep_id,
            selected_model_id = fy_obj$selected,
            selected_size = model_sizes[fy_obj$selected],
            true_model_id = true_model_id,
            correct = as.numeric(fy_obj$selected == true_model_id),
            acv_size = fy_obj$set_size,
            tie_size_min = NA_real_,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L

          rows[[row_id]] <- data.frame(
            method = "cvc_pmax",
            alpha = cfg$alpha,
            model_name = model_name,
            sigma2 = sigma2,
            n = n,
            rep = rep_id,
            selected_model_id = cvc_obj$selected_pmax,
            selected_size = model_sizes[cvc_obj$selected_pmax],
            true_model_id = true_model_id,
            correct = as.numeric(cvc_obj$selected_pmax == true_model_id),
            acv_size = cvc_obj$acv_size,
            tie_size_min = length(cvc_obj$min_size_models),
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L

          frac_correct <- mean(cvc_obj$min_size_models == true_model_id)
          rows[[row_id]] <- data.frame(
            method = "cvc_frac",
            alpha = cfg$alpha,
            model_name = model_name,
            sigma2 = sigma2,
            n = n,
            rep = rep_id,
            selected_model_id = NA_integer_,
            selected_size = min(model_sizes[cvc_obj$min_size_models]),
            true_model_id = true_model_id,
            correct = frac_correct,
            acv_size = cvc_obj$acv_size,
            tie_size_min = length(cvc_obj$min_size_models),
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L

          if (combo_count %% 20L == 0L || combo_count == combo_total) {
            elapsed <- round(proc.time()[3] - t0, 1)
            message(sprintf("[sim1] progress %d/%d, elapsed %.1fs", combo_count, combo_total, elapsed))
          }
        }
      }
    }
  }

  raw_df <- do.call(rbind, rows)
  raw_path <- file.path("results", paste0(cfg$save_prefix, "_raw.csv"))
  write.csv(raw_df, raw_path, row.names = FALSE)

  sm <- summarise_sim1_results(raw_df)
  cvc_path <- file.path("results", paste0(cfg$save_prefix, "_cvc_summary.csv"))
  other_path <- file.path("results", paste0(cfg$save_prefix, "_other_summary.csv"))
  write.csv(sm$cvc, cvc_path, row.names = FALSE)
  write.csv(sm$other, other_path, row.names = FALSE)

  tab1 <- make_table1_compare(sm$cvc, alpha = cfg$alpha, cvc_method = "cvc_pmax")
  tab1_path <- file.path("results", paste0(cfg$save_prefix, "_table1_compare.csv"))
  write.csv(tab1, tab1_path, row.names = FALSE)

  fig_path <- file.path("figures", paste0(cfg$save_prefix, "_selection_rate.png"))
  plot_sim1_selection(sm$cvc, sm$other, out_path = fig_path, alpha = cfg$alpha)
  nonrej_fig_path <- file.path("figures", paste0(cfg$save_prefix, "_nonrejected_count.png"))
  plot_sim1_nonrejected(raw_df, out_path = nonrej_fig_path, alpha = cfg$alpha)

  list(
    raw_path = raw_path,
    cvc_summary_path = cvc_path,
    other_summary_path = other_path,
    table1_path = tab1_path,
    fig_path = fig_path,
    nonrejected_fig_path = nonrej_fig_path
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  exp_name <- arg_value(args, "exp", default = "sim1")
  smoke <- arg_flag(args, "smoke", default = FALSE)

  if (!identical(exp_name, "sim1")) {
    stop("Currently only --exp=sim1 is implemented.")
  }

  cfg <- get_sim1_config(smoke = smoke)
  out <- run_simulation1(cfg)

  message("[done] outputs:")
  message(" - ", out$raw_path)
  message(" - ", out$cvc_summary_path)
  message(" - ", out$other_summary_path)
  message(" - ", out$table1_path)
  message(" - ", out$fig_path)
  message(" - ", out$nonrejected_fig_path)
}

main()
