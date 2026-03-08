source("configs/sim1_ncv_config.R")
source("R/ncv_sim1.R")

append_ncv0 <- function(cfg) {
  raw_path <- file.path("results", paste0(cfg$save_prefix, "_raw.csv"))
  if (!file.exists(raw_path)) stop("Base raw file not found: ", raw_path)

  raw <- read.csv(raw_path)
  # keep non-ncv0 rows if rerun
  raw <- raw[raw$method != "ncv0", ]

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

          ncv <- run_one_ncv(dat$X, dat$y, models, cfg, seed_base = seed_rep + 3000L)
          mid <- ncv$ncv0

          rows[[rid]] <- data.frame(
            method = "ncv0",
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

          if (cnt %% 50L == 0L || cnt == total) {
            message(sprintf("[append_ncv0][%s] %d/%d", cfg$save_prefix, cnt, total))
          }
        }
      }
    }
  }

  new_rows <- do.call(rbind, rows)
  all_raw <- rbind(raw, new_rows)
  all_raw <- all_raw[order(all_raw$model_name, all_raw$sigma2, all_raw$n, all_raw$rep, all_raw$method), ]
  write.csv(all_raw, raw_path, row.names = FALSE)

  sm <- summarise_ncv_results(all_raw)
  sum_path <- file.path("results", paste0(cfg$save_prefix, "_summary.csv"))
  write.csv(sm, sum_path, row.names = FALSE)

  fig_path <- file.path("figures", paste0(cfg$save_prefix, "_selection_rate.png"))
  plot_ncv_selection(sm, fig_path)

  list(raw_path = raw_path, summary_path = sum_path, fig_path = fig_path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- arg_flag(args, "smoke", default = FALSE)
  cfg <- get_sim1_ncv_config(smoke = smoke)
  out <- append_ncv0(cfg)
  message("[done] outputs:")
  message(" - ", out$raw_path)
  message(" - ", out$summary_path)
  message(" - ", out$fig_path)
}

main()
