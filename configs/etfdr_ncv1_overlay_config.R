get_etfdr_ncv1_overlay_config <- function(smoke = FALSE) {
  cfg <- list(
    seed = 20260310L,
    n = 600L,
    d = 200L,
    d1 = 20L,
    beta_signal = 1.0,
    sigma = 1.0,
    rho = 0.8,
    outer_folds = 5L,
    inner_folds = 4L,
    lambda_len = 80L,
    reps_ncv = 80L,
    base_prefix = "etfdr_fig1",
    save_prefix = "etfdr_fig1_ncv1_overlay"
  )

  if (isTRUE(smoke)) {
    cfg$reps_ncv <- 20L
    cfg$base_prefix <- "etfdr_fig1_smoke"
    cfg$save_prefix <- "etfdr_fig1_smoke_ncv1_overlay"
  }

  cfg
}
