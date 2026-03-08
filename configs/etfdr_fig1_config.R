get_etfdr_fig1_config <- function(smoke = FALSE) {
  cfg <- list(
    seed = 20260309L,
    n = 600L,
    d = 200L,
    d1 = 20L,
    beta_signal = 1.0,
    sigma = 1.0,
    rho = 0.8,
    folds = 5L,
    lambda_len = 80L,
    fdr_c = 0.1,
    # CVC settings for lambda selection comparison
    cvc_alpha = 0.05,
    cvc_alpha_prime = 0.005,
    cvc_B = 120L,
    reps_main = 200L,
    reps_cvc_compare = 80L,
    save_prefix = "etfdr_fig1"
  )

  if (isTRUE(smoke)) {
    cfg$reps_main <- 50L
    cfg$reps_cvc_compare <- 20L
    cfg$cvc_B <- 80L
    cfg$save_prefix <- "etfdr_fig1_smoke"
  }

  cfg
}
