get_sim1_ncv_config <- function(smoke = FALSE) {
  cfg <- list(
    seed = 20260309L,
    outer_folds = 5L,
    inner_folds = 4L,
    cvc_folds = 5L,
    cvc_alpha = 0.05,
    cvc_alpha_prime = 0.005,
    cvc_B = 200L,
    sigma2_grid = c(1, 4),
    n_grid = c(40L, 80L, 160L, 320L, 640L),
    repetitions = 100L,
    beta_list = list(
      model1 = c(2, 0, 0, 4, 0),
      model2 = c(2, 9, 0, 4, 8)
    ),
    noise = "gaussian",
    save_prefix = "sim1_ncv"
  )

  if (isTRUE(smoke)) {
    cfg$n_grid <- c(40L, 80L)
    cfg$repetitions <- 10L
    cfg$cvc_B <- 80L
    cfg$save_prefix <- "sim1_ncv_smoke"
  }

  cfg
}
