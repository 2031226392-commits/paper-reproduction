get_sim1_config <- function(smoke = FALSE) {
  cfg <- list(
    seed = 20260308L,
    folds = 5L,
    B = 200L,
    alpha_prime = 0.005,
    alpha = 0.05,
    sigma2_grid = c(1, 4),
    n_grid = c(40L, 80L, 160L, 320L, 640L),
    repetitions = 100L,
    beta_list = list(
      model1 = c(2, 0, 0, 4, 0),
      model2 = c(2, 9, 0, 4, 8)
    ),
    noise = "gaussian",
    save_prefix = "sim1_alpha005"
  )

  if (isTRUE(smoke)) {
    cfg$n_grid <- c(40L, 80L)
    cfg$repetitions <- 10L
    cfg$B <- 60L
    cfg$save_prefix <- "sim1_smoke_alpha005"
  }

  cfg
}
