make_folds <- function(n, v = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  idx <- sample.int(n)
  fold_id <- integer(n)
  for (i in seq_along(idx)) {
    fold_id[idx[i]] <- ((i - 1L) %% v) + 1L
  }
  fold_id
}

generate_sim1_data <- function(n, beta, sigma2 = 1, noise = "gaussian") {
  p <- length(beta)
  x <- matrix(rnorm(n * (p - 1L)), nrow = n, ncol = p - 1L)
  X <- cbind(1, x)

  eps <- switch(
    noise,
    gaussian = rnorm(n, sd = sqrt(sigma2)),
    t3 = rt(n, df = 3) * sqrt(sigma2 / 3),
    stop("Unsupported noise type: ", noise)
  )

  y <- as.numeric(X %*% beta + eps)
  list(X = X, y = y)
}
