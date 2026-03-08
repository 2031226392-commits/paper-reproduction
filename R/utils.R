set_project_seed <- function(seed) {
  set.seed(as.integer(seed))
}

enumerate_subset_models <- function() {
  # Intercept is always included; subset over x2..x5 (4 predictors -> 16 models)
  base_idx <- 2:5
  models <- vector("list", length = 16)
  k <- 1L
  for (mask in 0:15) {
    bits <- as.logical(intToBits(mask))[1:4]
    active <- c(1L, base_idx[bits])
    models[[k]] <- sort(active)
    k <- k + 1L
  }
  models
}

model_label <- function(active_idx) {
  paste0("x", paste(active_idx, collapse = "_"))
}

model_size <- function(active_idx) {
  length(active_idx)
}

find_true_model_id <- function(models, beta) {
  target <- which(beta != 0)
  hit <- which(vapply(models, function(m) identical(m, target), logical(1)))
  if (length(hit) != 1L) {
    stop("Cannot uniquely identify true model in candidate list.")
  }
  hit
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

arg_flag <- function(args, name, default = FALSE) {
  key <- paste0("--", name)
  if (!any(startsWith(args, key))) return(default)
  value <- sub(paste0("^", key, "=?"), "", args[startsWith(args, key)][1])
  if (identical(value, key) || nchar(value) == 0) return(TRUE)
  tolower(value) %in% c("1", "true", "yes", "y")
}

arg_value <- function(args, name, default = NULL) {
  key <- paste0("--", name, "=")
  idx <- which(startsWith(args, key))
  if (length(idx) == 0) return(default)
  sub(key, "", args[idx[1]])
}
