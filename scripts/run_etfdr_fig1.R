source("R/utils.R")
source("R/evaluate_model.R")
source("R/etfdr_fig1.R")
source("configs/etfdr_fig1_config.R")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- arg_flag(args, "smoke", default = FALSE)

  # Ensure user library (where glmnet was installed) is visible.
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  suppressPackageStartupMessages(library(glmnet))

  cfg <- get_etfdr_fig1_config(smoke = smoke)
  out <- run_etfdr_fig1(cfg)

  message("[done] outputs:")
  message(" - ", out$curve_path)
  message(" - ", out$sel_path)
  message(" - ", out$fig1_path)
  message(" - ", out$cmp_fig_path)
}

main()
