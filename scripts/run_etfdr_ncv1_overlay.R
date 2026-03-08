source("R/utils.R")
source("configs/etfdr_ncv1_overlay_config.R")
source("R/etfdr_ncv1_overlay.R")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- arg_flag(args, "smoke", default = FALSE)

  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  suppressPackageStartupMessages(library(glmnet))

  cfg <- get_etfdr_ncv1_overlay_config(smoke = smoke)
  out <- run_ncv1_overlay(cfg)

  message("[done] outputs:")
  message(" - ", out$sel_ncv)
  message(" - ", out$sel_all)
  message(" - ", out$fig)
}

main()
