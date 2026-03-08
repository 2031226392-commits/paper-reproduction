source("configs/sim1_ncv_config.R")
source("R/ncv_sim1.R")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- arg_flag(args, "smoke", default = FALSE)

  cfg <- get_sim1_ncv_config(smoke = smoke)
  out <- run_sim1_ncv(cfg)

  message("[done] outputs:")
  message(" - ", out$raw_path)
  message(" - ", out$summary_path)
  message(" - ", out$fig_path)
}

main()
