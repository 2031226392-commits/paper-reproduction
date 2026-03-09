source("R/etfdr_ncv1_overlay.R")
source("R/utils.R")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  smoke <- arg_flag(args, "smoke", default = FALSE)

  base_prefix <- if (smoke) "etfdr_fig1_smoke" else "etfdr_fig1"
  overlay_prefix <- if (smoke) "etfdr_fig1_smoke_ncv1_overlay" else "etfdr_fig1_ncv1_overlay"

  curve_path <- file.path("results", paste0(base_prefix, "_curves.csv"))
  sel_path <- file.path("results", paste0(overlay_prefix, "_all_selected.csv"))

  if (!file.exists(curve_path) || !file.exists(sel_path)) {
    stop("Missing existing CSVs. Need both base curves and overlay selections.")
  }

  curve_base <- read.csv(curve_path)
  sel_all <- read.csv(sel_path)

  curve_overlay <- build_overlay_curves_from_base(curve_base)
  out_curve <- file.path("results", paste0(overlay_prefix, "_curve_overlay.csv"))
  out_fig <- file.path("figures", paste0(overlay_prefix, "_figure.png"))

  write.csv(curve_overlay, out_curve, row.names = FALSE)
  plot_overlay(curve_overlay, sel_all, out_fig)

  message("[done] redraw outputs:")
  message(" - ", out_curve)
  message(" - ", out_fig)
}

main()
