source("R/utils.R")

plot_with_cv <- function(curve_df, sel_all, out_path) {
  png(out_path, width = 1700, height = 760, res = 120)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 4))

  if (!("fdp_mean" %in% names(curve_df))) {
    curve_df$fdp_mean <- curve_df$true_fdr_mean
  }

  mse_cols <- c(cv = "#1F77B4", cvc = "#FF7F0E", ncv0 = "#7F7F7F", ncv1 = "#2CA02C")
  fdr_cols <- mse_cols
  methods <- c("cv", "cvc", "ncv0", "ncv1")
  draw_order <- c("cvc", "ncv0", "ncv1", "cv")

  for (scenario in c("independent", "correlated")) {
    cur <- curve_df[curve_df$scenario == scenario, ]
    cur <- cur[order(cur$lambda, decreasing = TRUE), ]
    x <- sort(unique(cur$lambda), decreasing = TRUE)

    y_mse_max <- max(cur$mse_mean, na.rm = TRUE)
    plot(x, rep(NA_real_, length(x)),
         type = "n", xlab = "lambda", ylab = "MSE", log = "x",
         ylim = c(0, y_mse_max * 1.05),
         main = ifelse(scenario == "independent", "Independent", "Correlated (AR(0.8))"))

    for (m in draw_order) {
      tmp <- cur[cur$method == m, ]
      tmp <- tmp[order(tmp$lambda, decreasing = TRUE), ]
      lwd <- ifelse(m == "cv", 3.2, 2.2)
      lines(tmp$lambda, tmp$mse_mean, lwd = lwd, col = mse_cols[m], lty = 1)
      if (m == "cv") {
        idx_mark <- unique(round(seq(1, nrow(tmp), length.out = min(20, nrow(tmp)))))
        points(tmp$lambda[idx_mark], tmp$mse_mean[idx_mark], pch = 16, cex = 0.45, col = mse_cols[m])
      }
    }

    par(new = TRUE)
    plot(x, rep(NA_real_, length(x)), type = "n", axes = FALSE, xlab = "", ylab = "", log = "x", ylim = c(0, 1))
    for (m in draw_order) {
      tmp <- cur[cur$method == m, ]
      tmp <- tmp[order(tmp$lambda, decreasing = TRUE), ]
      lwd <- ifelse(m == "cv", 2.2, 1.6)
      lines(tmp$lambda, tmp$true_fdr_mean, lwd = lwd, col = fdr_cols[m], lty = 2)
      lines(tmp$lambda, tmp$est_fdr_mean, lwd = lwd, col = fdr_cols[m], lty = 3)
      lines(tmp$lambda, tmp$fdp_mean, lwd = lwd, col = fdr_cols[m], lty = 4)
      if (m == "cv") {
        idx_mark <- unique(round(seq(1, nrow(tmp), length.out = min(20, nrow(tmp)))))
        points(tmp$lambda[idx_mark], tmp$fdp_mean[idx_mark], pch = 16, cex = 0.35, col = fdr_cols[m])
      }
    }
    axis(4)
    mtext("FDR / FDP", side = 4, line = 2)

    sub <- sel_all[sel_all$scenario == scenario, ]
    mean_by_method <- aggregate(lambda ~ method, data = sub, FUN = mean)
    xline <- function(method, col, lty, lwd = 2) {
      val <- mean_by_method$lambda[mean_by_method$method == method]
      if (length(val) == 1L) abline(v = val, col = col, lty = lty, lwd = lwd)
    }
    xline("cv_1se", "#1F77B4", 2, 2.2)
    xline("ncv0_1se", "#7F7F7F", 2, 2.2)
    xline("ncv1_1se", "#2CA02C", 2, 2.2)

    legend("topleft",
           legend = c("cv (highlight)", "cvc", "ncv0", "ncv1",
                      "MSE", "FDR_true", "FDR_hat", "FDP",
                      "cv_1se", "ncv0_1se", "ncv1_1se"),
           col = c("#1F77B4", "#FF7F0E", "#7F7F7F", "#2CA02C",
                   "#111111", "#111111", "#111111", "#111111",
                   "#1F77B4", "#7F7F7F", "#2CA02C"),
           lty = c(1, 1, 1, 1, 1, 2, 3, 4, 2, 2, 2),
           lwd = c(3.2, 2.2, 2.2, 2.2, 2.2, 1.8, 1.8, 1.8, 2.2, 2.2, 2.2),
           cex = 0.78, bg = "white")
  }

  dev.off()
}

main <- function() {
  curve_path <- file.path("results", "etfdr_fig1_smoke_ncv1_overlay_curve_overlay.csv")
  sel_path <- file.path("results", "etfdr_fig1_smoke_ncv1_overlay_all_selected.csv")
  out_path <- file.path("figures", "etfdr_fig1_smoke_ncv1_overlay_with_cv_figure.png")

  if (!file.exists(curve_path) || !file.exists(sel_path)) {
    stop("Missing smoke results. Please run smoke first.")
  }

  curve_df <- read.csv(curve_path)
  sel_all <- read.csv(sel_path)
  plot_with_cv(curve_df, sel_all, out_path)

  message("[done] ", out_path)
}

main()
