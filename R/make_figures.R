get_table1_reference <- function() {
  expand.grid(
    alpha = c(0.01, 0.05, 0.1, 0.2),
    model_name = c("model1", "model2"),
    sigma2 = c(1, 4),
    n = c(40, 80, 160, 320, 640),
    stringsAsFactors = FALSE
  ) -> ref

  ref$paper_rate <- NA_real_

  fill_rates <- function(alpha, model_name, sigma2, vals) {
    idx <- ref$alpha == alpha & ref$model_name == model_name & ref$sigma2 == sigma2
    ref$paper_rate[idx] <<- vals
  }

  fill_rates(0.01, "model1", 1, c(0, 0.76, 1, 1, 1))
  fill_rates(0.01, "model1", 4, c(0, 0.27, 0.98, 1, 1))
  fill_rates(0.01, "model2", 1, c(0, 0.29, 0.98, 1, 1))
  fill_rates(0.01, "model2", 4, c(0, 0, 0.24, 0.96, 1))

  fill_rates(0.05, "model1", 1, c(0.11, 0.97, 1, 1, 1))
  fill_rates(0.05, "model1", 4, c(0.03, 0.69, 1, 1, 1))
  fill_rates(0.05, "model2", 1, c(0, 0.76, 1, 1, 1))
  fill_rates(0.05, "model2", 4, c(0, 0.07, 0.66, 1, 1))

  fill_rates(0.1, "model1", 1, c(0.32, 1, 1, 1, 1))
  fill_rates(0.1, "model1", 4, c(0.22, 0.86, 1, 1, 1))
  fill_rates(0.1, "model2", 1, c(0, 0.92, 1, 1, 1))
  fill_rates(0.1, "model2", 4, c(0, 0.19, 0.8, 1, 1))

  fill_rates(0.2, "model1", 1, c(0.9, 1, 1, 1, 1))
  fill_rates(0.2, "model1", 4, c(0.6, 0.99, 0.99, 1, 1))
  fill_rates(0.2, "model2", 1, c(0.34, 0.92, 1, 1, 0.98))
  fill_rates(0.2, "model2", 4, c(0.09, 0.4, 0.85, 1, 1))

  ref
}

make_table1_compare <- function(cvc_summary, alpha = 0.05, cvc_method = "cvc_pmax") {
  ref <- get_table1_reference()
  cur <- cvc_summary[cvc_summary$method == cvc_method & cvc_summary$alpha == alpha,
                     c("model_name", "sigma2", "n", "alpha", "selection_rate")]
  out <- merge(cur, ref, by = c("alpha", "model_name", "sigma2", "n"), all.x = TRUE)
  out$diff <- out$selection_rate - out$paper_rate
  out[order(out$model_name, out$sigma2, out$n), ]
}

plot_sim1_selection <- function(cvc_summary, other_summary, out_path, alpha = 0.05) {
  safe_dir_create(dirname(out_path))
  png(filename = out_path, width = 1300, height = 900, res = 120)

  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

  for (model_name in c("model1", "model2")) {
    for (sigma2 in c(1, 4)) {
      sub_cvc <- cvc_summary[cvc_summary$model_name == model_name & cvc_summary$sigma2 == sigma2 & cvc_summary$alpha == alpha, ]
      sub_other <- other_summary[other_summary$model_name == model_name & other_summary$sigma2 == sigma2, ]
      x <- sort(unique(sub_other$n))

      plot(NA, xlim = range(x), ylim = c(0, 1), xlab = "n", ylab = "Correct Selection Rate",
           main = paste0(model_name, ", sigma2=", sigma2, ", α=", alpha))
      grid()

      cvc_pmax <- sub_cvc[sub_cvc$method == "cvc_pmax", ]
      cvc_pmax <- cvc_pmax[order(cvc_pmax$n), ]
      lines(cvc_pmax$n, cvc_pmax$selection_rate, type = "b", pch = 16, col = "#1F77B4", lwd = 2)

      cvc_frac <- sub_cvc[sub_cvc$method == "cvc_frac", ]
      cvc_frac <- cvc_frac[order(cvc_frac$n), ]
      lines(cvc_frac$n, cvc_frac$selection_rate, type = "b", pch = 15, col = "#FF7F0E", lwd = 2)

      cv <- sub_other[sub_other$method == "cv", ]
      cv <- cv[order(cv$n), ]
      lines(cv$n, cv$selection_rate, type = "b", pch = 17, col = "#2CA02C", lwd = 2, lty = 2)

      bic <- sub_other[sub_other$method == "bic", ]
      bic <- bic[order(bic$n), ]
      lines(bic$n, bic$selection_rate, type = "b", pch = 18, col = "#D62728", lwd = 2, lty = 3)

      legend("bottomright",
             legend = c("cvc_pmax", "cvc_frac", "cv", "bic"),
             col = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"),
             lty = c(1, 1, 2, 3),
             pch = c(16, 15, 17, 18),
             cex = 0.85,
             bty = "n")
    }
  }

  dev.off()
}

plot_sim1_nonrejected <- function(raw_df, out_path, alpha = 0.05) {
  safe_dir_create(dirname(out_path))
  png(filename = out_path, width = 1300, height = 900, res = 120)
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

  sub <- raw_df[
    (raw_df$method %in% c("cvc_pmax", "fy")) &
      (is.na(raw_df$alpha) | raw_df$alpha == alpha),
    c("model_name", "sigma2", "n", "method", "acv_size")
  ]
  sm_mean <- aggregate(acv_size ~ model_name + sigma2 + n + method, data = sub, FUN = mean)
  sm_sd <- aggregate(acv_size ~ model_name + sigma2 + n + method, data = sub, FUN = sd)
  names(sm_mean)[names(sm_mean) == "acv_size"] <- "mean_acv_size"
  names(sm_sd)[names(sm_sd) == "acv_size"] <- "sd_acv_size"
  sm <- merge(sm_mean, sm_sd, by = c("model_name", "sigma2", "n", "method"))
  sm$sd_acv_size[is.na(sm$sd_acv_size)] <- 0

  for (model_name in c("model1", "model2")) {
    for (sigma2 in c(1, 4)) {
      cur <- sm[sm$model_name == model_name & sm$sigma2 == sigma2, ]
      x <- sort(unique(cur$n))
      y_max <- max(cur$mean_acv_size + cur$sd_acv_size, na.rm = TRUE)

      plot(NA, xlim = range(x), ylim = c(0, max(1, y_max * 1.1)),
           xlab = "n", ylab = "Mean # Non-Rejected Models",
           main = paste0(model_name, ", sigma2=", sigma2, ", α=", alpha))
      grid()

      cvc <- cur[cur$method == "cvc_pmax", ]
      cvc <- cvc[order(cvc$n), ]
      cvc_upper <- cvc$mean_acv_size + cvc$sd_acv_size
      cvc_lower <- pmax(0, cvc$mean_acv_size - cvc$sd_acv_size)
      polygon(c(cvc$n, rev(cvc$n)), c(cvc_upper, rev(cvc_lower)),
              col = adjustcolor("#1F77B4", alpha.f = 0.18), border = NA)
      lines(cvc$n, cvc$mean_acv_size, type = "b", pch = 16, col = "#1F77B4", lwd = 2)

      fy <- cur[cur$method == "fy", ]
      fy <- fy[order(fy$n), ]
      fy_upper <- fy$mean_acv_size + fy$sd_acv_size
      fy_lower <- pmax(0, fy$mean_acv_size - fy$sd_acv_size)
      polygon(c(fy$n, rev(fy$n)), c(fy_upper, rev(fy_lower)),
              col = adjustcolor("#D62728", alpha.f = 0.16), border = NA)
      lines(fy$n, fy$mean_acv_size, type = "b", pch = 17, col = "#D62728", lwd = 2, lty = 2)

      legend("topright",
             legend = c("cvc (mean ± 1 SD)", "fy (mean ± 1 SD)"),
             col = c("#1F77B4", "#D62728"),
             lty = c(1, 2),
             pch = c(16, 17),
             cex = 0.9,
             bty = "n")
    }
  }

  dev.off()
}
