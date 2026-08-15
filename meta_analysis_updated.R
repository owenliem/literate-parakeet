#!/usr/bin/env Rscript

# Reproducible aggregate-data meta-analysis for complete early ctDNA clearance
# during immune checkpoint inhibitor-based therapy in advanced NSCLC.
#
# Statistical model:
# - effect measure: log hazard ratio
# - within-study SE reconstructed from the reported 95% CI
# - random-effects model with REML heterogeneity estimator
# - modified Knapp-Hartung inference via metafor test = "adhoc"
# - Riley-type 95% prediction interval

options(stringsAsFactors = FALSE)

if (!requireNamespace("metafor", quietly = TRUE)) {
  install.packages("metafor", repos = "https://cloud.r-project.org")
}
library(metafor)

input_file <- "updated_ctdna_effects_20260815.csv"
out_dir <- "r_meta_outputs_20260815"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(input_file, check.names = FALSE)
dat$yi <- log(dat$hr)
dat$sei <- (log(dat$upper) - log(dat$lower)) / (2 * qnorm(0.975))
dat$vi <- dat$sei^2

fit_outcome <- function(subdat, outcome_name) {
  res <- rma.uni(
    yi = yi,
    vi = vi,
    data = subdat,
    method = "REML",
    test = "adhoc",
    level = 95,
    slab = study
  )
  pred <- predict(res, transf = exp, predtype = "Riley")
  w <- weights(res)

  main <- data.frame(
    outcome = outcome_name,
    k = res$k,
    n_total = sum(subdat$analysis_n),
    pooled_hr = exp(as.numeric(res$b[1])),
    ci_low = exp(res$ci.lb),
    ci_high = exp(res$ci.ub),
    tau2 = res$tau2,
    I2_percent = res$I2,
    Q = res$QE,
    Q_df = res$k - 1,
    Q_p = res$QEp,
    prediction_interval_low = pred$pi.lb,
    prediction_interval_high = pred$pi.ub,
    test = "Modified Knapp-Hartung (metafor test='adhoc')",
    tau2_estimator = "REML",
    prediction_interval = "Riley-type",
    stringsAsFactors = FALSE
  )

  study_results <- data.frame(
    outcome = outcome_name,
    study = subdat$study,
    analysis_n = subdat$analysis_n,
    hr = subdat$hr,
    ci_low = subdat$lower,
    ci_high = subdat$upper,
    log_hr = subdat$yi,
    se_log_hr = subdat$sei,
    inverse_variance_weight_percent = as.numeric(w),
    assessment = subdat$assessment,
    source = subdat$source,
    stringsAsFactors = FALSE
  )

  loo_rows <- vector("list", nrow(subdat))
  for (i in seq_len(nrow(subdat))) {
    d2 <- subdat[-i, , drop = FALSE]
    r2 <- rma.uni(yi = yi, vi = vi, data = d2, method = "REML", test = "adhoc", level = 95)
    p2 <- predict(r2, transf = exp, predtype = "Riley")
    loo_rows[[i]] <- data.frame(
      outcome = outcome_name,
      omitted_study = subdat$study[i],
      k = r2$k,
      n_remaining = sum(d2$analysis_n),
      pooled_hr = exp(as.numeric(r2$b[1])),
      ci_low = exp(r2$ci.lb),
      ci_high = exp(r2$ci.ub),
      I2_percent = r2$I2,
      tau2 = r2$tau2,
      prediction_interval_low = p2$pi.lb,
      prediction_interval_high = p2$pi.ub,
      stringsAsFactors = FALSE
    )
  }
  loo <- do.call(rbind, loo_rows)

  # Sensitivity analysis excluding the figure-derived Fei estimate.
  d_no_fei <- subdat[!grepl("^Fei", subdat$study), , drop = FALSE]
  r_no_fei <- rma.uni(yi = yi, vi = vi, data = d_no_fei, method = "REML", test = "adhoc", level = 95)
  p_no_fei <- predict(r_no_fei, transf = exp, predtype = "Riley")
  sensitivity <- data.frame(
    outcome = outcome_name,
    analysis = "Exclude Fei et al. 2025 figure-derived estimate",
    k = r_no_fei$k,
    n_total = sum(d_no_fei$analysis_n),
    pooled_hr = exp(as.numeric(r_no_fei$b[1])),
    ci_low = exp(r_no_fei$ci.lb),
    ci_high = exp(r_no_fei$ci.ub),
    I2_percent = r_no_fei$I2,
    tau2 = r_no_fei$tau2,
    prediction_interval_low = p_no_fei$pi.lb,
    prediction_interval_high = p_no_fei$pi.ub,
    stringsAsFactors = FALSE
  )

  # Publication-quality forest plot.
  make_forest <- function(filename, type = c("png", "tiff")) {
    type <- match.arg(type)
    if (type == "png") {
      png(filename, width = 2600, height = 1800, res = 300)
    } else {
      tiff(filename, width = 2600, height = 1800, res = 300, compression = "lzw")
    }
    oldpar <- par(mar = c(5, 5, 4, 2), family = "sans")
    on.exit({par(oldpar); dev.off()}, add = TRUE)
    forest(
      res,
      atransf = exp,
      at = log(c(0.05, 0.10, 0.20, 0.50, 1.00, 2.00)),
      xlim = c(-6.8, 3.2),
      alim = log(c(0.05, 2.00)),
      ilab = cbind(subdat$analysis_n, sprintf("%.1f", w)),
      ilab.xpos = c(-4.35, -3.10),
      cex = 0.95,
      psize = 1.15,
      refline = 0,
      xlab = "Hazard ratio (complete early ctDNA clearance vs persistent ctDNA)",
      mlab = sprintf("Random-effects model (REML; modified Hartung-Knapp), I^2 = %.0f%%", res$I2),
      header = c("Study", "HR [95% CI]")
    )
    text(-4.35, res$k + 2, "N", font = 2, cex = 0.95)
    text(-3.10, res$k + 2, "Weight, %", font = 2, cex = 0.95)
    title(main = ifelse(outcome_name == "PFS",
                        "Complete early ctDNA clearance and progression-free survival",
                        "Complete early ctDNA clearance and overall survival"),
          cex.main = 1.15)
  }

  prefix <- tolower(outcome_name)
  make_forest(file.path(out_dir, paste0("forest_", prefix, ".png")), "png")
  make_forest(file.path(out_dir, paste0("forest_", prefix, ".tiff")), "tiff")

  list(model = res, main = main, studies = study_results, loo = loo, sensitivity = sensitivity)
}

results <- lapply(c("PFS", "OS"), function(outcome_name) {
  fit_outcome(dat[dat$outcome == outcome_name, , drop = FALSE], outcome_name)
})
names(results) <- c("PFS", "OS")

write.csv(do.call(rbind, lapply(results, `[[`, "main")),
          file.path(out_dir, "meta_analysis_results_R.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "studies")),
          file.path(out_dir, "study_effects_and_weights_R.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "loo")),
          file.path(out_dir, "leave_one_out_results_R.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "sensitivity")),
          file.path(out_dir, "sensitivity_results_R.csv"), row.names = FALSE)
write.csv(dat, file.path(out_dir, "analysis_input_with_log_effects.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "R_session_info.txt"))

summary_lines <- c(
  sprintf("Analysis run: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Model: inverse-variance random effects; REML tau-squared; modified Knapp-Hartung; Riley-type prediction interval.",
  "",
  apply(do.call(rbind, lapply(results, `[[`, "main")), 1, function(x) {
    sprintf("%s: %s studies, %s patients; HR %.3f (95%% CI %.3f-%.3f); I2 %.1f%%; 95%% PI %.3f-%.3f.",
            x[["outcome"]], x[["k"]], x[["n_total"]], as.numeric(x[["pooled_hr"]]),
            as.numeric(x[["ci_low"]]), as.numeric(x[["ci_high"]]), as.numeric(x[["I2_percent"]]),
            as.numeric(x[["prediction_interval_low"]]), as.numeric(x[["prediction_interval_high"]]))
  })
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
