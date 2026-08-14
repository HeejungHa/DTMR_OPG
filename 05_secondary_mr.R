#!/usr/bin/env Rscript

# Secondary MR: circulating OPG -> fasting glucose, fasting insulin, and HbA1c.

# -----------------------------------------------------------------------------
# 1. Settings and exposure data
# -----------------------------------------------------------------------------

data_dir <- "data"
result_dir <- file.path("results", "secondary_mr")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

exposure <- TwoSampleMR::read_exposure_data(
  file.path(data_dir, "primary_OPG_exposure.csv"),
  sep = ",",
  snp_col = "rsid", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ", pval_col = "LOG10P", log_pval = TRUE,
  samplesize_col = "N"
)
exposure$exposure <- "Circulating OPG"

traits <- data.frame(
  key = c("fasting_glucose", "fasting_insulin", "hba1c"),
  name = c("Fasting glucose", "Fasting insulin", "HbA1c"),
  file = c("MAGIC1000G_FG_EUR.tsv", "MAGIC1000G_FI_EUR.tsv", "MAGIC1000G_HbA1c_EUR.tsv"),
  stringsAsFactors = FALSE
)

mr_methods <- c("mr_ivw_fe", "mr_weighted_median", "mr_egger_regression")
z975 <- stats::qnorm(0.975)

# -----------------------------------------------------------------------------
# 2. Analysis function
# -----------------------------------------------------------------------------

run_secondary_mr <- function(key, trait_name, filename) {
  out_dir <- file.path(result_dir, key)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  outcome <- TwoSampleMR::read_outcome_data(
    file.path(data_dir, filename),
    snps = exposure$SNP,
    sep = "\t",
    snp_col = "variant", beta_col = "beta", se_col = "standard_error",
    eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
    other_allele_col = "other_allele", pval_col = "p_value",
    samplesize_col = "sample_size", chr_col = "chromosome",
    pos_col = "base_pair_location"
  )
  outcome$outcome <- trait_name

  dat <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
  dat <- dat[dat$mr_keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No instruments remained for ", trait_name, ".")

  mr_results <- TwoSampleMR::mr(dat, method_list = mr_methods)
  mr_results$ci_lower <- mr_results$b - z975 * mr_results$se
  mr_results$ci_upper <- mr_results$b + z975 * mr_results$se

  heterogeneity <- TwoSampleMR::mr_heterogeneity(dat, method_list = "mr_ivw_fe")
  if (nrow(heterogeneity) > 0) {
    heterogeneity$I2_percent <- with(
      heterogeneity,
      ifelse(is.finite(Q) & Q > 0, pmax(0, (Q - Q_df) / Q) * 100, NA_real_)
    )
  }

  write.csv(dat, file.path(out_dir, paste0(key, "_harmonised_data.csv")), row.names = FALSE)
  write.csv(mr_results, file.path(out_dir, paste0(key, "_mr_results.csv")), row.names = FALSE)
  write.csv(
    TwoSampleMR::mr_pleiotropy_test(dat),
    file.path(out_dir, paste0(key, "_mr_egger_intercept.csv")),
    row.names = FALSE
  )
  write.csv(heterogeneity, file.path(out_dir, paste0(key, "_heterogeneity.csv")), row.names = FALSE)
  write.csv(
    TwoSampleMR::mr_singlesnp(dat, all_method = "mr_ivw_fe"),
    file.path(out_dir, paste0(key, "_single_snp.csv")),
    row.names = FALSE
  )
  write.csv(
    TwoSampleMR::mr_leaveoneout(dat),
    file.path(out_dir, paste0(key, "_leave_one_out.csv")),
    row.names = FALSE
  )

  if (all(c("samplesize.exposure", "samplesize.outcome") %in% names(dat)) &&
      all(is.finite(dat$samplesize.exposure)) && all(is.finite(dat$samplesize.outcome))) {
    write.csv(
      TwoSampleMR::directionality_test(dat),
      file.path(out_dir, paste0(key, "_steiger_directionality.csv")),
      row.names = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# 3. Run analyses
# -----------------------------------------------------------------------------

for (i in seq_len(nrow(traits))) {
  run_secondary_mr(traits$key[i], traits$name[i], traits$file[i])
}

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
message("Secondary MR analyses complete.")
