#!/usr/bin/env Rscript

# Primary three-sample MR: circulating OPG -> type 2 diabetes.
# Run after 01_data_processing.R.

# -----------------------------------------------------------------------------
# 1. Settings and data
# -----------------------------------------------------------------------------

data_dir <- "data"
result_dir <- file.path("results", "primary_mr")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

exposure_file <- file.path(data_dir, "primary_OPG_exposure.csv")
outcome_file <- file.path(data_dir, "primary_T2D_outcome.csv")
mr_methods <- c("mr_ivw_fe", "mr_weighted_median", "mr_egger_regression")

exposure <- TwoSampleMR::read_exposure_data(
  exposure_file,
  sep = ",",
  snp_col = "rsid", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ", pval_col = "LOG10P", log_pval = TRUE,
  samplesize_col = "N"
)
exposure$exposure <- "Circulating OPG"

outcome <- TwoSampleMR::read_outcome_data(
  outcome_file,
  snps = exposure$SNP,
  sep = ",",
  snp_col = "SNP", beta_col = "Beta", se_col = "SE",
  eaf_col = "EAF", effect_allele_col = "EffectAllele",
  other_allele_col = "NonEffectAllele", pval_col = "Pval",
  ncase_col = "Ncases", ncontrol_col = "Ncontrols",
  samplesize_col = "Neff", chr_col = "Chromsome", pos_col = "Position"
)
outcome$outcome <- "Type 2 diabetes"

# -----------------------------------------------------------------------------
# 2. Harmonization and instrument strength
# -----------------------------------------------------------------------------

dat <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
dat <- dat[dat$mr_keep, , drop = FALSE]
if (nrow(dat) == 0) stop("No instruments remained after harmonization.")

write.csv(dat, file.path(result_dir, "primary_harmonised_data.csv"), row.names = FALSE)

f_stats <- data.frame(
  SNP = dat$SNP,
  effect_allele = dat$effect_allele.exposure,
  other_allele = dat$other_allele.exposure,
  eaf = dat$eaf.exposure,
  beta = dat$beta.exposure,
  se = dat$se.exposure,
  p_value = dat$pval.exposure,
  F_statistic = (dat$beta.exposure / dat$se.exposure)^2
)
if ("samplesize.exposure" %in% names(dat)) f_stats$sample_size <- dat$samplesize.exposure
write.csv(f_stats, file.path(result_dir, "primary_instrument_F_statistics.csv"), row.names = FALSE)

mr_input <- MendelianRandomization::mr_input(
  bx = dat$beta.exposure, bxse = dat$se.exposure,
  by = dat$beta.outcome, byse = dat$se.outcome,
  snps = dat$SNP
)

overall_f <- MendelianRandomization::mr_ivw(mr_input, model = "fixed")@Fstat
write.csv(
  data.frame(n_instruments = nrow(dat), F_statistic = as.numeric(overall_f)),
  file.path(result_dir, "primary_overall_F_statistic.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. MR and sensitivity analyses
# -----------------------------------------------------------------------------

mr_results <- TwoSampleMR::mr(dat, method_list = mr_methods)
write.csv(mr_results, file.path(result_dir, "primary_mr_results.csv"), row.names = FALSE)
write.csv(
  TwoSampleMR::generate_odds_ratios(mr_results),
  file.path(result_dir, "primary_mr_results_odds_ratios.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::mr_pleiotropy_test(dat),
  file.path(result_dir, "primary_mr_egger_intercept.csv"),
  row.names = FALSE
)
write.csv(
  TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_ivw_fe", "mr_egger_regression")),
  file.path(result_dir, "primary_heterogeneity.csv"),
  row.names = FALSE
)
write.csv(
  TwoSampleMR::mr_singlesnp(dat, all_method = "mr_ivw_fe"),
  file.path(result_dir, "primary_single_snp.csv"),
  row.names = FALSE
)
write.csv(
  TwoSampleMR::mr_leaveoneout(dat),
  file.path(result_dir, "primary_leave_one_out.csv"),
  row.names = FALSE
)

if (all(c("samplesize.exposure", "samplesize.outcome") %in% names(dat)) &&
    all(is.finite(dat$samplesize.exposure)) && all(is.finite(dat$samplesize.outcome))) {
  write.csv(
    TwoSampleMR::directionality_test(dat),
    file.path(result_dir, "primary_steiger_directionality.csv"),
    row.names = FALSE
  )
}

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
message("Primary MR complete.")
