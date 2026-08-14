#!/usr/bin/env Rscript

# Alternative two-sample MR: circulating OPG -> type 2 diabetes.
#
# Workflow:
#   1) Select genome-wide significant OPG cis-pQTLs from UKB-PPP
#   2) LD clumping (r2 = 0.1, EUR)
#   3) Harmonize with T2D and perform radial MR
#   4) Exclude rs1485286 identified in the original radial MR analysis
#   5) CREATE alternative_OPG_exposure_rs1485286_excluded.csv
#   6) Re-read that file and run the final alternative two-sample MR
#
# RadialMR installation, if needed:
# remotes::install_github("WSpiller/RadialMR")

required_packages <- c(
  "TwoSampleMR",
  "MendelianRandomization",
  "RadialMR",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------

data_dir <- "data"
result_dir <- file.path("results", "alternative_two_sample_mr")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Raw UKB-PPP OPG pQTL data used for the alternative analysis.
opg_file <- file.path(
  data_dir,
  "TNFRSF11B_O00300_OID20735_v1_Inflammation_chr8.csv"
)

# Prepared T2D outcome data.
t2d_file <- file.path(data_dir, "alternative_T2D_outcome.csv")

# This file is CREATED in Section 5 below and then used in Section 6.
excluded_exposure_file <- file.path(
  result_dir,
  "alternative_OPG_exposure_rs1485286_excluded.csv"
)

if (!file.exists(opg_file)) stop("OPG pQTL file not found: ", opg_file)
if (!file.exists(t2d_file)) stop("T2D outcome file not found: ", t2d_file)

# Coordinates used in the original alternative two-sample MR workflow.
gene_start <- 118923557
gene_end <- 118951885
cis_window_bp <- 100000

p_threshold <- 5e-8
clump_r2 <- 0.1
outlier_snp <- "rs1485286"

mr_methods <- c(
  "mr_ivw_mre",
  "mr_weighted_median",
  "mr_egger_regression"
)

# -----------------------------------------------------------------------------
# 2. Construct alternative OPG instruments
# -----------------------------------------------------------------------------

opg_raw <- read.csv(
  opg_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "rsid", "BETA", "SE", "ALLELE1", "ALLELE0",
  "A1FREQ", "LOG10P", "N", "GENPOS"
)
missing_cols <- setdiff(required_cols, names(opg_raw))
if (length(missing_cols) > 0) {
  stop("Missing OPG column(s): ", paste(missing_cols, collapse = ", "))
}

# LOG10P is -log10(P).
opg_raw$P <- 10^(-opg_raw$LOG10P)

# Genome-wide significant variants within the OPG cis +/-100 kb region.
opg_cis <- opg_raw[
  opg_raw$P < p_threshold &
    opg_raw$GENPOS >= gene_start - cis_window_bp &
    opg_raw$GENPOS <= gene_end + cis_window_bp,
  ,
  drop = FALSE
]
opg_cis <- opg_cis[!duplicated(opg_cis$rsid), , drop = FALSE]

if (nrow(opg_cis) == 0) {
  stop("No OPG variants remained after P-value and cis-region filtering.")
}

exposure <- TwoSampleMR::format_data(
  opg_cis,
  type = "exposure",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P",
  samplesize_col = "N"
)
exposure$exposure <- "Circulating OPG - alternative analysis"

exposure_clumped <- TwoSampleMR::clump_data(
  exposure,
  clump_r2 = clump_r2,
  pop = "EUR"
)

if (nrow(exposure_clumped) == 0) {
  stop("No OPG instruments remained after LD clumping.")
}

write.csv(
  exposure_clumped,
  file.path(result_dir, "alternative_OPG_exposure_before_outlier_exclusion.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Read T2D outcome and harmonize before radial MR
# -----------------------------------------------------------------------------

outcome <- TwoSampleMR::read_outcome_data(
  filename = t2d_file,
  snps = exposure_clumped$SNP,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EffectAllele",
  other_allele_col = "NonEffectAllele",
  pval_col = "Pval",
  ncase_col = "Ncases",
  ncontrol_col = "Ncontrols",
  samplesize_col = "Neff",
  chr_col = "Chromsome",
  pos_col = "Position",
  log_pval = FALSE
)
outcome$outcome <- "Type 2 diabetes"

dat <- TwoSampleMR::harmonise_data(
  exposure_clumped,
  outcome,
  action = 2
)
dat <- dat[dat$mr_keep, , drop = FALSE]

if (nrow(dat) == 0) {
  stop("No instruments remained after harmonization.")
}

write.csv(
  dat,
  file.path(result_dir, "alternative_harmonised_before_outlier_exclusion.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::mr_heterogeneity(dat),
  file.path(result_dir, "alternative_heterogeneity_before_outlier_exclusion.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Radial MR
# -----------------------------------------------------------------------------

# Same radial-analysis settings used in the original code.
dat_clean <- dat

radial_in <- RadialMR::format_radial(
  BXG = dat_clean$beta.exposure,
  BYG = dat_clean$beta.outcome,
  seBXG = dat_clean$se.exposure,
  seBYG = dat_clean$se.outcome,
  RSID = dat_clean$SNP
)

ivw <- RadialMR::ivw_radial(
  radial_in,
  alpha = 0.05,
  weights = 3,
  tol = 1e-4,
  summary = TRUE
)

egg <- RadialMR::egger_radial(
  radial_in,
  alpha = 0.05,
  weights = 3,
  summary = TRUE
)

saveRDS(ivw, file.path(result_dir, "radial_ivw_results.rds"))
saveRDS(egg, file.path(result_dir, "radial_egger_results.rds"))

if (is.data.frame(ivw$outliers)) {
  write.csv(
    ivw$outliers,
    file.path(result_dir, "radial_ivw_outliers.csv"),
    row.names = FALSE
  )
}

if (is.data.frame(egg$outliers)) {
  write.csv(
    egg$outliers,
    file.path(result_dir, "radial_egger_outliers.csv"),
    row.names = FALSE
  )
}

p_ivw <- RadialMR::plot_radial(
  ivw,
  radial_scale = TRUE,
  show_outliers = FALSE,
  scale_match = TRUE
)

p_egg <- RadialMR::plot_radial(
  egg,
  radial_scale = TRUE,
  show_outliers = FALSE,
  scale_match = TRUE
)

p_both <- RadialMR::plot_radial(
  c(ivw, egg),
  radial_scale = TRUE,
  show_outliers = FALSE,
  scale_match = TRUE
)

ggplot2::ggsave(
  file.path(result_dir, "mr_radial_ivw.png"),
  p_ivw,
  width = 6,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  file.path(result_dir, "mr_radial_egg.png"),
  p_egg,
  width = 6,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  file.path(result_dir, "mr_radial_both.png"),
  p_both,
  width = 6,
  height = 5,
  dpi = 600
)

# -----------------------------------------------------------------------------
# 5. Exclude rs1485286 and CREATE the final exposure CSV
# -----------------------------------------------------------------------------

# rs1485286 was the radial-MR outlier removed in the original analysis.
if (!outlier_snp %in% exposure_clumped$SNP) {
  stop(outlier_snp, " was not present in the clumped OPG instrument set.")
}

exposure_no_outlier <- exposure_clumped[
  exposure_clumped$SNP != outlier_snp,
  ,
  drop = FALSE
]

if (nrow(exposure_no_outlier) == 0) {
  stop("No OPG instruments remained after excluding ", outlier_snp, ".")
}

# Export a clean input-format CSV so that the exact exposure dataset used in
# the final alternative MR analysis is preserved and can be inspected directly.
alternative_exposure_export <- data.frame(
  SNP = exposure_no_outlier$SNP,
  Beta = exposure_no_outlier$beta.exposure,
  SE = exposure_no_outlier$se.exposure,
  EffectAllele = exposure_no_outlier$effect_allele.exposure,
  NonEffectAllele = exposure_no_outlier$other_allele.exposure,
  EAF = exposure_no_outlier$eaf.exposure,
  Pval = exposure_no_outlier$pval.exposure,
  N = exposure_no_outlier$samplesize.exposure,
  stringsAsFactors = FALSE
)

# >>> THIS LINE CREATES alternative_OPG_exposure_rs1485286_excluded.csv <<<
write.csv(
  alternative_exposure_export,
  excluded_exposure_file,
  row.names = FALSE
)

message("Created: ", excluded_exposure_file)

# -----------------------------------------------------------------------------
# 6. Re-read the created exposure CSV and harmonize for the final MR
# -----------------------------------------------------------------------------

# Re-reading the file makes the analysis path explicit and verifies that the
# saved CSV itself is the exposure dataset used for the final analysis.
exposure_final <- TwoSampleMR::read_exposure_data(
  filename = excluded_exposure_file,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "EffectAllele",
  other_allele_col = "NonEffectAllele",
  eaf_col = "EAF",
  pval_col = "Pval",
  samplesize_col = "N",
  log_pval = FALSE
)
exposure_final$exposure <- "Circulating OPG - alternative analysis"

outcome_final <- TwoSampleMR::read_outcome_data(
  filename = t2d_file,
  snps = exposure_final$SNP,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EffectAllele",
  other_allele_col = "NonEffectAllele",
  pval_col = "Pval",
  ncase_col = "Ncases",
  ncontrol_col = "Ncontrols",
  samplesize_col = "Neff",
  chr_col = "Chromsome",
  pos_col = "Position",
  log_pval = FALSE
)
outcome_final$outcome <- "Type 2 diabetes"

dat_final <- TwoSampleMR::harmonise_data(
  exposure_final,
  outcome_final,
  action = 2
)
dat_final <- dat_final[dat_final$mr_keep, , drop = FALSE]

if (nrow(dat_final) == 0) {
  stop("No instruments remained after final harmonization.")
}

write.csv(
  dat_final,
  file.path(result_dir, "alternative_harmonised_data.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Instrument strength
# -----------------------------------------------------------------------------

f_stats <- data.frame(
  SNP = dat_final$SNP,
  effect_allele = dat_final$effect_allele.exposure,
  other_allele = dat_final$other_allele.exposure,
  eaf = dat_final$eaf.exposure,
  beta = dat_final$beta.exposure,
  se = dat_final$se.exposure,
  p_value = dat_final$pval.exposure,
  F_statistic = (dat_final$beta.exposure / dat_final$se.exposure)^2
)

write.csv(
  f_stats,
  file.path(result_dir, "alternative_instrument_F_statistics.csv"),
  row.names = FALSE
)

mr_input <- MendelianRandomization::mr_input(
  bx = dat_final$beta.exposure,
  bxse = dat_final$se.exposure,
  by = dat_final$beta.outcome,
  byse = dat_final$se.outcome,
  snps = dat_final$SNP
)

overall_f <- MendelianRandomization::mr_ivw(
  mr_input,
  model = "random"
)@Fstat

write.csv(
  data.frame(
    n_instruments = nrow(dat_final),
    F_statistic = as.numeric(overall_f)
  ),
  file.path(result_dir, "alternative_overall_F_statistic.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Alternative two-sample MR
# -----------------------------------------------------------------------------

mr_results <- TwoSampleMR::mr(
  dat_final,
  method_list = mr_methods
)

write.csv(
  mr_results,
  file.path(result_dir, "alternative_mr_results.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::generate_odds_ratios(mr_results),
  file.path(result_dir, "alternative_mr_results_odds_ratios.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. Sensitivity analyses
# -----------------------------------------------------------------------------

write.csv(
  TwoSampleMR::mr_pleiotropy_test(dat_final),
  file.path(result_dir, "alternative_mr_egger_intercept.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::mr_heterogeneity(
    dat_final,
    method_list = c("mr_ivw_mre", "mr_egger_regression")
  ),
  file.path(result_dir, "alternative_heterogeneity.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::mr_singlesnp(
    dat_final,
    all_method = "mr_ivw_mre"
  ),
  file.path(result_dir, "alternative_single_snp.csv"),
  row.names = FALSE
)

write.csv(
  TwoSampleMR::mr_leaveoneout(dat_final),
  file.path(result_dir, "alternative_leave_one_out.csv"),
  row.names = FALSE
)

if (
  all(c("samplesize.exposure", "samplesize.outcome") %in% names(dat_final)) &&
    all(is.finite(dat_final$samplesize.exposure)) &&
    all(is.finite(dat_final$samplesize.outcome))
) {
  write.csv(
    TwoSampleMR::directionality_test(dat_final),
    file.path(result_dir, "alternative_steiger_directionality.csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 10. Reproducibility information
# -----------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "sessionInfo.txt")
)

message("Alternative two-sample MR complete.")
