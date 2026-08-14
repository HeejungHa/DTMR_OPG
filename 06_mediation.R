#!/usr/bin/env Rscript

# Two-step MR mediation: OPG -> HOMA-B -> type 2 diabetes.
# Run after 02_primary_mr.R.

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------

data_dir <- "data"
result_dir <- file.path("results", "mediation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

opg_file <- file.path(data_dir, "primary_OPG_exposure.csv")
homa_b_file <- file.path(data_dir, "mediation_HOMA_B_GWAS.tsv")
t2d_file <- file.path(data_dir, "mediation_T2D_GWAS_EUR.csv")
primary_mr_file <- file.path("results", "primary_mr", "primary_mr_results.csv")

homa_b_gwas_id <- "ieu-b-117"
homa_b_sample_size <- 36466
mr_methods <- c("mr_ivw_fe", "mr_weighted_median", "mr_egger_regression")
z975 <- stats::qnorm(0.975)

extract_ivw_fe <- function(x, label) {
  keep <- grepl("inverse variance weighted.*fixed", x$method, ignore.case = TRUE)
  if (sum(keep) != 1) stop("Could not identify the fixed-effect IVW estimate for ", label, ".")
  x[keep, , drop = FALSE]
}

save_sensitivity <- function(dat, prefix) {
  write.csv(
    TwoSampleMR::mr_pleiotropy_test(dat),
    file.path(result_dir, paste0(prefix, "_mr_egger_intercept.csv")),
    row.names = FALSE
  )
  write.csv(
    TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_ivw_fe", "mr_egger_regression")),
    file.path(result_dir, paste0(prefix, "_heterogeneity.csv")),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. Path a: OPG -> HOMA-B
# -----------------------------------------------------------------------------

opg <- TwoSampleMR::read_exposure_data(
  opg_file,
  sep = ",",
  snp_col = "rsid", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ", pval_col = "LOG10P", log_pval = TRUE,
  samplesize_col = "N"
)
opg$exposure <- "Circulating OPG"

homa_b_outcome <- TwoSampleMR::read_outcome_data(
  homa_b_file,
  snps = opg$SNP,
  sep = "\t",
  snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
  eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
  other_allele_col = "other_allele", pval_col = "p_value"
)
homa_b_outcome$outcome <- "HOMA-B"
homa_b_outcome$samplesize.outcome <- homa_b_sample_size

path_a_dat <- TwoSampleMR::harmonise_data(opg, homa_b_outcome, action = 2)
path_a_dat <- path_a_dat[path_a_dat$mr_keep, , drop = FALSE]
if (nrow(path_a_dat) == 0) stop("No instruments remained for OPG -> HOMA-B.")

path_a_results <- TwoSampleMR::mr(path_a_dat, method_list = mr_methods)
path_a_ivw <- extract_ivw_fe(path_a_results, "OPG -> HOMA-B")

write.csv(path_a_dat, file.path(result_dir, "path_a_OPG_to_HOMA_B_harmonised.csv"), row.names = FALSE)
write.csv(path_a_results, file.path(result_dir, "path_a_OPG_to_HOMA_B_mr_results.csv"), row.names = FALSE)
save_sensitivity(path_a_dat, "path_a_OPG_to_HOMA_B")

# -----------------------------------------------------------------------------
# 3. Path b: HOMA-B -> T2D
# -----------------------------------------------------------------------------

jwt <- ieugwasr::get_opengwas_jwt()
if (!nzchar(jwt)) stop("Set OPENGWAS_JWT in .Renviron before running this script.")

homa_b_exposure <- TwoSampleMR::extract_instruments(
  outcomes = homa_b_gwas_id,
  p1 = 5e-8,
  clump = TRUE,
  r2 = 0.1,
  kb = 10000,
  opengwas_jwt = jwt,
  force_server = FALSE
)
homa_b_exposure$exposure <- "HOMA-B"

t2d <- read.csv(t2d_file, stringsAsFactors = FALSE, check.names = FALSE)
required_t2d <- c(
  "Chromsome", "Position", "EffectAllele", "NonEffectAllele",
  "Beta", "SE", "EAF", "Pval", "Ncases", "Ncontrols", "Neff"
)
missing_t2d <- setdiff(required_t2d, names(t2d))
if (length(missing_t2d) > 0) stop("Missing T2D column(s): ", paste(missing_t2d, collapse = ", "))

homa_b_exposure$chr.exposure <- as.character(homa_b_exposure$chr.exposure)
t2d$Chromsome <- as.character(t2d$Chromsome)
matched <- merge(
  homa_b_exposure,
  t2d,
  by.x = c("chr.exposure", "pos.exposure"),
  by.y = c("Chromsome", "Position"),
  all = FALSE
)
if (nrow(matched) == 0) stop("No HOMA-B instruments matched T2D by chromosome/position.")

path_b_outcome <- TwoSampleMR::format_data(
  data.frame(
    SNP = matched$SNP,
    Beta = matched$Beta,
    SE = matched$SE,
    EAF = matched$EAF,
    EffectAllele = matched$EffectAllele,
    NonEffectAllele = matched$NonEffectAllele,
    Pval = matched$Pval,
    Ncases = matched$Ncases,
    Ncontrols = matched$Ncontrols,
    Neff = matched$Neff,
    Chromsome = matched$chr.exposure,
    Position = matched$pos.exposure
  ),
  type = "outcome",
  snp_col = "SNP", beta_col = "Beta", se_col = "SE", eaf_col = "EAF",
  effect_allele_col = "EffectAllele", other_allele_col = "NonEffectAllele",
  pval_col = "Pval", ncase_col = "Ncases", ncontrol_col = "Ncontrols",
  samplesize_col = "Neff", chr_col = "Chromsome", pos_col = "Position"
)
path_b_outcome$outcome <- "Type 2 diabetes"

path_b_dat <- TwoSampleMR::harmonise_data(homa_b_exposure, path_b_outcome, action = 2)
path_b_dat <- path_b_dat[path_b_dat$mr_keep, , drop = FALSE]
if (nrow(path_b_dat) == 0) stop("No instruments remained for HOMA-B -> T2D.")

path_b_results <- TwoSampleMR::mr(path_b_dat, method_list = mr_methods)
path_b_ivw <- extract_ivw_fe(path_b_results, "HOMA-B -> T2D")

write.csv(homa_b_exposure, file.path(result_dir, "path_b_HOMA_B_instruments.csv"), row.names = FALSE)
write.csv(path_b_dat, file.path(result_dir, "path_b_HOMA_B_to_T2D_harmonised.csv"), row.names = FALSE)
write.csv(path_b_results, file.path(result_dir, "path_b_HOMA_B_to_T2D_mr_results.csv"), row.names = FALSE)
write.csv(
  TwoSampleMR::generate_odds_ratios(path_b_results),
  file.path(result_dir, "path_b_HOMA_B_to_T2D_odds_ratios.csv"),
  row.names = FALSE
)
save_sensitivity(path_b_dat, "path_b_HOMA_B_to_T2D")

# -----------------------------------------------------------------------------
# 4. Indirect effect and mediation proportion
# -----------------------------------------------------------------------------

total_results <- read.csv(primary_mr_file, check.names = FALSE)
total_ivw <- extract_ivw_fe(total_results, "OPG -> T2D")

a <- path_a_ivw$b[[1]]
se_a <- path_a_ivw$se[[1]]
b <- path_b_ivw$b[[1]]
se_b <- path_b_ivw$se[[1]]
c_total <- total_ivw$b[[1]]
se_total <- total_ivw$se[[1]]

if (!all(is.finite(c(a, se_a, b, se_b, c_total, se_total)))) {
  stop("Non-finite MR estimate or standard error encountered.")
}
if (c_total == 0) stop("Total effect is zero; mediation proportion is undefined.")

indirect <- a * b
var_indirect <- b^2 * se_a^2 + a^2 * se_b^2
se_indirect <- sqrt(var_indirect)
ci_indirect <- indirect + c(-1, 1) * z975 * se_indirect
p_indirect <- 2 * stats::pnorm(-abs(indirect / se_indirect))

# Mediation proportion = indirect effect / total effect.
# A confidence interval for the mediation proportion is not calculated here.
prop <- indirect / c_total

mediation_summary <- data.frame(
  measure = c(
    "Path a: OPG -> HOMA-B",
    "Path b: HOMA-B -> T2D",
    "Total effect: OPG -> T2D",
    "Indirect effect",
    "Mediation proportion (%)"
  ),
  estimate = c(a, b, c_total, indirect, 100 * prop),
  se = c(se_a, se_b, se_total, se_indirect, NA_real_),
  ci_lower = c(a - z975 * se_a, b - z975 * se_b, c_total - z975 * se_total, ci_indirect[1], NA_real_),
  ci_upper = c(a + z975 * se_a, b + z975 * se_b, c_total + z975 * se_total, ci_indirect[2], NA_real_),
  p_value = c(path_a_ivw$pval[[1]], path_b_ivw$pval[[1]], total_ivw$pval[[1]], p_indirect, NA_real_)
)

write.csv(mediation_summary, file.path(result_dir, "mediation_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
message("Mediation analysis complete.")
