#!/usr/bin/env Rscript

# Colocalization of the OPG cis-pQTL and T2D GWAS signals using coloc.abf.

# -----------------------------------------------------------------------------
# 1. Settings and data
# -----------------------------------------------------------------------------

data_dir <- "data"
result_dir <- file.path("results", "colocalization")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

exposure <- read.csv(file.path(data_dir, "coloc_OPG_pQTL_chr8.csv"), check.names = FALSE)
outcome <- read.csv(file.path(data_dir, "coloc_T2D_GWAS_chr8.csv"), check.names = FALSE)

cis_window_bp <- 100000
n_exposure <- 33649
n_case <- 242283
n_control <- 1569734
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

required_exp <- c("gene_start", "gene_end", "CHROM", "POS38", "POS19", "rsid", "BETA", "SE", "LOG10P")
required_out <- c("Chromsome", "Position", "Beta", "SE", "Pval")
if (length(setdiff(required_exp, names(exposure))) > 0) stop("Required exposure columns are missing.")
if (length(setdiff(required_out, names(outcome))) > 0) stop("Required outcome columns are missing.")

# -----------------------------------------------------------------------------
# 2. Define cis region and align variants
# -----------------------------------------------------------------------------

gene_start <- unique(stats::na.omit(exposure$gene_start))
gene_end <- unique(stats::na.omit(exposure$gene_end))
if (length(gene_start) != 1 || length(gene_end) != 1) stop("Expected one OPG gene interval.")

# Cis region is defined in GRCh38; exposure-outcome matching uses GRCh37.
exposure <- exposure[
  exposure$POS38 >= gene_start - cis_window_bp &
    exposure$POS38 <= gene_end + cis_window_bp,
  ,
  drop = FALSE
]
exposure <- exposure[!duplicated(exposure$rsid), , drop = FALSE]
exposure$P <- 10^(-exposure$LOG10P)
exposure$chr_pos <- paste(exposure$CHROM, exposure$POS19, sep = "_")
outcome$chr_pos <- paste(outcome$Chromsome, outcome$Position, sep = "_")

variant_map <- exposure[, c("rsid", "chr_pos")]
outcome <- merge(variant_map, outcome, by = "chr_pos", all = FALSE)
outcome <- outcome[!duplicated(outcome$rsid), , drop = FALSE]

common <- intersect(exposure$rsid, outcome$rsid)
if (length(common) < 2) stop("Fewer than two overlapping variants were available.")

exposure <- exposure[match(common, exposure$rsid), , drop = FALSE]
outcome <- outcome[match(common, outcome$rsid), , drop = FALSE]
if (!identical(exposure$rsid, outcome$rsid)) stop("SNP alignment failed.")

write.csv(exposure, file.path(result_dir, "coloc_exposure_variants.csv"), row.names = FALSE)
write.csv(outcome, file.path(result_dir, "coloc_outcome_variants.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Colocalization
# -----------------------------------------------------------------------------

dataset_exposure <- list(
  type = "quant",
  beta = exposure$BETA,
  varbeta = exposure$SE^2,
  pvalues = exposure$P,
  snp = exposure$rsid,
  position = exposure$POS19,
  sdY = 1,
  N = n_exposure
)

dataset_outcome <- list(
  type = "cc",
  beta = outcome$Beta,
  varbeta = outcome$SE^2,
  pvalues = outcome$Pval,
  snp = outcome$rsid,
  position = outcome$Position,
  N = n_case + n_control,
  s = n_case / (n_case + n_control)
)

coloc::check_dataset(dataset_exposure)
coloc::check_dataset(dataset_outcome)

fit <- coloc::coloc.abf(
  dataset1 = dataset_exposure,
  dataset2 = dataset_outcome,
  p1 = p1, p2 = p2, p12 = p12
)

write.csv(
  data.frame(statistic = names(fit$summary), value = as.numeric(fit$summary)),
  file.path(result_dir, "coloc_abf_summary.csv"),
  row.names = FALSE
)
write.csv(fit$results, file.path(result_dir, "coloc_abf_variant_results.csv"), row.names = FALSE)
saveRDS(fit, file.path(result_dir, "coloc_abf_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
message("Colocalization analysis complete.")
