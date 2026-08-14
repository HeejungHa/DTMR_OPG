#!/usr/bin/env Rscript

# Prepare the primary OPG and T2D datasets used in the MR analyses.
# Run from the repository root.

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------

data_dir <- "data"

ukb_file <- file.path(data_dir, "UKB_100kb_0_1_only_for_SCALLOP.csv")
t2d_raw_file <- file.path(data_dir, "EUR_Metal_LDSC-CORR_Neff.v2.txt")

opg_gwas_id <- "ebi-a-GCST90012002"
target_chr <- 8
opg_gene_start <- 119935796
opg_gene_end <- 119964439
cis_window_bp <- 100000

p_threshold <- 5e-8
clump_r2 <- 0.1
clump_kb <- 10000

# -----------------------------------------------------------------------------
# 2. Select cis-pQTL instruments from SCALLOP
# -----------------------------------------------------------------------------

jwt <- ieugwasr::get_opengwas_jwt()
if (!nzchar(jwt)) stop("Set OPENGWAS_JWT in .Renviron before running this script.")

scallop <- TwoSampleMR::extract_instruments(
  outcomes = opg_gwas_id,
  p1 = p_threshold,
  clump = TRUE,
  r2 = clump_r2,
  kb = clump_kb,
  opengwas_jwt = jwt,
  force_server = FALSE
)

scallop_cis <- scallop[
  scallop$chr.exposure == target_chr &
    scallop$pos.exposure >= opg_gene_start - cis_window_bp &
    scallop$pos.exposure <= opg_gene_end + cis_window_bp,
  ,
  drop = FALSE
]

if (nrow(scallop_cis) == 0) stop("No SCALLOP instruments remained in the OPG cis region.")
message("SCALLOP cis instruments: ", nrow(scallop_cis))
write.csv(scallop_cis, file.path(data_dir, "SCALLOP_OPG_cis_instruments.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Prepare UKB-PPP exposure data
# -----------------------------------------------------------------------------

ukb <- read.csv(ukb_file, stringsAsFactors = FALSE, check.names = FALSE)
required_ukb <- c("rsid", "BETA", "SE", "ALLELE1", "ALLELE0", "A1FREQ", "LOG10P", "N")
missing_ukb <- setdiff(required_ukb, names(ukb))
if (length(missing_ukb) > 0) stop("Missing UKB-PPP column(s): ", paste(missing_ukb, collapse = ", "))

ukb <- ukb[!duplicated(ukb$rsid), , drop = FALSE]
idx <- match(scallop_cis$SNP, ukb$rsid)
primary_exposure <- ukb[idx[!is.na(idx)], , drop = FALSE]

if (nrow(primary_exposure) == 0) stop("No SCALLOP instruments matched the UKB-PPP data.")
message("Matched UKB-PPP instruments: ", nrow(primary_exposure))
write.csv(primary_exposure, file.path(data_dir, "primary_OPG_exposure.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Prepare T2DGGI outcome data
# -----------------------------------------------------------------------------

t2d <- data.table::fread(t2d_raw_file, data.table = FALSE)
chr_col <- intersect(c("Chromsome", "Chromosome"), names(t2d))[1]
if (is.na(chr_col)) stop("Could not identify the chromosome column in the T2DGGI file.")
if (!"Position" %in% names(t2d)) stop("Missing Position column in the T2DGGI file.")

required_t2d <- c("EffectAllele", "NonEffectAllele", "Beta", "SE", "EAF", "Pval", "Ncases", "Ncontrols", "Neff")
missing_t2d <- setdiff(required_t2d, names(t2d))
if (length(missing_t2d) > 0) stop("Missing T2DGGI column(s): ", paste(missing_t2d, collapse = ", "))

instrument_map <- scallop_cis[
  match(primary_exposure$rsid, scallop_cis$SNP),
  c("SNP", "pos.exposure"),
  drop = FALSE
]
names(instrument_map)[2] <- "Position"

outcome_chr <- t2d[t2d[[chr_col]] == target_chr, , drop = FALSE]
primary_outcome <- merge(instrument_map, outcome_chr, by = "Position", all = FALSE)

if (nrow(primary_outcome) == 0) stop("No OPG instruments matched the T2DGGI data by position.")
message("Matched T2DGGI instruments: ", nrow(primary_outcome))
if (chr_col != "Chromsome") names(primary_outcome)[names(primary_outcome) == chr_col] <- "Chromsome"

write.csv(primary_outcome, file.path(data_dir, "primary_T2D_outcome.csv"), row.names = FALSE)

message("Data processing complete.")
