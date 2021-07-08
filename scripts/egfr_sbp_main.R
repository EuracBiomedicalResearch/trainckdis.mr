# eGFR and SBP main script ---------------------------------------------
library(TwoSampleMR)

# 1) Part one: all 256 index SNPs --------------------------------------
# 1.1) Load eGFR (exposure) --------------------------------------------
exp_data_dir <- paste(
  "data",
  "eGFR_EA_index_SNPs_ST11.txt",
  sep = "/"
)

exp_data_all <- read_exposure_data(
  filename = exp_data_dir,
  clump = TRUE,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  pval_col = "pvalue",
  samplesize_col = "n",
  gene_col = "gene",
  min_pval = 0,
  log_pval = FALSE
)

str(exp_data_all)
head(exp_data_all)

# 1.2) Load SBP (outcome) data -----------------------------------------
out_data_dir <- paste(
  "data",
  "Neale_SBP_eGFR_SNPs_res.txt",
  sep = "/"
)

out_data_all <- read_outcome_data(
  filename = out_data_dir,
  snps = exp_data_all$SNP,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "AF",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
  samplesize_col = "n_complete_samples",
  min_pval = 0,
  log_pval = FALSE
)

str(out_data_all)
head(out_data_all)

# 1.3) Harmonise datasets ----------------------------------------------
db <- harmonise_data(exp_data_all, out_data_all)

# 1.4) Apply Steiger filtering -----------------------------------------
steiger_data <- steiger_filtering(db)

mr_data <- db[steiger_data$steiger_dir, ]

# 1.5) Run MR analysis -------------------------------------------------
# Main function
mr_1 <- mr(mr_data)

mr_1

# Define MR methods
mr_method_list()

mr_2 <- mr(
  mr_data,
  method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression")
)

mr_2

# 1.6) Display results -------------------------------------------------
# Scatter plot
mr_scatter_plot(mr_2, mr_data)

# Forest plot
res_single <- mr_singlesnp(mr_data)
mr_forest_plot(res_single)

# 1.7) Sensitivity analyses --------------------------------------------
# Presence of heterogeneity
mr_heterogeneity(mr_data)

# Horizontal pleiotropy with MR-Egger intercept
mr_pleiotropy_test(mr_data)

# MR robust methods
mr_3 <- mr(
  mr_data,
  method_list = c(
    "mr_ivw", "mr_ivw_mre", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"
  )
)

mr_3

mr_scatter_plot(mr_3, mr_data)

# 2) Part two: all 256 index SNPs --------------------------------------
# 2.1) Load eGFR (exposure) --------------------------------------------
exp_data_dir <- paste(
  "data",
  "eGFR_EA_index_SNPs_ST11_BUN_filtered.txt",
  sep = "/"
)

exp_data_bun <- read_exposure_data(
  filename = exp_data_dir,
  clump = TRUE,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "effect_allele",
  other_allele_col = "non_effect_allele",
  pval_col = "pvalue",
  samplesize_col = "n",
  gene_col = "gene",
  min_pval = 0,
  log_pval = FALSE
)

str(exp_data_bun)
head(exp_data_bun)

# 2.2) Harmonise datasets ----------------------------------------------
db <- harmonise_data(exp_data_bun, out_data_all)

# 2.3) Apply Steiger filtering -----------------------------------------
steiger_data <- steiger_filtering(db)

mr_data <- db[steiger_data$steiger_dir, ]

# 2.4) Run MR analysis -------------------------------------------------
# Main function
mr_1 <- mr(mr_data)

mr_1

# Define MR methods
mr_method_list()

mr_2 <- mr(
  mr_data,
  method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression")
)

mr_2

# 2.5) Display results -------------------------------------------------
# Scatter plot
mr_scatter_plot(mr_2, mr_data)

# Forest plot
res_single <- mr_singlesnp(mr_data)
mr_forest_plot(res_single)

# 2.6) Sensitivity analyses --------------------------------------------
# Presence of heterogeneity
mr_heterogeneity(mr_data)

# Horizontal pleiotropy with MR-Egger intercept
mr_pleiotropy_test(mr_data)

# MR robust methods
mr_3 <- mr(
  mr_data,
  method_list = c(
    "mr_ivw", "mr_ivw_mre", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"
  )
)

mr_3

mr_scatter_plot(mr_3, mr_data)
