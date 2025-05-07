### Identify TMRs within 50KB in CPGEA normal, CPGEA tumour and MCRPC datasets

# Load required packages
library(methodical)
library(dplyr)
library(ggplot2)

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(10)

# Load methylation-correlation results for CPGEA normal. Took 1 minute. 
system.time({cpgea_normal_cor_results = readRDS("cpgea_normal_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for CPGEA normal samples. Took 20 minutes with 1 core.
system.time({cpgea_normal_tmrs_50kb = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_normal_tmrs_50kb, "tmr_granges/cpgea_normal_tmrs_50kb.rds")
rm(cpgea_normal_cor_results); gc()

# Load methylation-correlation results for CPGEA tumour. Took 3 minutes. 
system.time({cpgea_tumour_cor_results = readRDS("cpgea_tumour_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for CPGEA tumour samples. Took 25 minutes with 1 core.
system.time({cpgea_tumour_tmrs_50kb = findTMRs(correlation_list = cpgea_tumour_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_tumour_tmrs_50kb, "tmr_granges/cpgea_tumour_tmrs_50kb.rds")
rm(cpgea_tumour_cor_results); gc()

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("mcrpc_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for MCRPC samples. Took 17 minutes with 1 core.
system.time({mcrpc_tmrs_50kb = findTMRs(correlation_list = mcrpc_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(mcrpc_tmrs_50kb, "tmr_granges/mcrpc_tmrs_50kb.rds")
rm(mcrpc_cor_results); gc()