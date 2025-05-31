### Identify TMRs in CPGEA normal, CPGEA tumour and MCRPC datasets

# Load required packages
library(methodical)
library(dplyr)
library(ggplot2)

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(1)

# Load methylation-correlation results for CPGEA normal. Took 3 minutes. 
system.time({cpgea_normal_cor_results = readRDS("meth_transcript_cors/cpgea_normal_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA normal samples. Took 20 minutes with 1 core.
system.time({cpgea_normal_tmrs = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_normal_tmrs, "tmr_granges/cpgea_normal_tmrs_with_repeats.rds")
rm(cpgea_normal_cor_results); gc()

# Load methylation-correlation results for CPGEA tumour. Took 3 minutes. 
system.time({cpgea_tumour_cor_results = readRDS("meth_transcript_cors/cpgea_tumour_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA tumour samples. Took 25 minutes with 1 core.
system.time({cpgea_tumour_tmrs = findTMRs(correlation_list = cpgea_tumour_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_tumour_tmrs, "tmr_granges/cpgea_tumour_tmrs_with_repeats.rds")
rm(cpgea_tumour_cor_results); gc()

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")})

# Find TMRs for MCRPC samples. Took 17 minutes with 1 core.
system.time({mcrpc_tmrs = findTMRs(correlation_list = mcrpc_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(mcrpc_tmrs, "tmr_granges/mcrpc_tmrs_with_repeats.rds")
rm(mcrpc_cor_results); gc()

tmr_list = list(
  cpgea_normal_tmrs_negative = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Negative"],
  cpgea_normal_tmrs_positive = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Positive"],
  cpgea_tumour_tmrs_negative = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Negative"],
  cpgea_tumour_tmrs_positive = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Positive"],
  mcrpc_tmrs_negative = mcrpc_tmrs[mcrpc_tmrs$direction == "Negative"],
  mcrpc_tmrs_positive = mcrpc_tmrs[mcrpc_tmrs$direction == "Positive"]
)

# Get repeat ranges for hg38
repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

# Save all TMRs without repeats
cpgea_normal_tmrs = subsetByOverlaps(cpgea_normal_tmrs, repeat_ranges, invert = T)
cpgea_tumour_tmrs = subsetByOverlaps(cpgea_tumour_tmrs, repeat_ranges, invert = T)
mcrpc_tmrs = subsetByOverlaps(mcrpc_tmrs, repeat_ranges, invert = T)

# Save CAGE-supported TMRs without repeats
saveRDS(cpgea_normal_tmrs, "tmr_granges/cpgea_normal_tmrs.rds")
saveRDS(cpgea_tumour_tmrs, "tmr_granges/cpgea_tumour_tmrs.rds")
saveRDS(mcrpc_tmrs, "tmr_granges/mcrpc_tmrs.rds")

# Create a list of TMRs which do not overlap repeats and save
tmr_list = list(
  cpgea_normal_tmrs_negative = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Negative"],
  cpgea_normal_tmrs_positive = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Positive"],
  cpgea_tumour_tmrs_negative = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Negative"],
  cpgea_tumour_tmrs_positive = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Positive"],
  mcrpc_tmrs_negative = mcrpc_tmrs[mcrpc_tmrs$direction == "Negative"],
  mcrpc_tmrs_positive = mcrpc_tmrs[mcrpc_tmrs$direction == "Positive"]
)
saveRDS(tmr_list, "tmr_granges/tmr_list.rds")
