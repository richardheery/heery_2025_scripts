### Identify TMRs in CPGEA normal, CPGEA tumour and MCRPC datasets

# Load required packages
library(methodical)
library(dplyr)
library(ggplot2)

# Get repeat ranges for hg38
unmappable_regions = readRDS("../auxillary_data/low_mappability_regions.rds")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(1)

# Load methylation-correlation results for CPGEA normal. Took 3 minutes. 
system.time({cpgea_normal_cor_results = readRDS("meth_transcript_cors/cpgea_normal_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA normal samples with even a single CpG sites. Took 20 minutes with 1 core.
system.time({cpgea_normal_tmrs_1_or_more_cpg = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, min_meth_sites = 1, BPPARAM = bpparam)})
cpgea_normal_tmrs_1_or_more_cpg = subsetByOverlaps(cpgea_normal_tmrs_1_or_more_cpg, unmappable_regions, invert = T)
saveRDS(cpgea_normal_tmrs_1_or_more_cpg, "tmr_granges/cpgea_normal_tmrs_1_or_more_cpg_unfiltered.rds")

# Find TMRs for CPGEA normal samples. Took 20 minutes with 1 core.
system.time({cpgea_normal_tmrs = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_normal_tmrs, "tmr_granges/cpgea_normal_tmrs_unfiltered.rds")
rm(cpgea_normal_cor_results); gc()

# Load methylation-correlation results for CPGEA tumour. Took 3 minutes. 
system.time({cpgea_tumour_cor_results = readRDS("meth_transcript_cors/cpgea_tumour_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA tumour samples. Took 25 minutes with 1 core.
system.time({cpgea_tumour_tmrs = findTMRs(correlation_list = cpgea_tumour_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_tumour_tmrs, "tmr_granges/cpgea_tumour_tmrs_unfiltered.rds")
rm(cpgea_tumour_cor_results); gc()

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")})

# Find TMRs for MCRPC samples. Took 17 minutes with 1 core.
system.time({mcrpc_tmrs = findTMRs(correlation_list = mcrpc_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(mcrpc_tmrs, "tmr_granges/mcrpc_tmrs_unfiltered.rds")
rm(mcrpc_cor_results); gc()

# Load all unfiltered TMRs within 5 kb
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs_unfiltered.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs_unfiltered.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs_unfiltered.rds")

# Remove TMRs overlapping unmappable regions and save
cpgea_normal_tmrs = subsetByOverlaps(cpgea_normal_tmrs, unmappable_regions, invert = T)
cpgea_tumour_tmrs = subsetByOverlaps(cpgea_tumour_tmrs, unmappable_regions, invert = T)
mcrpc_tmrs = subsetByOverlaps(mcrpc_tmrs, unmappable_regions, invert = T)
saveRDS(cpgea_normal_tmrs, "tmr_granges/cpgea_normal_tmrs.rds")
saveRDS(cpgea_tumour_tmrs, "tmr_granges/cpgea_tumour_tmrs.rds")
saveRDS(mcrpc_tmrs, "tmr_granges/mcrpc_tmrs.rds")

# Create a list of filtered TMRs and save
tmr_list = list(
  cpgea_normal_tmrs_negative = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Negative"],
  cpgea_normal_tmrs_positive = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Positive"],
  cpgea_tumour_tmrs_negative = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Negative"],
  cpgea_tumour_tmrs_positive = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Positive"],
  mcrpc_tmrs_negative = mcrpc_tmrs[mcrpc_tmrs$direction == "Negative"],
  mcrpc_tmrs_positive = mcrpc_tmrs[mcrpc_tmrs$direction == "Positive"]
)
saveRDS(tmr_list, "tmr_granges/tmr_list.rds")

### Find TMRs within 50KB of the transcript bodies

# Load methylation-correlation results for CPGEA normal. Took 1 minute. 
system.time({cpgea_normal_cor_results = readRDS("meth_transcript_cors/cpgea_normal_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for CPGEA normal samples. Took 20 minutes with 1 core.
system.time({cpgea_normal_tmrs_50kb = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_normal_tmrs_50kb, "tmr_granges/cpgea_normal_tmrs_50kb_unfiltered.rds")
rm(cpgea_normal_cor_results); gc()

# Load methylation-correlation results for CPGEA tumour. Took 3 minutes. 
system.time({cpgea_tumour_cor_results = readRDS("meth_transcript_cors/cpgea_tumour_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for CPGEA tumour samples. Took 25 minutes with 1 core.
system.time({cpgea_tumour_tmrs_50kb = findTMRs(correlation_list = cpgea_tumour_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(cpgea_tumour_tmrs_50kb, "tmr_granges/cpgea_tumour_tmrs_50kb_unfiltered.rds")
rm(cpgea_tumour_cor_results); gc()

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("meth_transcript_cors/mcrpc_whole_gene_body_correlations_50kb.rds")})

# Find TMRs for MCRPC samples. Took 17 minutes with 1 core.
system.time({mcrpc_tmrs_50kb = findTMRs(correlation_list = mcrpc_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
saveRDS(mcrpc_tmrs_50kb, "tmr_granges/mcrpc_tmrs_50kb_unfiltered.rds")
rm(mcrpc_cor_results); gc()

# Load all unfiltered TMRs within 50 kb
cpgea_normal_tmrs_50kb = readRDS("tmr_granges/cpgea_normal_tmrs_50kb_unfiltered.rds")
cpgea_tumour_tmrs_50kb = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb_unfiltered.rds")
mcrpc_tmrs_50kb = readRDS("tmr_granges/mcrpc_tmrs_50kb_unfiltered.rds")

# Remove TMRs overlapping unmappable regions and save
cpgea_normal_tmrs_50kb = subsetByOverlaps(cpgea_normal_tmrs_50kb, unmappable_regions, invert = T)
cpgea_tumour_tmrs_50kb = subsetByOverlaps(cpgea_tumour_tmrs_50kb, unmappable_regions, invert = T)
mcrpc_tmrs_50kb = subsetByOverlaps(mcrpc_tmrs_50kb, unmappable_regions, invert = T)
saveRDS(cpgea_normal_tmrs_50kb, "tmr_granges/cpgea_normal_tmrs_50kb.rds")
saveRDS(cpgea_tumour_tmrs_50kb, "tmr_granges/cpgea_tumour_tmrs_50kb.rds")
saveRDS(mcrpc_tmrs_50kb, "tmr_granges/mcrpc_tmrs_50kb.rds")