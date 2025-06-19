
# Load required packages
library(methodical)
library(foreach)

# Get repeat ranges for hg38
unmappable_regions = readRDS("../auxillary_data/low_mappability_regions.rds")


#
offset_length = c(0, 5, 10, 20, 40, 60, 80, 100)
smoothing_factor = c(0, 0.1, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
parameters = setNames(expand.grid(offset_length, smoothing_factor), c("offset_length", "smoothing_factor"))

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")})

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(4)

system.time({tmr_list = lapply(seq.int(nrow(parameters)[1]), function(x)
  subsetByOverlaps(findTMRs(correlation_list = mcrpc_cor_results, p_adjust_method = "fdr", p_value_threshold = 0.05, 
    offset_length = parameters$offset_length[x], smoothing_factor = parameters$smoothing_factor[x], BPPARAM = bpparam), unmappable_regions))})


