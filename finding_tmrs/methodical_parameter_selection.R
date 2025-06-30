
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
bpparam = BiocParallel::MulticoreParam(20)

system.time({foreach(i = seq.int(nrow(parameters)), .packages = c("methodical", "GenomicRanges")) %do% {
  
  print(i)
  tmrs = findTMRs(correlation_list = mcrpc_cor_results, p_adjust_method = "fdr", p_value_threshold = 0.05, 
    offset_length = parameters$offset_length[i], smoothing_factor = parameters$smoothing_factor[i], BPPARAM = bpparam)
  
  tmrs = subsetByOverlaps(tmrs, unmappable_regions, invert = T)
  saveRDS(tmrs, paste0("methodical_parameter_tmrs/tmrs_", paste(unlist(parameters[i, ]), collapse = "_"), ".rds"))
  
}})

# Read in TMRs for all combinations of parameters
tmr_parameter_combination_list = list.files("methodical_parameter_tmrs", full.names = T)
names(tmr_parameter_combination_list) = gsub("tmrs_", "", tools::file_path_sans_ext(basename(tmr_parameter_combination_list)))
tmr_parameter_combination_list = lapply(tmr_parameter_combination_list, readRDS)

# Load CPGEA tumour WGBS and RNA-seq data
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_transcript_counts.tsv.gz"), row.names = 1)
common_cpgea_tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)

# Make a function which will evaluate TMRs in CPGEA tumour samples
evaluate_tmrs_in_cpgea_tumours = function(tmrs){
  
  results = calculateRegionMethylationTranscriptCors(
    genomic_regions = tmrs, genomic_region_names = tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
    genomic_region_transcripts = tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
    cor_method = "s", BPPARAM = bpparam)
  
  results
  
}

#bpparam = BiocParallel::SerialParam()
bpparam = BiocParallel::MulticoreParam(12)
system.time({tmr_parameter_combination_list_evaluation_results = lapply(tmr_parameter_combination_list, evaluate_tmrs_in_cpgea_tumours)})
saveRDS(tmr_parameter_combination_list_evaluation_results, "tmr_parameter_combination_list_evaluation_results.rds")

bpparam = BiocParallel::SerialParam()
bpparam = BiocParallel::MulticoreParam(10)
system.time({test = calculateRegionMethylationTranscriptCors(
  genomic_regions = x, genomic_region_names = x$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = x$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BPPARAM = bpparam)})
