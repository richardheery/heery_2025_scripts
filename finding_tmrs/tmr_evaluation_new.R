# Evaluate the TMR correlations

# Load required packages
library(methodical)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")

# Load CPGEA methylation RSE and MCRPC RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
mcrpc_meth_rse =  HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")

# Get transcript counts for CPGEA and MCRPC
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_transcript_counts.tsv.gz"), row.names = 1)
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_transcript_counts.tsv.gz"), row.names = 1)

# Get CPGEA normal and tumour samples and MCRPC samples with both WGBS and RNA-seq data
common_cpgea_normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_cpgea_tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Get TMR GRanges
# cpgea_normal_tmrs = rtracklayer::import.bed("tmr_bed_files/normal_prostate_tmrs.bed")
# cpgea_normal_tmrs$ID = gsub("_.*", "", cpgea_normal_tmrs$tmr_name)
cpgea_normal_tmrs = readRDS("new_tmr_granges/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("new_tmr_granges/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("new_tmr_granges/mcrpc_tmrs.rds")

### Evaluate TMR correlations in normal prostate samples

#bpparam = BiocParallel::MulticoreParam(workers = 4)
bpparam = BiocParallel::SerialParam()

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BPPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_normal_samples, "tmr_evaluation_tables_new/cpgea_normal_tmrs_correlations_normal_samples.tsv.gz")

# # Filter for correlated and uncorrelated TMRs
# correlated_tmrs = filter(cpgea_normal_tmrs_correlations_normal_samples, q_val < 0.05)
# uncorrelated_tmrs = filter(cpgea_normal_tmrs_correlations_normal_samples, q_val > 0.05)
# boxplot(mean_tmr_coverage[uncorrelated_tmrs$genomic_region_name], mean_tmr_coverage[correlated_tmrs$genomic_region_name])
# 
# # Get the mean coverage for correlated and uncorrelated TMRs
# system.time({tmr_coverage = summarizeRegionMethylation(meth_rse = cpgea_meth_rse, assay = 2, genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, BPPARAM = bpparam)})
# system.time({tmr_meth = summarizeRegionMethylation(meth_rse = cpgea_meth_rse, assay = 1, genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, BPPARAM = bpparam)})
# mean_tmr_coverage = rowMeans(tmr_coverage, na.rm = T)
# boxplot(mean_tmr_coverage[uncorrelated_tmrs$genomic_region_name], mean_tmr_coverage[correlated_tmrs$genomic_region_name])

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_normal_samples, "tmr_evaluation_tables_new/cpgea_tumour_tmrs_correlations_normal_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_normal_samples, "tmr_evaluation_tables_new/mcrpc_tmrs_correlations_normal_samples.tsv.gz")

### Evaluate TMR correlations in prostate tumour samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_new/cpgea_normal_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_new/cpgea_tumour_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_new/mcrpc_tmrs_correlations_tumour_samples.tsv.gz")

### Evaluate TMR correlations in metastases samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_new/cpgea_normal_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_new/cpgea_tumour_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_new/mcrpc_tmrs_correlations_metastases_samples.tsv.gz")