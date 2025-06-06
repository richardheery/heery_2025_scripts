# Find TMRs in Roadmap samples

# Load required packages
library(methodical)
library(plotR)
library(patchwork)
source("../auxillary_scripts/tmr_plot_functions.R")

# Get CAGE-supported MANE transcripts
cage_supported_mane_transcripts = readRDS("../auxillary_data/cage_supported_mane_pc_transcript_ids.rds")

# Get TSS Granges for TSS and transcripts and subset for cage_supported_mane_transcripts
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")[cage_supported_mane_transcripts]
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")[cage_supported_mane_transcripts]

# Rename using gene ID
names(tss_gr) = gsub("\\.[0-9]*", "", tss_gr$gene_id)
names(transcripts_gr) = gsub("\\.[0-9]*", "", transcripts_gr$gene_id)

# Expand transcripts_gr
transcripts_gr = methodical:::expand_granges(transcripts_gr, 5000, 5000)

# Load Roadmap meth RSE and gene counts
roadmap_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/roadmap_meth_rse_hg38/")
gene_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/roadmap_gene_estimated_gene_counts.tsv.gz"), row.names = 1)

# Find common samples
common_samples = intersect(names(gene_counts), colnames(roadmap_meth_rse))

# Create a bpparm object
bpparam = BiocParallel::MulticoreParam(workers = 3) 

# Calculate methylation-transcription correlations. Took 22 minutes with 10 cores.  
system.time({transcript_meth_cors_roadmap = calculateMethSiteTranscriptCors(meth_rse = roadmap_meth_rse , 
  transcript_expression_table = gene_counts, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", samples_subset = common_samples, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_roadmap, "meth_transcript_cors/roadmap_whole_gene_body_correlations.rds")
transcript_meth_cors_roadmap = readRDS("meth_transcript_cors/roadmap_whole_gene_body_correlations.rds")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(5)

# Find TMRs for correlations. Took 1 minute. 
system.time({roadmap_tmrs = findTMRs(correlation_list = transcript_meth_cors_roadmap, 
  p_adjust_method = "fdr", p_value_threshold = 0.1, BPPARAM = bpparam)})

# Get repeat ranges for hg38 and subset TMRs for those not overlapping repeats
repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")
roadmap_tmrs = subsetByOverlaps(roadmap_tmrs, repeat_ranges, invert = T)
saveRDS(roadmap_tmrs, "tmr_granges/roadmap_tmrs.rds")