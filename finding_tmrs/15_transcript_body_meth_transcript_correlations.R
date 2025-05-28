# Calculate methylation-transcription correlations for CPGEA normal, CPGEA tumour and MCRPC samples

# Load required packages
library(methodical)
library(dplyr)

# Get TSS Granges
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")[names(tss_gr)]

# Expand transcripts_gr 5 KB upstream and downstream
transcripts_gr = methodical:::expand_granges(transcripts_gr, 5000, 5000)

# Load CPGEA methylation RSE and transcript counts
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/cpgea_normalized_kallisto_pcg_counts.tsv.gz"), row.names = 1)

# Mask CpGs with less than 10 reads covering them
assay(cpgea_meth_rse, 2)[is.na(assay(cpgea_meth_rse, 2))] = 0
assay(cpgea_meth_rse, 1)[assay(cpgea_meth_rse, 2) < 10] = NA

# Get CPGEA normal and tumour samples
normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)

# Create a bpparam object
bpparam = BiocParallel::SnowParam(workers = 5, type = "SOCK")

# Calculate methylation-transcription correlations for CPGEA normal samples. Took 34 minutes hours with 5 cores.  
system.time({transcript_meth_cors_cpgea_normal_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = normal_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_normal_samples_5kb, "meth_transcript_cors/cpgea_normal_whole_gene_body_correlations.rds")

# Calculate methylation-transcription correlations for CPGEA tumour samples. Took 42 minutes with 5 cores.  
system.time({transcript_meth_cors_cpgea_tumour_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = tumour_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_tumour_samples_5kb, "meth_transcript_cors/cpgea_tumour_whole_gene_body_correlations.rds")

### Calculate correlations for MCRPC
 
# Load MCRPC methylation RSE and transcript counts
mcrpc_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_transcript_counts.tsv.gz"), row.names = 1)

# Mask CpGs with less than 10 reads covering them
assay(mcrpc_meth_rse, 2)[is.na(assay(mcrpc_meth_rse, 2))] = 0
assay(mcrpc_meth_rse, 1)[assay(mcrpc_meth_rse, 2) < 10] = NA

# Get mcrpc normal and tumour samples
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Calculate methylation-transcription correlations for MCRPC samples. Took 29 hours with 5 cores.  
system.time({transcript_meth_cors_mcrpc_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = mcrpc_meth_rse, 
  transcript_expression_table = mcrpc_kallisto_deseq2_counts, samples_subset = common_mcrpc_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_mcrpc_samples_5kb, "meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")

### Repeat with 50 KB upstream and downstream

# Expand transcripts_gr_50kb 50 KB upstream and downstream
transcripts_gr_50kb = methodical:::expand_granges(transcripts_gr, 50000, 50000)

# Calculate methylation-transcription correlations for CPGEA normal samples. Took 45 minutes hours with 5 cores.  
system.time({transcript_meth_cors_cpgea_normal_samples_50kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = normal_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr_50kb, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_normal_samples_50kb, "meth_transcript_cors/cpgea_normal_whole_gene_body_correlations_50kb.rds")
rm(transcript_meth_cors_cpgea_normal_samples_50kb); gc()

# Calculate methylation-transcription correlations for CPGEA tumour samples. Took 42 minutes with 5 cores.  
system.time({transcript_meth_cors_cpgea_tumour_samples_50kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = tumour_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr_50kb, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_tumour_samples_50kb, "meth_transcript_cors/cpgea_tumour_whole_gene_body_correlations_50kb.rds")
rm(transcript_meth_cors_cpgea_tumour_samples_50kb); gc()

### Calculate correlations for MCRPC

# Load MCRPC methylation RSE and transcript counts
mcrpc_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_transcript_counts.tsv.gz"), row.names = 1)

# Get mcrpc normal and tumour samples
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Calculate methylation-transcription correlations for MCRPC samples. Took 29 hours with 5 cores.  
system.time({transcript_meth_cors_mcrpc_samples_50kb = calculateMethSiteTranscriptCors(meth_rse = mcrpc_meth_rse, 
  transcript_expression_table = mcrpc_kallisto_deseq2_counts, samples_subset = common_mcrpc_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr_50kb, 
  cor_method = "spearman", min_number_complete_pairs = 30, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_mcrpc_samples_50kb, "meth_transcript_cors/mcrpc_whole_gene_body_correlations_50kb.rds")
rm(transcript_meth_cors_cpgea_tumour_samples_50kb); gc()
""